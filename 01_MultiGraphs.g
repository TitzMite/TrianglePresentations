LoadPackage("GRAPE");
LoadPackage("KBMAG");

TrianglesMultiGraph := function(multigraph)
    local edges, outgoing, triangles, e, a, b, c;
    edges := multigraph[2];
    outgoing := List([1..multigraph[1]], i -> []);
    triangles := [];
    # Index outgoing edges and record singleton triangles.
    for e in [1..Length(edges)] do
        Add(outgoing[edges[e][1]], e);
        if edges[e][1] = edges[e][2] then
            Add(triangles, [e]);
        fi;
    od;
    # Follow a: u -> v, b: v -> w, c: w -> u.
    for a in [1..Length(edges)] do
        for b in outgoing[edges[a][2]] do
            if b > a then
                for c in outgoing[edges[b][2]] do
                    if c > a and c <> b and
                       edges[c][2] = edges[a][1] then
                        Add(triangles, [a,b,c]);
                    fi;
                od;
            fi;
        od;
    od;
    return triangles;
end;

# 1. Enumerates all triangles of the multigraph, treating parallel edges as distinct.
# 2. Applies light pruning: removes edges in no triangle and selects triangles with at least two uniquely contained edges.
# 3. Applies hard pruning: selects triangles forced by any uniquely contained edge.
# 4. Repeats each pruning phase until no further change occurs.
# 5. Selecting a triangle counts its edges and excludes all triangles overlapping it.
# 6. Computes the exact maximum additional score on the residual multigraph by backtracking.
# 7. Branches over triangles covering a chosen edge, and over leaving that edge uncovered.
# 8. Returns the total number of covered edges without modifying the input.
# 9. The score equals the total edge count exactly when a perfect cover exists.
# 10. Otherwise, it is a lower bound that may depend on the randomized pruning order.
ScoreMultiGraph := function(multigraph)
    local triangles, active_edges, active_triangles, score,
          light, changed, pending, e, possible, t, tri,
          select_triangle, best, search;
    #   
    triangles := TrianglesMultiGraph(multigraph);
    active_edges := [1..Length(multigraph[2])];
    active_triangles := [1..Length(triangles)];
    score := 0;
    # Select a triangle during pruning.
    select_triangle := function(t)
        local tri;
        tri := triangles[t];
        score := score + Length(tri);
        active_edges := Difference(active_edges, tri);
        active_triangles := Filtered(
            active_triangles,
            u -> Intersection(triangles[u], tri) = []);
    end;
    # First light pruning, then hard pruning.
    for light in [true, false] do
        repeat
            changed := false;
            pending := ShallowCopy(active_edges);
            while pending <> [] do
                e := Random(pending);
                RemoveSet(pending, e);
                if e in active_edges then
                    possible := Filtered(
                        active_triangles,
                        t -> e in triangles[t]);
                    if Length(possible) = 0 then
                        RemoveSet(active_edges, e);
                        changed := true;
                    elif Length(possible) = 1 then
                        t := possible[1];
                        tri := triangles[t];
                        # Light pruning requires a second edge
                        # belonging to this triangle only.
                        if not light or ForAny(
                            Filtered(
                                tri, f -> f <> e),
                                f ->
                                    Number(active_triangles,
                                    u -> f in triangles[u]) = 1
                            ) then
                            select_triangle(t);
                            changed := true;
                        fi;
                    fi;
                fi;
            od;
        until not changed;
    od;
    if active_edges = [] then
        return score;
    fi;
    # Exact maximum-score search on the residual problem.
    best := score;
    search := function(edges, available, current)
        local counts, usable, i, minimum, e, possible, t, tri, remaining, upper;
        if current > best then
            best := current;
        fi;
        if available = [] then
            return;
        fi;
        # Discard edges that no remaining triangle can cover.
        counts := List(edges,e -> Number(available, t -> e in triangles[t]));
        usable := Filtered([1..Length(edges)], i -> counts[i] > 0);
        edges := edges{usable};
        counts := counts{usable};
        # Even covering every remaining edge cannot improve best.
        upper := current + Length(edges);
        if upper <= best then
            return;
        fi;
        # Branch on an edge with the fewest available triangles.
        minimum := Minimum(counts);
        e := edges[Position(counts, minimum)];
        possible := Filtered(available, t -> e in triangles[t]);
        for t in possible do
            tri := triangles[t];
            remaining := Filtered(
                available,
                u -> Intersection(triangles[u], tri) = []
            );
            search(
                Difference(edges, tri),
                remaining,
                current + Length(tri)
            );
            if best = upper then
                return;
            fi;
        od;
        # A maximum partial cover may leave e uncovered.
        search(
            Difference(edges, [e]),
            Filtered(available, t -> not (e in triangles[t])),
            current);
    end;
    search(active_edges, active_triangles, score);
    return best;
end;

# One pruning pass.
# Returns:
# [remaining_edges, remaining_triangle_IDs, forced_triangles,
#  changed, potentially_perfect]
PruneMultiGraph := function(triangles, active_edges, active_triangles)
    local edges, forced_triangles, changed, e, possible, t, tri;
    edges := ShallowCopy(active_edges);
    forced_triangles := [];
    changed := false;
    for e in active_edges do
        if e in edges then
            possible := Filtered(
                active_triangles,
                t -> e in triangles[t]
            );
            if Length(possible) = 0 then
                return [edges, active_triangles, forced_triangles,
                    changed, false];
            fi;
            if Length(possible) = 1 then
                t := possible[1];
                tri := triangles[t];
                Add(forced_triangles, tri);
                edges := Difference(edges, tri);
                active_triangles := Filtered(
                    active_triangles,
                    u -> Intersection(triangles[u], tri) = []
                );
                changed := true;
            fi;
        fi;
    od;
    return [edges, active_triangles, forced_triangles,changed, true];
end;

# 1. Enumerates the multigraph's triangles and repeatedly selects forced triangles.
# 2. Rejects a branch if an uncovered edge belongs to no available triangle.
# 3. Otherwise, chooses an edge with the fewest available triangles and branches over those choices.
# 4. Each choice removes the covered edges and excludes overlapping triangles.
# 5. Records a cover whenever all edges are covered, returning all covers in your triangle representation.
AllPerfectCovers := function(multigraph)
    local triangles, covers, search;
    triangles := TrianglesMultiGraph(multigraph);
    covers := [];
    search := function(active_edges, active_triangles, chosen)
        local outcome, counts, e, possible, t, tri;
        # Select forced triangles.
        repeat
            outcome := PruneMultiGraph(triangles, active_edges, active_triangles);
            active_edges := outcome[1];
            active_triangles := outcome[2];
            Append(chosen, outcome[3]);
            if not outcome[5] then
                return;
            fi;
        until not outcome[4];
        if active_edges = [] then
            Add(covers, ShallowCopy(chosen));
            return;
        fi;
        # Branch on an edge with the fewest available triangles.
        counts := List(active_edges,e -> Number(active_triangles, t -> e in triangles[t]));
        e := active_edges[Position(counts, Minimum(counts))];
        possible := Filtered(active_triangles, t -> e in triangles[t]);
        for t in possible do
            tri := triangles[t];
            search(
                Difference(active_edges, tri),
                Filtered(active_triangles, u -> Intersection(triangles[u], tri) = []),
                Concatenation(chosen, [tri]));
        od;
    end;
    search([1..Length(multigraph[2])],[1..Length(triangles)],[]);
    return covers;
end;