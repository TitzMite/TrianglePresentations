# A directed multigraph is stored as
#
#     [n, edges]
#
# where edges is a list of ordered pairs [i,j].
# Parallel edges are allowed.
# The identity of an edge is its POSITION in edges.

TrianglesMultiGraph := function(multigraph)
    local edges, m, triangles, e, abc,
          a, b, c, ea, eb, ec;
    edges := multigraph[2];
    m := Length(edges);
    triangles := [];
    # Singleton loops, following the convention used for
    # triangle presentations [i,i,i].
    for e in [1..m] do
        if edges[e][1] = edges[e][2] then
            Add(triangles, [e]);
        fi;
    od;
    # Triangles made from three distinct edge copies.
    for abc in Combinations([1..m], 3) do
        a := abc[1];
        b := abc[2];
        c := abc[3];
        ea := edges[a];
        eb := edges[b];
        ec := edges[c];
        # With edge a fixed first, there are only two
        # possible cyclic orders.
        if (ea[2] = eb[1] and eb[2] = ec[1] and ec[2] = ea[1])
            or
            (ea[2] = ec[1] and ec[2] = eb[1] and eb[2] = ea[1])
        then
            Add(triangles, abc);
        fi;
    od;
    return triangles;
end;


# One pruning pass.
#
# active_edges     = edge IDs still uncovered
# active_triangles = triangle IDs still available
#
# Returns
#
# [
#     new_active_edges,
#     new_active_triangles,
#     forced_triangles,
#     changed,
#     potentially_perfect
# ]

PruneMultiGraph := function(triangles,active_edges,active_triangles)
    local edges, forced_triangles, changed, e, possible, t, tri;
    edges := ShallowCopy(active_edges);
    forced_triangles := [];
    changed := false;
    for e in active_edges do
        # e may already have been covered by a triangle
        # forced earlier in this pass.
        if e in edges then
            possible := Filtered(active_triangles, t -> e in triangles[t]);
            # No available triangle contains e.
            if Length(possible) = 0 then
                return [edges,active_triangles,forced_triangles,changed,false];
            fi;
            # Exactly one triangle contains e, so it is forced.
            if Length(possible) = 1 then
                t := possible[1];
                tri := triangles[t];
                Add(forced_triangles, tri);
                # The edges of this triangle are now covered.
                edges := Difference(edges, tri);
                # Any triangle meeting one of these edges
                # can no longer be used.
                active_triangles := Filtered(
                    active_triangles,
                    x -> Intersection(triangles[x], tri) = []
                );
                changed := true;
            fi;
        fi;
    od;
    return [edges,active_triangles,forced_triangles,changed,true];
end;


# Repeatedly prune until:
#
#   - impossibility is detected, or
#   - every edge is covered, or
#   - no more triangles are forced.
#
# Returns
#
# [remaining_edges,remaining_triangles,forced_triangles,potentially_perfect]
#

HardPruneMultiGraph := function(multigraph)
    local triangles,active_edges,active_triangles,forced_triangles,outcome;
    triangles := TrianglesMultiGraph(multigraph);
    active_edges := [1..Length(multigraph[2])];
    active_triangles := [1..Length(triangles)];
    forced_triangles := [];
    repeat
        outcome := PruneMultiGraph(triangles,active_edges,active_triangles);
        active_edges := outcome[1];
        active_triangles := outcome[2];
        Append(forced_triangles, outcome[3]);
        if not outcome[5] then
            return [
                active_edges,
                active_triangles,
                forced_triangles,
                false
            ];
        fi;
    until not outcome[4];
    return [active_edges,active_triangles,forced_triangles,true];
end;

CanCompletePerfectTriangleCoverMultiGraph := function(multigraph, outcome)
    local triangles, active_edges, active_triangles, residual_triangles, weights, k,
          edges_in_trianglegraph, trianglegraph, covers, score;
    # Hard pruning already proved impossibility
    if not outcome[4] then
        return false;
    fi;
    active_edges := outcome[1];
    # Hard pruning already found a perfect cover
    if active_edges = [] then
        return true;
    fi;
    triangles := TrianglesMultiGraph(multigraph);
    active_triangles := outcome[2];
    residual_triangles := triangles{active_triangles};
    weights := List(residual_triangles, t -> Size(t));
    k := Size(residual_triangles);
    # Two triangle vertices are adjacent iff the triangles
    # are edge-disjoint.
    edges_in_trianglegraph :=
        Filtered(
            Combinations([1..k], 2),
            x -> Intersection(residual_triangles[x[1]],residual_triangles[x[2]]) = []
        );
    Append(edges_in_trianglegraph,List(edges_in_trianglegraph, x -> [x[2],x[1]]));
    trianglegraph := EdgeOrbitsGraph(Group(()),edges_in_trianglegraph,k);
    # We need to cover exactly all remaining edges.
    score := Size(active_edges);
    covers :=
        CompleteSubgraphsOfGivenSize(
            trianglegraph,
            score,
            2,
            false,
            false,
            weights
        );
    return covers <> [];
end;

IsPerfectMultiGraph := function(multigraph)
    local prune_outcome;
    prune_outcome := HardPruneMultiGraph(multigraph);
    # Pruning proved that no perfect triangle cover exists.
    if not prune_outcome[4] then
        return false;
    fi;
    # Pruning already covered every edge.
    if prune_outcome[1] = [] then
        return true;
    fi;
    # A nontrivial core remains; use the expensive exact test.
    return CanCompletePerfectTriangleCoverMultiGraph(
        multigraph,
        prune_outcome
    );
end;

IsPerfectMultiGraphFast := function(multigraph)
    local triangles,recurse;
    triangles := TrianglesMultiGraph(multigraph);
    recurse := function(active_edges,active_triangles)
        local outcome,e,possible,counts,mincount,t,tri;
        repeat
            outcome := PruneMultiGraph(
                triangles,
                active_edges,
                active_triangles
            );
            active_edges := outcome[1];
            active_triangles := outcome[2];
            if not outcome[5] then
                return false;
            fi;
        until not outcome[4];
        if active_edges = [] then
            return true;
        fi;
        counts := List(
            active_edges,
            e -> Length(
                Filtered(
                    active_triangles,
                    t -> e in triangles[t]
                )
            )
        );
        mincount := Minimum(counts);
        e := active_edges[Position(counts,mincount)];
        possible := Filtered(
            active_triangles,
            t -> e in triangles[t]
        );
        for t in possible do
            tri := triangles[t];
            if recurse(
                Difference(active_edges,tri),
                Filtered(
                    active_triangles,
                    x -> Intersection(triangles[x],tri) = []
                )
            ) then
                return true;
            fi;
        od;
        return false;
    end;
    return recurse(
        [1..Length(multigraph[2])],
        [1..Length(triangles)]
    );
end;

PerfectTriangleCoversMultiGraph := function(multigraph)
    local outcome, triangles, active_edges, active_triangles,
          forced_triangles, residual_triangles, weights, k,
          edges_in_trianglegraph, trianglegraph, covers, score;
    outcome := HardPruneMultiGraph(multigraph);
    if not outcome[4] then
        return [];
    fi;
    active_edges := outcome[1];
    active_triangles := outcome[2];
    forced_triangles := outcome[3];
    if active_edges = [] then
        return [forced_triangles];
    fi;
    triangles := TrianglesMultiGraph(multigraph);
    residual_triangles := triangles{active_triangles};
    weights := List(residual_triangles, t -> Size(t));
    k := Size(residual_triangles);
    edges_in_trianglegraph :=
        Filtered(
            Combinations([1..k],2),
            x -> Intersection(
                residual_triangles[x[1]],
                residual_triangles[x[2]]
            ) = []
        );
    Append(
        edges_in_trianglegraph,
        List(edges_in_trianglegraph, x -> [x[2],x[1]])
    );
    trianglegraph := EdgeOrbitsGraph(
        Group(()),
        edges_in_trianglegraph,
        k
    );
    score := Size(active_edges);
    covers := CompleteSubgraphsOfGivenSize(
        trianglegraph,
        score,
        2,
        false,
        false,
        weights
    );
    return List(
        covers,
        cov -> Concatenation(
            forced_triangles,
            List(cov, i -> residual_triangles[i])
        )
    );
end;

PerfectTriangleCoversMultiGraphFast := function(multigraph)
    local triangles,covers,recurse;
    triangles := TrianglesMultiGraph(multigraph);
    covers := [];
    recurse := function(active_edges,active_triangles,chosen)
        local outcome,e,possible,counts,mincount,t,tri;
        repeat
            outcome := PruneMultiGraph(
                triangles,
                active_edges,
                active_triangles
            );
            active_edges := outcome[1];
            active_triangles := outcome[2];
            Append(chosen,outcome[3]);
            if not outcome[5] then
                return;
            fi;
        until not outcome[4];
        if active_edges = [] then
            Add(covers,ShallowCopy(chosen));
            return;
        fi;
        counts := List(
            active_edges,
            e -> Length(
                Filtered(
                    active_triangles,
                    t -> e in triangles[t]
                )
            )
        );
        mincount := Minimum(counts);
        e := active_edges[Position(counts,mincount)];
        possible := Filtered(
            active_triangles,
            t -> e in triangles[t]
        );
        for t in possible do
            tri := triangles[t];
            recurse(
                Difference(active_edges,tri),
                Filtered(
                    active_triangles,
                    x -> Intersection(triangles[x],tri) = []
                ),
                Concatenation(chosen,[tri])
            );
        od;
    end;
    recurse(
        [1..Length(multigraph[2])],
        [1..Length(triangles)],
        []
    );
    return covers;
end;


SingleActionMultiGraph := function(multigraph, tau)
    local n, edges;
    n := multigraph[1];
    edges := multigraph[2];
    return [n,List(edges, e -> [e[1], e[2]^tau])];
end;
