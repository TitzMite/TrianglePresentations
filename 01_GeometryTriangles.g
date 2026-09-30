LoadPackage("GRAPE");
LoadPackage("KBMAG");

SingleAction := function(geo, tau)
    local geocopy, newgeo, n, lines;
    n := geo[1];
    lines := StructuralCopy(geo[2]);
    newgeo := [geo[1], List([1..Size(lines)], i -> lines[i^(tau^(-1))])];
    return newgeo;
end;

IsEmptyLineDatum := function(line_datum)
    return ForAll(line_datum, l -> l = []);
end;

#An arrow is an integer point [i,j], such that the point i
#is contained in the line j
ArrowsLineDatum := function(line_datum)
    local arrows, i, j, m;
    m := Size(line_datum);
    arrows := [];
    for j in [1..m] do
        for i in line_datum[j] do
            Add(arrows, [i,j]);
        od;
    od;
    arrows := SSortedList(arrows);
    return arrows;
end;

#this needs to be fast
#returns all triangles of line_datum, that contain e
TrianglesArrow := function(line_datum, e)
    local f, triangles, t;
    triangles := [];
    #if f is in line_datum[e[1]] then [f, e[1]] is an arrow
    for f in line_datum[e[1]] do
        #if e[2] is in line_datum[f], then [e[2], f] is an arrow
        if e[2] in line_datum[f] then
            if e[1] = f  and f = e[2] then
                t := [e];
                Add(triangles, t);
            else
                t := [[e[1],e[2]], [e[2],f], [f, e[1]]];
                Add(triangles, t);
            fi;
        fi;
    od;
    return triangles;
end;

SortedTriangle := function(t)
    local cyclics, rho;
    rho := (1,2,3);
    if Size(t) = 1 then
        return t;
    else
        cyclics := [t, Permuted(t, rho), Permuted(t, rho^2)];
        Sort(cyclics);
        return First(cyclics);
    fi;
end;

TrianglesLineDatum := function(line_datum)
    local triangles, arrows, e, local_triangles, t;
    triangles := [];
    arrows := ArrowsLineDatum(line_datum);
    for e in arrows do
        local_triangles := TrianglesArrow(line_datum, e);
        for t in local_triangles do
            if ForAll(triangles, t2 -> not IsEqualSet(t,t2)) then
                Add(triangles, t);
            fi;
        od;
    od;
    triangles := List(triangles, t->SortedTriangle(t));
    triangles := SSortedList(triangles);
    return triangles;
end;

TriangleGraphLineDatum := function(line_datum)
    local triangles, weights, k, edges_in_trianglegraph, trianglegraph;
    ###
    triangles := TrianglesLineDatum(line_datum);
    weights := List(triangles, t->Size(t));
    k := Size(triangles);
    ###
    edges_in_trianglegraph := Filtered(Combinations([1..k],2), x -> Intersection(triangles[x[1]],triangles[x[2]]) = []);
    Append(edges_in_trianglegraph, List(edges_in_trianglegraph, x->[x[2],x[1]]));
    ###
    trianglegraph := EdgeOrbitsGraph(Group(()), edges_in_trianglegraph, k);
    return [trianglegraph, weights, triangles];
end;

#returns all edge-disjoint triangle families witnessing the score
OptimalDisjointTriangleFamiliesLineDatum := function(line_datum)
    local arrows,
    trianglegraph_datum, trianglegraph, weights, triangles,
    score, maximal_covers,
    triangle_families;
    ###
    arrows := ArrowsLineDatum(line_datum);
    ###
    trianglegraph_datum := TriangleGraphLineDatum(line_datum);
    trianglegraph := trianglegraph_datum[1];
    weights := trianglegraph_datum[2];
    triangles := trianglegraph_datum[3];
    ###
    score := Size(arrows);
    ###
    maximal_covers := CompleteSubgraphsOfGivenSize(trianglegraph, score, 2, false, false, weights);
    while maximal_covers = [] do
        score := score - 1;
        maximal_covers := CompleteSubgraphsOfGivenSize(trianglegraph, score, 2, false, false, weights);
    od;
    triangle_families := List(maximal_covers, cov->List(cov, i->triangles[i]));
    return triangle_families;
end;

OptimalScore := function(geo)
    local optimal_score;
    optimal_score := Size(ArrowsLineDatum(geo[2]));
    return optimal_score;
end;

OptimalDisjointTriangleFamilies := function(geo)
    local linedatum;
    linedatum := StructuralCopy(geo[2]);
    return OptimalDisjointTriangleFamiliesLineDatum(linedatum);
end;

ScoreTriangleFamily := function(triangles)
    local size, t;
    size := 0;
    for t in triangles do
        size :=  size + Size(t);
    od;
    return size;
end;

#next function only check if the geo is perfect
HardPruneLineDatumToCheckPerfectness := function(line_datum)
    local
    deleted_triangles, changed, m, perm,
    i, j, k, points, triangles, a, t, potentially_perfect;
    ###
    deleted_triangles := [];
    changed := false;
    potentially_perfect := true;
    m := Size(line_datum);
    perm := Random(SymmetricGroup(m));
    for k in [1..m] do
        j := k^perm;
        points := StructuralCopy(line_datum[j]);
        while not points = [] do
            i := Random(points);
            RemoveSet(points, i);
            triangles := TrianglesArrow(line_datum, [i,j]);
            if Size(triangles) = 0 then
                RemoveSet(line_datum[j], i);
                changed := true;
                potentially_perfect := false;
                break;
            elif Size(triangles) = 1 then
                t := triangles[1];
                for a in t do
                    if [i,j] = a then
                        RemoveSet(line_datum[j], i);
                        changed := true;
                    else
                        if a[2] = j then
                            RemoveSet(points, a[1]);
                        fi;
                        RemoveSet(line_datum[a[2]], a[1]);
                    fi;
                od;
                Add(deleted_triangles, t);
            fi;
        od;
        if not potentially_perfect then
            break;
        fi;
    od;
    return [deleted_triangles, changed, potentially_perfect];
end;

IsPerfectGeo := function(geo)
    local line_datum, triangles, arrows, prune_outcome, more_triangles,
    potentially_perfect;
    line_datum :=StructuralCopy(geo[2]);
    triangles := [];
    potentially_perfect := true;
    repeat
        prune_outcome := HardPruneLineDatumToCheckPerfectness(line_datum);
        if prune_outcome[3] = false then
            potentially_perfect := false;
            break;
        fi;
        Append(triangles, prune_outcome[1]);
    until prune_outcome[2] = false;
    if not potentially_perfect then
        return false;
    else
        if not IsEmptyLineDatum(line_datum) then
            Print("Calculating this score is more expensive...\n");
            more_triangles := OptimalDisjointTriangleFamiliesLineDatum(line_datum)[1];
            Append(triangles, more_triangles);
        fi;
        return OptimalScore(geo) = ScoreTriangleFamily(triangles);
    fi;
end;