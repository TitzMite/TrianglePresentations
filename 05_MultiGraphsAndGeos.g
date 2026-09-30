
MatchedGeometryData := function(geo, pre_image_function, tau, multigraph)
    local g, newgeo, newmultigraph;
    # Lift the permutation of A-orbits to the original point set.
    g := pre_image_function(tau);
    # Apply the corresponding line matching upstairs.
    newgeo := SingleAction(geo, g);
    # Apply the same orbit matching to the quotient multigraph.
    newmultigraph := SingleActionMultiGraph(multigraph, tau);
    return [newgeo, newmultigraph];
end;

#computes lifts
LiftMultiGraphEdges := function(geo,A,multigraph)
    local orbits, inc_orbits, edge_orbits, used, e, candidates, o;
    orbits := Orbits(A,[1..geo[1]]);
    inc_orbits := Orbits(
        A,
        Concatenation(
            List([1..geo[1]], l -> List(geo[2][l], p -> [p,l]))
        ),
        OnTuples
    );
    edge_orbits := [];
    used := [];
    for e in multigraph[2] do
        candidates := Filtered([1..Length(inc_orbits)], i ->
            not i in used and
            inc_orbits[i][1][1] in orbits[e[1]] and
            inc_orbits[i][1][2] in orbits[e[2]]
        );
        o := candidates[1];
        Add(used,o);
        Add(edge_orbits,inc_orbits[o]);
    od;
    return [orbits,edge_orbits];
end;

# only use for cyclic groups
LiftEquationSystem := function(geo,A,multigraph,cover)
    local data, orbits, edge_orbits, els, offsets, e, p, l, a, b, equations, t, row, int_grp_elm, gen, n;
    n := Size(A);
    gen := MinimalGeneratingSet(A)[1];
    int_grp_elm := function(a)
        return First([0..n-1], i -> gen^i = a);
    end;
    data := LiftMultiGraphEdges(geo,A,multigraph);
    orbits := data[1];
    edge_orbits := data[2];
    els := Elements(A);
    offsets := [];
    for e in [1..Length(edge_orbits)] do
        p := edge_orbits[e][1][1];
        l := edge_orbits[e][1][2];
        a := First(
            els,
            x -> Minimum(orbits[multigraph[2][e][1]])^x = p
        );
        b := First(
            els,
            x -> Minimum(orbits[multigraph[2][e][2]])^x = l
        );
        Add(offsets,int_grp_elm(a^-1*b));
    od;
    equations := [];
    for t in cover do
        row := List([1..multigraph[1]],i -> 0);
        for e in t do
            row[multigraph[2][e][2]] :=
                row[multigraph[2][e][2]] + 1;
        od;
        Add(equations,[row,Sum(List(t,e -> offsets[e])) mod n]);
    od;
    return equations;
end;

SolutionsModN := function(M,b,n)
    local snf, D, P, Q, rank, m, r, c, options,
          i, d, g, d1, c1, n1, base, ys;
    m := Length(M);
    r := Length(M[1]);
    snf := SmithNormalFormIntegerMatTransforms(M);
    D := snf.normal;
    P := snf.rowtrans;
    Q := snf.coltrans;
    rank := snf.rank;
    c := List(P, row -> ScalarProduct(row,b));
    c := List(c, x -> x mod n);
    for i in [rank+1..m] do
        if c[i] <> 0 then
            return [];
        fi;
    od;
    options := [];
    for i in [1..rank] do
        d := D[i][i];
        g := GcdInt(d,n);
        if c[i] mod g <> 0 then
            return [];
        fi;
        d1 := d/g;
        c1 := c[i]/g;
        n1 := n/g;
        if n1 = 1 then
            base := 0;
        else
            base := (c1*PowerModInt(d1,-1,n1)) mod n1;
        fi;
        Add(options,List([0..g-1], k -> (base+k*n1) mod n));
    od;
    for i in [rank+1..r] do
        Add(options,[0..n-1]);
    od;
    ys := [[]];
    for i in [1..r] do
        ys := Concatenation(
            List(ys, y -> List(options[i], a -> Concatenation(y,[a])))
        );
    od;
    return List(
        ys,
        y -> List(Q, row -> ScalarProduct(row,y) mod n)
    );
end;

TrianglePresentationFromSolution := function(geo,A,multigraph,cover,solution)
    local gen, data, edge_orbits, shifted_orbits, e, t, ord, a, b, c, triangle_family, triangle_presentation;
    gen := MinimalGeneratingSet(A)[1];
    data := LiftMultiGraphEdges(geo,A,multigraph);
    edge_orbits := data[2];
    shifted_orbits := List([1..Length(edge_orbits)], e ->
        List(edge_orbits[e], a ->
            [
                a[1],
                a[2]^(gen^(-solution[multigraph[2][e][2]]))
            ]
        )
    );
    triangle_family := [];
    for t in cover do
        if Length(t) = 1 then
            for a in shifted_orbits[t[1]] do
                Add(triangle_family,[a]);
            od;
        else
            ord := First(PermutationsList(t), q ->
                multigraph[2][q[1]][2] = multigraph[2][q[2]][1] and
                multigraph[2][q[2]][2] = multigraph[2][q[3]][1] and
                multigraph[2][q[3]][2] = multigraph[2][q[1]][1]
            );
            for a in shifted_orbits[ord[1]] do
                b := First(
                    shifted_orbits[ord[2]],
                    x -> x[1] = a[2]
                );
                c := First(
                    shifted_orbits[ord[3]],
                    x -> x[1] = b[2] and x[2] = a[1]
                );
                Add(triangle_family,[a,b,c]);
            od;
        fi;
    od;
    triangle_presentation := [];
    for t in triangle_family do
        Append(
            triangle_presentation,
            PresentationFormsTriangle(t)
        );
    od;
    return SSortedList(triangle_presentation);
end;
