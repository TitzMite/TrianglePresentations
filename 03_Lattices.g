
IsTorsionFreeTrianglePresentation := function(triangle_presentation)
    return ForAll(triangle_presentation, t-> not( t[1] = t[2] and t[2] = t[3]));
end;

TypeRotatingLattice := function(triangle_presentation)
    local n, points, f, gens, rels;
    points := List(triangle_presentation, t->t[1]);
    points := SSortedList(points);
    n := Size(points);
    f := FreeGroup(n);
    gens := GeneratorsOfGroup(f);
    rels := List(triangle_presentation, t-> gens[t[1]]*gens[t[2]]*gens[t[3]]);
    return f/rels;
end;

RewritingSystemTypeRotatingLattice := function(lattice)
    local rws, n, list, i, confluent;
    n := Size(GeneratorsOfGroup(lattice));
    rws := KBMAGRewritingSystem(lattice);
    list := [];
    for i in [1..n] do
        Append(list, [i, i+n]);
    od;
    ReorderAlphabetOfKBMAGRewritingSystem(rws, PermList(list));
    confluent := KnuthBendix(rws);
    if not confluent then
        Print("Knuth-Bendix Algorithm failed.\n");
    fi;
    return rws;
end;

CayleyBall := function(G, rws, R, A)
    local F, gens, vertices, N, edges, i, k, w, v, j, a, phi,
    imgs, permimgs, ballperms, Aball,gamma;
    F := FreeGroupOfFpGroup(G);
    gens := GeneratorsOfGroup(F);
    # Vertices = shortlex normal forms of length <= R
    vertices := EnumerateReducedWords(rws, 0, R);
    Sort(vertices);
    N := Length(vertices);
    ####################################################################
    # Induce A from the generators to the vertices of the ball
    ####################################################################
    ballperms := [];
    for a in GeneratorsOfGroup(A) do
        # Automorphism of the free group:
        #     x_i -> x_(i^a)
        imgs := List([1..Length(gens)], i -> gens[i^a]);
        phi := GroupHomomorphismByImages(F, F, gens, imgs);
        # Compute induced permutation on the normal forms in B_R
        permimgs := [];
        for i in [1..N] do
            w := Image(phi, vertices[i]);
            w := ReducedWord(rws, w);
            j := PositionSorted(vertices, w);
            if j > N or vertices[j] <> w then
                Error("A does not preserve the Cayley ball");
            fi;
            permimgs[i] := j;
        od;
        Add(ballperms, PermList(permimgs));
    od;
    Aball := Group(ballperms);
    ####################################################################
    # Edges
    ####################################################################
    edges := [];
    for i in [1..N] do
        w := vertices[i];
        for k in [1..Length(gens)] do
            v := ReducedWord(rws, w * gens[k]);
            if Length(v) <= R then
                j := PositionSorted(vertices, v);
                if j <= N and vertices[j] = v and i <> j then
                    Add(edges, [i, j]);
                fi;
            fi;
        od;
    od;
    ####################################################################
    # GRAPE graph with Aball stored as its group
    ####################################################################
    gamma := UnderlyingGraph(EdgeOrbitsGraph(Aball, edges, N));
    return gamma;
end;

AutGroupBall2 := function(gamma)
    local H, B1, delta, innerColours, A1, colours, outer;
    # Known automorphisms fixing the centre.
    H := Stabilizer(gamma.group, 1);
    # Vertices of B_1, in the original numbering.
    B1 := Concatenation([1], Adjacency(gamma, 1));
    # Retain the induced action of the known group.
    # Vertex i of delta corresponds to B1[i] in gamma.
    delta := InducedSubgraph(gamma, B1, H);
    # Compute the full automorphism group of B_1 fixing its centre.
    # Since B1 is sorted, its first vertex is the original vertex 1.
    innerColours := Filtered(
        [[1], Difference([1..Length(B1)], [1])],
        C -> Length(C) > 0
    );
    A1 := AutGroupGraph(rec(
        graph := delta,
        colourClasses := innerColours
    ));
    # Translate the A1-orbits back to gamma's vertex numbering.
    colours := List(
        Orbits(A1, [1..Length(B1)]),
        O -> Set(O, i -> B1[i])
    );
    # All vertices outside B_1 receive one additional colour.
    outer := Difference([1..gamma.order], B1);
    if Length(outer) > 0 then
        Add(colours, outer);
    fi;
    # Keep the known centre stabilizer associated with the graph.
    return AutGroupGraph(rec(
        graph := NewGroupGraph(H, gamma),
        colourClasses := colours
    ));
end;

AutGroupBallFixingInner := function(gamma, R, S)
    local ball, frontier, next, r, v, colours, outer, H;
    if not IsInt(R) or not IsInt(S) or S < 0 or R < S then
        Error("Require integers 0 <= S <= R");
    fi;
    # Find the S-ball around vertex 1.
    ball := [1];
    frontier := [1];
    for r in [1..S] do
        next := [];
        for v in frontier do
            UniteSet(next, Adjacency(gamma, v));
        od;
        frontier := Difference(next, ball);
        UniteSet(ball, frontier);
    od;
    # Fix every vertex of the S-ball individually.
    colours := List(ball, v -> [v]);
    outer := Difference([1..gamma.order], ball);
    if Length(outer) > 0 then
        Add(colours, outer);
    fi;
    # Retain the compatible part of the graph's known group.
    H := Stabilizer(gamma.group, ball, OnTuples);
    return AutGroupGraph(rec(
        graph := NewGroupGraph(H, gamma),
        colourClasses := colours
    ));
end;