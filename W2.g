
Read("01_GeometryTriangles.g");
Read("02_TrianglePresentations.g");
Read("03_Lattices.g");
Read("04_MultiGraphs.g");
Read("05_MultiGraphsAndGeos.g");

W2 := [ 15,
    [ [ 1, 6, 14 ], [ 2, 7, 15 ], [ 3, 8, 11 ],
      [ 4, 9, 12 ], [ 5, 10, 13 ], [ 2, 3, 10 ],
      [ 3, 4, 6 ], [ 4, 5, 7 ], [ 1, 5, 8 ],
      [ 1, 2, 9 ], [ 6, 13, 15 ], [ 7, 11, 14 ],
      [ 8, 12, 15 ], [ 9, 11, 13 ], [ 10, 12, 14 ] ] ];

# C5 action on W2
s :=
    (1,2,3,4,5)
    (6,7,8,9,10)
    (11,12,13,14,15)
    (16,17,18,19,20)
    (21,22,23,24,25)
    (26,27,28,29,30);

A := Group(s);

# centralizer of <s> restricted to the 15 points:
# C5 wr Sym(3)
C := Group(
    # C5 acting on the first orbit
    (1,2,3,4,5),
    # generators of Sym3 acting on the three C5-orbits
    (1,6)(2,7)(3,8)(4,9)(5,10),
    (1,6,11)
    (2,7,12)
    (3,8,13)
    (4,9,14)
    (5,10,15)
);

Sym3 := SymmetricGroup(3);

gensC := GeneratorsOfGroup(C);

phi := GroupHomomorphismByImages(
    C,
    Sym3,
    gensC,
    [(),(1,2),(1,2,3)]
);

PreimageInC := function(tau)
    return PreImagesRepresentative(phi,tau);
end;

W2MultiGraph :=
    [ 3,
    [   [1,1], [1,2], [1,2],
        [2,1], [2,2], [2,3],
        [3,1], [3,3], [3,3] ] ];

compute_W2_tps := function()
    local tau, multi_graph_covers, new_multi_graph, tps, new_geo,
          sols, equation_data, cover, M, b, solution;
    tps := [];
    for tau in Sym3 do
        new_multi_graph :=
            SingleActionMultiGraph(W2MultiGraph,tau);
        if IsPerfectMultiGraph(new_multi_graph) then
            new_geo := SingleAction(W2,PreimageInC(tau));
            multi_graph_covers := PerfectTriangleCoversMultiGraph(new_multi_graph);
            for cover in multi_graph_covers do
                equation_data :=
                    LiftEquationSystem(
                        new_geo,
                        A,
                        new_multi_graph,
                        cover
                    );
                M := List(equation_data,p -> p[1]);
                b := List(equation_data,p -> p[2]);
                sols := SolutionsModN(M,b,Size(A));
                for solution in sols do
                    Add(tps,
                        TrianglePresentationFromSolution(new_geo,A,
                            new_multi_graph,cover,solution)
                    );
                od;
            od;
        fi;
    od;
    if ForAny(tps,tp -> not IsValidTrianglePresentation(tp)) then
        Print("ERROR.\n");
        return fail;
    else
        return tps;
    fi;
end;
