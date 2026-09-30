
Read("01_GeometryTriangles.g");
Read("02_TrianglePresentations.g");
Read("03_Lattices.g");
Read("04_MultiGraphs.g");
Read("05_MultiGraphsAndGeos.g");

H2 := [ 63, [ [ 31, 49, 62 ], [ 32, 43, 63 ], [ 33, 44, 57 ], [ 34, 45, 58 ], 
      [ 35, 46, 59 ], [ 29, 47, 60 ], [ 30, 48, 61 ], [ 12, 49, 59 ], 
      [ 13, 43, 60 ], [ 14, 44, 61 ], [ 8, 45, 62 ], [ 9, 46, 63 ], 
      [ 10, 47, 57 ], [ 11, 48, 58 ], [ 22, 47, 59 ], [ 23, 48, 60 ], 
      [ 24, 49, 61 ], [ 25, 43, 62 ], [ 26, 44, 63 ], [ 27, 45, 57 ], 
      [ 28, 46, 58 ], [ 31, 51, 53 ], [ 32, 52, 54 ], [ 33, 53, 55 ], 
      [ 34, 54, 56 ], [ 35, 50, 55 ], [ 29, 51, 56 ], [ 30, 50, 52 ], 
      [ 6, 7, 12 ], [ 1, 7, 13 ], [ 1, 2, 14 ], [ 2, 3, 8 ], [ 3, 4, 9 ], 
      [ 4, 5, 10 ], [ 5, 6, 11 ], [ 11, 17, 53 ], [ 12, 18, 54 ], 
      [ 13, 19, 55 ], [ 14, 20, 56 ], [ 8, 21, 50 ], [ 9, 15, 51 ], 
      [ 10, 16, 52 ], [ 2, 17, 22 ], [ 3, 18, 23 ], [ 4, 19, 24 ], 
      [ 5, 20, 25 ], [ 6, 21, 26 ], [ 7, 15, 27 ], [ 1, 16, 28 ], 
      [ 26, 38, 41 ], [ 27, 39, 42 ], [ 28, 36, 40 ], [ 22, 37, 41 ], 
      [ 23, 38, 42 ], [ 24, 36, 39 ], [ 25, 37, 40 ], [ 21, 29, 36 ], 
      [ 15, 30, 37 ], [ 16, 31, 38 ], [ 17, 32, 39 ], [ 18, 33, 40 ], 
      [ 19, 34, 41 ], [ 20, 35, 42 ] ] ];


# C5 action on the quadrangle H2
s := (  1,  2,  3,  4,  5,  6,  7)(  8,  9, 10, 11, 12, 13, 14)( 15, 16, 17, 18,
      19, 20, 21)( 22, 23, 24, 25, 26, 27, 28)( 29, 30, 31, 32, 33, 34, 35)
    ( 36, 37, 38, 39, 40, 41, 42)( 43, 44, 45, 46, 47, 48, 49)
    ( 50, 51, 52, 53, 54, 55, 56)( 57, 58, 59, 60, 61, 62, 63)
    ( 64, 65, 66, 67, 68, 69, 70)( 71, 72, 73, 74, 75, 76, 77)
    ( 78, 79, 80, 81, 82, 83, 84)( 85, 86, 87, 88, 89, 90, 91)
    ( 92, 93, 94, 95, 96, 97, 98)( 99,100,101,102,103,104,105)
    (106,107,108,109,110,111,112)(113,114,115,116,117,118,119)
    (120,121,122,123,124,125,126);

A := Group(s);

# the centralizer of <s> restricted to the 63 points
# it has wreath product structure C7 wr Sym(9).
C := Group(
    # generates C7 on the first orbit
    (1,2,3,4,5,6,7),
    # together these two permutations generate Sym9 on the nine C7-orbits
    (1,8)(2,9)(3,10)(4,11)(5,12)(6,13)(7,14),
    (1,8,15,22,29,36,43,50,57)
    (2,9,16,23,30,37,44,51,58)
    (3,10,17,24,31,38,45,52,59)
    (4,11,18,25,32,39,46,53,60)
    (5,12,19,26,33,40,47,54,61)
    (6,13,20,27,34,41,48,55,62)
    (7,14,21,28,35,42,49,56,63)
);

Sym9 := SymmetricGroup(9);

gensC := GeneratorsOfGroup(C);

phi := GroupHomomorphismByImages(C,Sym9,gensC,
    [(),(1,2),(1,2,3,4,5,6,7,8,9)]);

PreimageInC := function(tau)
    return PreImagesRepresentative(phi,tau);
end;

H2MultiGraph :=
    [ 9,
    [   [1,5], [1,5], [1,7],
        [2,2], [2,5], [2,6],
        [3,6], [3,7], [3,9],
        [4,3], [4,7], [4,8],
        [5,1], [5,4], [5,9],
        [6,8], [6,8], [6,9],
        [7,1], [7,2], [7,3],
        [8,4], [8,4], [8,6],
        [9,1], [9,2], [9,3] ] ];

compute_H2_tps := function()
    local tau, multi_graph_covers, new_multi_graph, tps, new_geo,
    sols, equation_data, cover , M, b, solution;
    tps := [];
    for tau in Sym9 do
        new_multi_graph := SingleActionMultiGraph(H2MultiGraph, tau);
        if IsPerfectMultiGraph(new_multi_graph) then
            new_geo := SingleAction(H2, PreimageInC(tau));
            multi_graph_covers := PerfectTriangleCoversMultiGraph(new_multi_graph);
            for cover in multi_graph_covers do
                equation_data := LiftEquationSystem(new_geo,A,new_multi_graph,cover);
                M := List(equation_data, p->p[1]);
                b := List(equation_data, p->p[2]);
                sols := SolutionsModN(M,b, Size(A));
                for solution in sols do
                    Add(tps, TrianglePresentationFromSolution(new_geo,A,new_multi_graph,cover,solution));
                od;
            od;
        fi;
    od;
    if ForAny(tps, tp-> not IsValidTrianglePresentation(tp)) then
        Print("ERROR.\n");
        return fail;
    else
        return tps;
    fi;
end;

Read("H2_reps.g");

LatticesH2C7 := List(H2_reps, tp->TypeRotatingLattice(tp));

RewritingSystemsH2C7 := List(LatticesH2C7, pi -> RewritingSystemTypeRotatingLattice(pi));
