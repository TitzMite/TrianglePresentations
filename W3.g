
Read("01_GeometryTriangles.g");
Read("02_TrianglePresentations.g");
Read("03_Lattices.g");
Read("04_MultiGraphs.g");
Read("05_MultiGraphsAndGeos.g");

W3 := [ 40, [ [ 1, 5, 17, 27 ], [ 1, 2, 18, 28 ], [ 2, 3, 19, 29 ], [ 3, 4, 20, 30 ], [ 4, 5, 16, 26 ], [ 6, 15, 26, 33 ], 
      [ 7, 11, 27, 34 ], [ 8, 12, 28, 35 ], [ 9, 13, 29, 31 ], [ 10, 14, 30, 32 ], [ 9, 11, 16, 36 ], [ 10, 12, 17, 37 ], 
      [ 6, 13, 18, 38 ], [ 7, 14, 19, 39 ], [ 8, 15, 20, 40 ], [ 1, 21, 30, 31 ], [ 2, 22, 26, 32 ], [ 3, 23, 27, 33 ], 
      [ 4, 24, 28, 34 ], [ 5, 25, 29, 35 ], [ 6, 17, 19, 24 ], [ 7, 18, 20, 25 ], [ 8, 16, 19, 21 ], [ 9, 17, 20, 22 ], 
      [ 10, 16, 18, 23 ], [ 21, 25, 33, 37 ], [ 21, 22, 34, 38 ], [ 22, 23, 35, 39 ], [ 23, 24, 31, 40 ], 
      [ 24, 25, 32, 36 ], [ 1, 15, 36, 39 ], [ 2, 11, 37, 40 ], [ 3, 12, 36, 38 ], [ 4, 13, 37, 39 ], [ 5, 14, 38, 40 ], 
      [ 6, 11, 30, 35 ], [ 7, 12, 26, 31 ], [ 8, 13, 27, 32 ], [ 9, 14, 28, 33 ], [ 10, 15, 29, 34 ] ] ];


# C5 action on the quadrangle W3
s := (1,2,3,4,5)(6,7,8,9,10)(11,12,13,14,15)(16,17,18,19,20)(21,22,23,24,25)(26,27,28,29,30)(31,32,33,34,35)(36,37,38,39,
    40)(41,42,43,44,45)(46,47,48,49,50)(51,52,53,54,55)(56,57,58,59,60)(61,62,63,64,65)(66,67,68,69,70)(71,72,73,74,75)(76,
    77,78,79,80);

A := Group(s);

# the centralizer of <s> (restricted on points) in Sym(40)
# it has wreath product structure C5 wr Sym(8).
C := Group(
    #next perm generates a C5 acting on the first orbit
    (1,2,3,4,5), 
    #next two perms generate a Sym8 action on the eight orbits
    (1,6)(2,7)(3,8)(4,9)(5,10),
    (1,6,11,16,21,26,31,36) (2,7,12,17,22,27,32,37) (3,8,13,18,23,28,33,38) (4,9,14,19,24,29,34,39) (5,10,15,20,25,30,35,40) 
);

Sym8 := SymmetricGroup(8);

gensC := GeneratorsOfGroup(C);

phi := GroupHomomorphismByImages(
        C,Sym8,gensC,[(),(1,2),(1,2,3,4,5,6,7,8)]);

PreimageInC := function(tau)
    return PreImagesRepresentative(phi, tau);
end;

W3MultiGraph :=
    [ 8,
    [   [1,1], [1,1], [1,4], [1,7],
        [2,2], [2,3], [2,5], [2,8],
        [3,2], [3,3], [3,7], [3,8],
        [4,1], [4,3], [4,5], [4,5],
        [5,4], [5,5], [5,6], [5,6],
        [6,1], [6,2], [6,4], [6,8],
        [7,2], [7,4], [7,6], [7,8],
        [8,3], [8,6], [8,7], [8,7] ] ];

compute_W3_tps := function()
    local tau, multi_graph_covers, new_multi_graph, tps, new_geo,
    sols, equation_data, cover , M, b, solution;
    tps := [];
    for tau in Sym8 do
        new_multi_graph := SingleActionMultiGraph(W3MultiGraph, tau);
        if IsPerfectMultiGraph(new_multi_graph) then
            new_geo := SingleAction(W3, PreimageInC(tau));
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

Read("W3_reps.g");

LatticesW3C5 := List(W3_reps, tp->TypeRotatingLattice(tp));

RewritingSystemsW3C5 := List(LatticesW3C5, pi -> RewritingSystemTypeRotatingLattice(pi));

BallsW3C5 := List([1..Size(W3_reps)], i-> CayleyBall(LatticesW3C5[i], RewritingSystemsW3C5[i], 2, A));