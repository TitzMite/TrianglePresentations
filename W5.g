
Read("01_GeometryTriangles.g");
Read("02_TrianglePresentations.g");
Read("03_Lattices.g");
Read("04_MultiGraphs.g");
Read("05_MultiGraphsAndGeos.g");
Read("06_SearchMethods.g");


W5 := [ 156, [ [ 1, 11, 17, 69, 73, 88 ], [ 2, 12, 18, 70, 74, 89 ], [ 3, 13, 19, 71, 75, 90 ], 
      [ 1, 4, 20, 72, 76, 91 ], [ 2, 5, 21, 73, 77, 79 ], [ 3, 6, 22, 74, 78, 80 ], 
      [ 4, 7, 23, 66, 75, 81 ], [ 5, 8, 24, 67, 76, 82 ], [ 6, 9, 25, 68, 77, 83 ], 
      [ 7, 10, 26, 69, 78, 84 ], [ 8, 11, 14, 66, 70, 85 ], [ 9, 12, 15, 67, 71, 86 ], 
      [ 10, 13, 16, 68, 72, 87 ], [ 11, 13, 34, 56, 62, 113 ], [ 1, 12, 35, 57, 63, 114 ], 
      [ 2, 13, 36, 58, 64, 115 ], [ 1, 3, 37, 59, 65, 116 ], [ 2, 4, 38, 53, 60, 117 ], 
      [ 3, 5, 39, 54, 61, 105 ], [ 4, 6, 27, 55, 62, 106 ], [ 5, 7, 28, 56, 63, 107 ], 
      [ 6, 8, 29, 57, 64, 108 ], [ 7, 9, 30, 58, 65, 109 ], [ 8, 10, 31, 53, 59, 110 ], 
      [ 9, 11, 32, 54, 60, 111 ], [ 10, 12, 33, 55, 61, 112 ], [ 14, 39, 58, 98, 127, 137 ], 
      [ 15, 27, 59, 99, 128, 138 ], [ 16, 28, 60, 100, 129, 139 ], [ 17, 29, 61, 101, 130, 140 ], 
      [ 18, 30, 62, 102, 118, 141 ], [ 19, 31, 63, 103, 119, 142 ], [ 20, 32, 64, 104, 120, 143 ], 
      [ 21, 33, 65, 92, 121, 131 ], [ 22, 34, 53, 93, 122, 132 ], [ 23, 35, 54, 94, 123, 133 ], 
      [ 24, 36, 55, 95, 124, 134 ], [ 25, 37, 56, 96, 125, 135 ], [ 26, 38, 57, 97, 126, 136 ], 
      [ 17, 27, 58, 103, 125, 139 ], [ 18, 28, 59, 104, 126, 140 ], [ 19, 29, 60, 92, 127, 141 ], 
      [ 20, 30, 61, 93, 128, 142 ], [ 21, 31, 62, 94, 129, 143 ], [ 22, 32, 63, 95, 130, 131 ], 
      [ 23, 33, 64, 96, 118, 132 ], [ 24, 34, 65, 97, 119, 133 ], [ 25, 35, 53, 98, 120, 134 ], 
      [ 26, 36, 54, 99, 121, 135 ], [ 14, 37, 55, 100, 122, 136 ], [ 15, 38, 56, 101, 123, 137 ], 
      [ 16, 39, 57, 102, 124, 138 ], [ 8, 123, 131, 139, 144, 145 ], [ 9, 124, 132, 140, 145, 146 ], 
      [ 10, 125, 133, 141, 146, 147 ], [ 11, 126, 134, 142, 147, 148 ], 
      [ 12, 127, 135, 143, 148, 149 ], [ 13, 128, 131, 136, 149, 150 ], 
      [ 1, 129, 132, 137, 150, 151 ], [ 2, 130, 133, 138, 151, 152 ], [ 3, 118, 134, 139, 152, 153 ], 
      [ 4, 119, 135, 140, 153, 154 ], [ 5, 120, 136, 141, 154, 155 ], [ 6, 121, 137, 142, 155, 156 ], 
      [ 7, 122, 138, 143, 144, 156 ], [ 30, 35, 48, 79, 80, 100 ], [ 31, 36, 49, 80, 81, 101 ], 
      [ 32, 37, 50, 81, 82, 102 ], [ 33, 38, 51, 82, 83, 103 ], [ 34, 39, 52, 83, 84, 104 ], 
      [ 27, 35, 40, 84, 85, 92 ], [ 28, 36, 41, 85, 86, 93 ], [ 29, 37, 42, 86, 87, 94 ], 
      [ 30, 38, 43, 87, 88, 95 ], [ 31, 39, 44, 88, 89, 96 ], [ 27, 32, 45, 89, 90, 97 ], 
      [ 28, 33, 46, 90, 91, 98 ], [ 29, 34, 47, 79, 91, 99 ], [ 40, 50, 61, 113, 117, 144 ], 
      [ 41, 51, 62, 105, 114, 145 ], [ 42, 52, 63, 106, 115, 146 ], [ 40, 43, 64, 107, 116, 147 ], 
      [ 41, 44, 65, 108, 117, 148 ], [ 42, 45, 53, 105, 109, 149 ], [ 43, 46, 54, 106, 110, 150 ], 
      [ 44, 47, 55, 107, 111, 151 ], [ 45, 48, 56, 108, 112, 152 ], [ 46, 49, 57, 109, 113, 153 ], 
      [ 47, 50, 58, 110, 114, 154 ], [ 48, 51, 59, 111, 115, 155 ], [ 49, 52, 60, 112, 116, 156 ], 
      [ 25, 61, 66, 97, 115, 129 ], [ 26, 62, 67, 98, 116, 130 ], [ 14, 63, 68, 99, 117, 118 ], 
      [ 15, 64, 69, 100, 105, 119 ], [ 16, 65, 70, 101, 106, 120 ], [ 17, 53, 71, 102, 107, 121 ], 
      [ 18, 54, 72, 103, 108, 122 ], [ 19, 55, 73, 104, 109, 123 ], [ 20, 56, 74, 92, 110, 124 ], 
      [ 21, 57, 75, 93, 111, 125 ], [ 22, 58, 76, 94, 112, 126 ], [ 23, 59, 77, 95, 113, 127 ], 
      [ 24, 60, 78, 96, 114, 128 ], [ 1, 36, 45, 83, 141, 144 ], [ 2, 37, 46, 84, 142, 145 ], 
      [ 3, 38, 47, 85, 143, 146 ], [ 4, 39, 48, 86, 131, 147 ], [ 5, 27, 49, 87, 132, 148 ], 
      [ 6, 28, 50, 88, 133, 149 ], [ 7, 29, 51, 89, 134, 150 ], [ 8, 30, 52, 90, 135, 151 ], 
      [ 9, 31, 40, 91, 136, 152 ], [ 10, 32, 41, 79, 137, 153 ], [ 11, 33, 42, 80, 138, 154 ], 
      [ 12, 34, 43, 81, 139, 155 ], [ 13, 35, 44, 82, 140, 156 ], [ 67, 88, 92, 115, 122, 153 ], 
      [ 68, 89, 93, 116, 123, 154 ], [ 69, 90, 94, 117, 124, 155 ], [ 70, 91, 95, 105, 125, 156 ], 
      [ 71, 79, 96, 106, 126, 144 ], [ 72, 80, 97, 107, 127, 145 ], [ 73, 81, 98, 108, 128, 146 ], 
      [ 74, 82, 99, 109, 129, 147 ], [ 75, 83, 100, 110, 130, 148 ], [ 76, 84, 101, 111, 118, 149 ], 
      [ 77, 85, 102, 112, 119, 150 ], [ 78, 86, 103, 113, 120, 151 ], [ 66, 87, 104, 114, 121, 152 ], 
      [ 70, 83, 94, 107, 128, 153 ], [ 71, 84, 95, 108, 129, 154 ], [ 72, 85, 96, 109, 130, 155 ], 
      [ 73, 86, 97, 110, 118, 156 ], [ 74, 87, 98, 111, 119, 144 ], [ 75, 88, 99, 112, 120, 145 ], 
      [ 76, 89, 100, 113, 121, 146 ], [ 77, 90, 101, 114, 122, 147 ], [ 78, 91, 102, 115, 123, 148 ], 
      [ 66, 79, 103, 116, 124, 149 ], [ 67, 80, 104, 117, 125, 150 ], [ 68, 81, 92, 105, 126, 151 ], 
      [ 69, 82, 93, 106, 127, 152 ], [ 14, 20, 49, 51, 71, 133 ], [ 15, 21, 50, 52, 72, 134 ], 
      [ 16, 22, 40, 51, 73, 135 ], [ 17, 23, 41, 52, 74, 136 ], [ 18, 24, 40, 42, 75, 137 ], 
      [ 19, 25, 41, 43, 76, 138 ], [ 20, 26, 42, 44, 77, 139 ], [ 14, 21, 43, 45, 78, 140 ], 
      [ 15, 22, 44, 46, 66, 141 ], [ 16, 23, 45, 47, 67, 142 ], [ 17, 24, 46, 48, 68, 143 ], 
      [ 18, 25, 47, 49, 69, 131 ], [ 19, 26, 48, 50, 70, 132 ] ] ];

# C13 action on W4
s :=  (  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13)( 14, 15, 16, 17, 18, 19, 20, 21, 22,
      23, 24, 25, 26)( 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39)( 40, 41, 42, 43, 44, 45,
      46, 47, 48, 49, 50, 51, 52)( 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65)
    ( 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78)( 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89,
     90, 91)( 92, 93, 94, 95, 96, 97, 98, 99,100,101,102,103,104)(105,106,107,108,109,110,111,112,113,
     114,115,116,117)(118,119,120,121,122,123,124,125,126,127,128,129,130)(131,132,133,134,135,136,
     137,138,139,140,141,142,143)(144,145,146,147,148,149,150,151,152,153,154,155,156)
    (157,158,159,160,161,162,163,164,165,166,167,168,169)(170,171,172,173,174,175,176,177,178,179,180,
     181,182)(183,184,185,186,187,188,189,190,191,192,193,194,195)(196,197,198,199,200,201,202,203,
     204,205,206,207,208)(209,210,211,212,213,214,215,216,217,218,219,220,221)(222,223,224,225,226,
     227,228,229,230,231,232,233,234)(235,236,237,238,239,240,241,242,243,244,245,246,247)
    (248,249,250,251,252,253,254,255,256,257,258,259,260)(261,262,263,264,265,266,267,268,269,270,271,
     272,273)(274,275,276,277,278,279,280,281,282,283,284,285,286)(287,288,289,290,291,292,293,294,
     295,296,297,298,299)(300,301,302,303,304,305,306,307,308,309,310,311,312) ;

A := Group(s);

# centralizer of <s> restricted to the 156 points:
# C13 wr Sym(12)

porbits := Orbits(A,[1..156]);
c13first_images := [1..156];

for k in [0..12] do
    x := Minimum(porbits[1])^(s^k);
    c13first_images[x] := x^s;
od;

c13first := PermList(c13first_images);


swap12_images := [1..156];

for k in [0..12] do
    x := Minimum(porbits[1])^(s^k);
    y := Minimum(porbits[2])^(s^k);
    swap12_images[x] := y;
    swap12_images[y] := x;
od;

swap12 := PermList(swap12_images);

cycle12_images := [1..156];

for i in [1..12] do
    for k in [0..12] do
        x := Minimum(porbits[i])^(s^k);
        if i < 12 then
            y := Minimum(porbits[i+1])^(s^k);
        else
            y := Minimum(porbits[1])^(s^k);
        fi;
        cycle12_images[x] := y;
    od;
od;

cycle12 := PermList(cycle12_images);

gensC := [c13first,swap12,cycle12];

C := Group(gensC);

Sym12 := SymmetricGroup(12);

phi := GroupHomomorphismByImages(C,Sym12,gensC,
    [(),(1,2),(1,2,3,4,5,6,7,8,9,10,11,12)]);

PreimageInC := function(tau)
    return PreImagesRepresentative(phi,tau);
end;

W5MultiGraph :=
    [ 12,
    [   [1,1], [1,1], [1,2], [1,2], [1,5], [1,9],
        [2,1], [2,3], [2,4], [2,8], [2,12], [2,12],
        [3,2], [3,3], [3,4], [3,6], [3,6], [3,9],
        [4,6], [4,7], [4,7], [4,9], [4,12], [4,12],
        [5,2], [5,2], [5,3], [5,4], [5,7], [5,8],
        [6,1], [6,1], [6,8], [6,10], [6,11], [6,12],
        [7,1], [7,6], [7,6], [7,9], [7,10], [7,11],
        [8,3], [8,4], [8,6], [8,8], [8,10], [8,11],
        [9,2], [9,7], [9,7], [9,8], [9,10], [9,11],
        [10,3], [10,4], [10,5], [10,8], [10,10], [10,11],
        [11,3], [11,4], [11,5], [11,5], [11,9], [11,12],
        [12,5], [12,5], [12,7], [12,9], [12,10], [12,11] ] ];

###################################
#improving computation

#action of normalizer of A on the orbits
Naction := [ [ (), () ], [ (), (3,4)(10,11) ], 
  [ (2,9)(3,7)(5,6)(11,12), (1,2)(3,10,4,11)(7,12) ], 
  [ (2,9)(3,7)(5,6)(11,12), (1,2)(3,11,4,10)(7,12) ], 
  [ (1,4)(2,5)(3,11)(6,9)(7,12)(8,10), (1,7)(2,12)(5,6) ], 
  [ (1,4)(2,5)(3,11)(6,9)(7,12)(8,10), (1,7)(2,12)(3,4)(5,6)(10,11) ], 
  [ (1,4)(2,6)(3,12)(5,9)(7,11)(8,10), (1,12)(2,7)(3,10,4,11)(5,6) ], 
  [ (1,4)(2,6)(3,12)(5,9)(7,11)(8,10), (1,12)(2,7)(3,11,4,10)(5,6) ] ];

TauKey := function(tau)
    return List([1..12],i -> i^tau);
end;

IsNormalizerRepresentativeW5 := function(tau)
    local orbit;
    orbit := List(Naction,p -> p[2]^-1*tau*p[1]);
    return List([1..12],i -> i^tau) = Minimum(List(orbit,TauKey));
end;

compute_W5_tps := function()
    local E, i, tau, new_multi_graph, new_geo,
          covers, cover, eq, M, b, sols,
          solution, tp, tps;
    E := Enumerator(Sym12);
    tps := [];
    for i in [1..Size(Sym12)] do
        tau := E[i];
        if IsNormalizerRepresentativeW5(tau) then
            new_multi_graph :=
                SingleActionMultiGraph(W5MultiGraph,tau);
            if IsPerfectMultiGraphFast(new_multi_graph) then
                new_geo :=
                    SingleAction(W5,PreimageInC(tau));
                covers :=
                    PerfectTriangleCoversMultiGraphFast(new_multi_graph);
                for cover in covers do
                    eq :=
                        LiftEquationSystem(
                            new_geo,A,new_multi_graph,cover);
                    M := List(eq,x -> x[1]);
                    b := List(eq,x -> x[2]);
                    sols :=
                        SolutionsModN(M,b,Size(A));
                    for solution in sols do
                        tp :=
                            TrianglePresentationFromSolution(
                                new_geo,A,new_multi_graph,cover,solution);
                        if not IsValidTrianglePresentation(tp) then
                            Error("Invalid triangle presentation");
                        fi;
                        Add(tps,tp);
                    od;
                od;
            fi;
        fi;
    od;
    return tps;
end;

compute_W5_tps_parallel := function(number_kernels)
    local size,number_packages,cuts,packages,worker,results,tps;
    size := Size(Sym12);
    number_packages := 10*number_kernels;
    cuts := List([0..number_packages],i -> Int(i*size/number_packages));
    packages := List([1..number_packages],i -> [cuts[i]+1,cuts[i+1]]);
    worker := function(interval)
        local E,i,tau,new_multi_graph,new_geo,covers,cover,
              eq,M,b,sols,solution,tp,tps,reps,perfect,
              number_covers,number_solutions;
        E := Enumerator(Sym12);
        tps := [];
        reps := 0;
        perfect := 0;
        number_covers := 0;
        number_solutions := 0;
        for i in [interval[1]..interval[2]] do
            tau := E[i];
            if IsNormalizerRepresentativeW5(tau) then
                reps := reps+1;
                new_multi_graph :=
                    SingleActionMultiGraph(W5MultiGraph,tau);
                if IsPerfectMultiGraphFast(new_multi_graph) then
                    perfect := perfect+1;
                    new_geo :=
                        SingleAction(W5,PreimageInC(tau));
                    covers :=
                        PerfectTriangleCoversMultiGraphFast(new_multi_graph);
                    number_covers :=
                        number_covers+Length(covers);
                    for cover in covers do
                        eq :=
                            LiftEquationSystem(new_geo,A,new_multi_graph,cover);
                        M := List(eq,x -> x[1]);
                        b := List(eq,x -> x[2]);
                        sols :=
                            SolutionsModN(M,b,Size(A));
                        number_solutions :=
                            number_solutions+Length(sols);
                        for solution in sols do
                            tp :=
                                TrianglePresentationFromSolution(
                                    new_geo,A,new_multi_graph,
                                    cover,solution
                                );
                            if not IsValidTrianglePresentation(tp) then
                                Error("Invalid triangle presentation");
                            fi;
                            Add(tps,tp);
                        od;
                    od;
                fi;
            fi;
        od;
        return [tps,reps,perfect,number_covers,number_solutions];
    end;
    results := ParListByFork(
        packages,
        worker,
        rec(NumberJobs := number_kernels)
    );
    tps := Concatenation(List(results,x -> x[1]));
    Print(
        "Normalizer representatives: ",
        Sum(List(results,x -> x[2])),"\n"
    );
    Print(
        "Perfect multigraphs: ",
        Sum(List(results,x -> x[3])),"\n"
    );
    Print(
        "Perfect covers: ",
        Sum(List(results,x -> x[4])),"\n"
    );
    Print(
        "Lift solutions: ",
        Sum(List(results,x -> x[5])),"\n"
    );
    Print(
        "Triangle presentations: ",
        Length(tps),"\n"
    );
    return tps;
end;

find_W5tp := function(number_kernels)
    local x,tau,multi_graph,geo,covers,ncovers,block_size,print_interval,next_print,start,packages,j,lo,hi,worker,results,tp,checked;
    block_size := 5000;
    print_interval := 1000000;
    repeat
        #more kernels do not seem to help
        x := MultiGraphSearcher(W5MultiGraph,4);
        tau := x[1];
        multi_graph := x[2];
        Print("Perfect multigraph found at ",CurrentDateTimeString(),".\n");
        geo := SingleAction(W5,PreimageInC(tau));
        covers := PerfectTriangleCoversMultiGraphFast(multi_graph);
        ncovers := Length(covers);
        Print(ncovers," perfect covers found.\n");
        worker := function(interval)
            local i,cover,eq,M,b,sols;
            for i in [interval[1]..interval[2]] do
                cover := covers[i];
                eq := LiftEquationSystem(geo,A,multi_graph,cover);
                M := List(eq,x -> x[1]);
                b := List(eq,x -> x[2]);
                sols := SolutionsModN(M,b,Size(A));
                if sols <> [] then
                    return TrianglePresentationFromSolution(geo,A,multi_graph,cover,sols[1]);
                fi;
            od;
            return fail;
        end;
        tp := fail;
        start := 1;
        next_print := print_interval;
        while start <= ncovers and tp = fail do
            packages := [];
            for j in [1..number_kernels] do
                lo := start+(j-1)*block_size;
                if lo <= ncovers then
                    hi := Minimum(lo+block_size-1,ncovers);
                    Add(packages,[lo,hi]);
                fi;
            od;
            results := ParListByFork(packages,worker,rec(NumberJobs := number_kernels));
            tp := First(results,y -> y <> fail);
            if tp = fail then
                start := start+number_kernels*block_size;
                checked := Minimum(start-1,ncovers);
                if checked >= next_print then
                    Print("Checked ",checked," / ",ncovers," covers at ",CurrentDateTimeString(),".\n");
                    while next_print <= checked do
                        next_print := next_print+print_interval;
                    od;
                fi;
            fi;
        od;
    until tp <> fail;
    return tp;
end;
