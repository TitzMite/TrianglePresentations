
TrianglePresentationFromCover := function(geo, cover)
    local tp, t, triple;
    tp := [];
    for t in cover do
        if Length(t) = 1 then
            triple := [t[1][1], t[1][1], t[1][1]];
        else
            triple := List(t, e -> e[1]);
        fi;
        Append(tp, [
            triple,
            [triple[2], triple[3], triple[1]],
            [triple[3], triple[1], triple[2]]
        ]);
    od;
    return Set(tp);
end;

SingleActionMultiGraph := function(multigraph, tau)
    return [
        multigraph[1],
        List(multigraph[2], e -> [e[1], e[2]^tau])
    ];
end;

SingleActionMultiGraphWithLifts := function(multigraph, line_lifts, tau)
    local new_multigraph, new_line_lifts;
    new_multigraph := SingleActionMultiGraph(multigraph, tau);
    new_line_lifts := List([1..Length(line_lifts)],j -> line_lifts[j^(tau^-1)]);
    return [new_multigraph, new_line_lifts];
end;

SolutionsModN := function(M,b,n)
    local snf, D, P, Q, rank, m, r, c, options, i, d, g, d1, c1, n1, base, ys;
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
        ys := Concatenation(List(ys, y -> List(options[i], a -> Concatenation(y,[a]))));
    od;
    return List(ys,y -> List(Q, row -> ScalarProduct(row,y) mod n));
end;

ComputeAllSymmetricTp := function(geo, A, N, number_kernels)
    local gen, order, data, multigraph, point_lifts, line_lifts,
          edge_lifts, differences, r, npoints, S,
          point_orbit_index, line_orbit_index, j, k, action,
          size, number_packages, cuts, packages, worker, results, tps;
    # 
    if number_kernels < 1 then
        Error("number_kernels must be positive");
    fi;
    order := Size(A);
    if order = 1 then
        gen := ();
    else
        gen := MinimalGeneratingSet(A)[1];
    fi;
    data := QuotientData(geo, gen);
    multigraph := data[1];
    point_lifts := data[2];
    line_lifts := data[3];
    edge_lifts := data[4];
    differences := data[5];
    r := multigraph[1];
    npoints := geo[1];
    S := SymmetricGroup(r);
    # Identify the point and line orbits.
    point_orbit_index := [];
    line_orbit_index := [];
    for j in [1..r] do
        for k in [0..order-1] do
            point_orbit_index[point_lifts[j]^(gen^k)] := j;
            line_orbit_index[
                (npoints + line_lifts[j])^(gen^k) - npoints
            ] := j;
        od;
    od;
    # Action of N on orbit matchings.
    action := function(tau, h)
        local alpha, beta;
        alpha := PermList(List(point_lifts, p -> point_orbit_index[p^h]));
        beta := PermList(List(line_lifts,l -> line_orbit_index[(npoints + l)^h - npoints]));
        return beta^-1 * tau * alpha;
    end;
    # Split the enumerator into nonempty intervals.
    size := Size(S);
    number_packages := Minimum(size, 10*number_kernels);
    cuts := List([0..number_packages],i -> QuoInt(i*size, number_packages));
    packages := List([1..number_packages],i -> [cuts[i]+1, cuts[i+1]]);
    worker := function(interval)
        local E, i, tau, matched, covers, cover, system, solutions,
              solution, lifted, tps, reps, perfect,
              number_covers, number_solutions;
        E := Enumerator(S);
        tps := [];
        reps := 0;
        perfect := 0;
        number_covers := 0;
        number_solutions := 0;
        for i in [interval[1]..interval[2]] do
            tau := E[i];
            if tau = Minimum(Orbit(N, tau, action)) then
                reps := reps + 1;
                matched := SingleActionMultiGraphWithLifts(multigraph, line_lifts, tau);
                if ScoreMultiGraph(matched[1]) = Length(matched[1][2]) then
                    perfect := perfect + 1;
                    covers := AllPerfectCovers(matched[1]);
                    number_covers := number_covers + Length(covers);
                    for cover in covers do
                        system := LiftEquationSystem(matched[1], cover, differences, order);
                        if system[1] = [] then
                            solutions := Tuples([0..order-1], r);
                        else
                            solutions := SolutionsModN(system[1], system[2], order);
                        fi;
                        number_solutions := number_solutions + Length(solutions);
                        for solution in solutions do
                            lifted := LiftCover(matched[1], gen, point_lifts, matched[2], edge_lifts, cover, solution);
                            Add(tps, TrianglePresentationFromCover(lifted[1], lifted[2]));
                        od;
                    od;
                fi;
            fi;
        od;
        return [tps, reps, perfect, number_covers, number_solutions];
    end;
    if number_kernels = 1 then
        results := List(packages, worker);
    else
        results := ParListByFork(packages, worker,rec(NumberJobs := Minimum(number_kernels, number_packages)));
    fi;
    tps := Concatenation(List(results, x -> x[1]));
    Print("Normalizer representatives: ", Sum(List(results, x -> x[2])), "\n");
    Print("Perfect multigraphs: ", Sum(List(results, x -> x[3])), "\n");
    Print("Perfect covers: ", Sum(List(results, x -> x[4])), "\n");
    Print("Lift solutions: ", Sum(List(results, x -> x[5])), "\n");
    Print("Triangle presentations: ", Length(tps), "\n");
    return tps;
end;