Reset(GlobalMersenneTwister, CurrentDateTimeString());;

BestNeighbourMultiGraph := function(multigraph, tau, perms, checked)
    local best, best_score, move, new_tau, score;
    best := fail;
    best_score := -1;
    for move in perms do
        new_tau := tau * move;
        if not (new_tau in checked) then
            score := ScoreMultiGraph(
                SingleActionMultiGraph(multigraph, new_tau)
            );
            if score > best_score then
                best := new_tau;
                best_score := score;
            fi;
        fi;
    od;
    return [best, best_score];
end;

BestNeighbourMultiGraphParallel := function(multigraph, tau, perms, checked, number_kernels)
    local jobs, cuts, packages, worker, results, best_score;
    if perms = [] then
        return [fail, -1];
    fi;
    jobs := Minimum(number_kernels, Length(perms));
    if jobs = 1 then
        return BestNeighbourMultiGraph(multigraph, tau, perms, checked);
    fi;
    cuts := List([0..jobs], i -> Int(i * Length(perms) / jobs));
    packages := List([1..jobs], i -> perms{[cuts[i]+1..cuts[i+1]]});
    worker := function(i)
        return BestNeighbourMultiGraph(multigraph, tau, packages[i], checked);
    end;
    results := ParListByFork([1..jobs], worker, rec(NumberJobs := jobs));
    results := Filtered(results, x -> x[1] <> fail);
    if results = [] then
        return [fail, -1];
    fi;
    best_score := Maximum(List(results, x -> x[2]));
    return Random(Filtered(results, x -> x[2] = best_score));
end;

# Searches for a matching admitting a perfect quotient cover.
# Returns [tau, matched_multigraph].
MultiGraphSearcher := function(multigraph, number_kernels)
    local S, perms, perfect_score, limit, tau, graph, score,
          checked, best, very_best_score, steps_since_best,
          total_steps;
    # 
    S := SymmetricGroup(multigraph[1]);
    perms := List(Combinations([1..multigraph[1]], 2),c -> (c[1], c[2]));
    perfect_score := Length(multigraph[2]);
    limit := 500;
    while true do
        tau := Random(S);
        graph := SingleActionMultiGraph(multigraph, tau);
        score := ScoreMultiGraph(graph);
        checked := [];
        very_best_score := score;
        steps_since_best := 0;
        total_steps := 0;
        Print("New start with score ", score, ".\n");
        while true do
            if score = perfect_score then
                return [tau, graph];
            fi;
            if perms = [] then
                return fail;
            fi;
            if steps_since_best >= limit then
                break;
            fi;
            AddSet(checked, tau);
            best := BestNeighbourMultiGraphParallel(multigraph, tau, perms, checked, number_kernels);
            if best[1] = fail then
                break;
            fi;
            tau := best[1];
            score := best[2];
            graph := SingleActionMultiGraph(multigraph, tau);
            total_steps := total_steps + 1;
            steps_since_best := steps_since_best + 1;
            if score > very_best_score then
                very_best_score := score;
                steps_since_best := 0;
                Print("New best score ", score," after ", total_steps, " steps.\n");
            fi;
        od;
        Print("Restarting. Best score was ", very_best_score, ".\n");
    od;
end;

# Enumerates the covers and checks their lifting systems in parallel batches.
# Returns [cover, solution], or fail if no cover admits a solution.
FirstLiftableCoverParallel := function(multigraph, differences, n, number_kernels)
    local covers, ncovers, block_size, print_interval, next_print,
          worker, start, packages, j, lo, hi, results, result, checked;
    # 
    block_size := 5000;
    print_interval := 1000000;
    covers := AllPerfectCovers(multigraph);
    ncovers := Length(covers);
    Print(ncovers, " perfect covers found.\n");
    worker := function(interval)
        local i, cover, system, solutions;
        for i in [interval[1]..interval[2]] do
            cover := covers[i];
            system := LiftEquationSystem(multigraph, cover, differences, n);
            if system[1] = [] then
                solutions := [List([1..multigraph[1]], j -> 0)];
            else
                solutions := SolutionsModN(system[1], system[2], n);
            fi;
            if solutions <> [] then
                return [cover, solutions[1]];
            fi;
        od;
        return fail;
    end;
    start := 1;
    next_print := print_interval;
    while start <= ncovers do
        packages := [];
        for j in [1..number_kernels] do
            lo := start + (j-1)*block_size;
            if lo <= ncovers then
                hi := Minimum(lo + block_size - 1, ncovers);
                Add(packages, [lo, hi]);
            fi;
        od;
        if Length(packages) = 1 then
            results := [worker(packages[1])];
        else
            results := ParListByFork(packages, worker,rec(NumberJobs := Length(packages)));
        fi;
        result := First(results, x -> x <> fail);
        if result <> fail then
            return result;
        fi;
        start := start + number_kernels*block_size;
        checked := Minimum(start-1, ncovers);
        if checked >= next_print then
            Print("Checked ", checked, " / ", ncovers," covers at ", CurrentDateTimeString(), ".\n");
            while next_print <= checked do
                next_print := next_print + print_interval;
            od;
        fi;
    od;
    return fail;
end;

# Searches for and returns one triangle presentation.
TrianglePresentationSearcher := function(geo, A, number_kernels)
    local data, search_kernels, candidate, result, matched, lifted, gen;
    if Size(A) = 1 then
        gen := ();
    else
        gen := MinimalGeneratingSet(A)[1];
    fi;
    data := QuotientData(geo, gen);
    if data[1][1] < 20 then
        search_kernels := 4;
    else
        search_kernels := number_kernels;
    fi;
    while true do
        candidate := MultiGraphSearcher(data[1], search_kernels);
        if candidate = fail then
            return fail;
        fi;
        Print("Perfect multigraph found at ",CurrentDateTimeString(), ".\n");
        result := FirstLiftableCoverParallel(candidate[2], data[5], Order(gen), number_kernels);
        if result <> fail then
            matched := SingleActionMultiGraphWithLifts(data[1], data[3], candidate[1]);
            lifted := LiftCover(matched[1], gen,data[2], matched[2],
                data[4],result[1], result[2]);
            return TrianglePresentationFromCover(lifted[1], lifted[2]);
        fi;
        if data[1][1] <= 1 then
            return fail;
        fi;
        Print("No cover admits a solution. Restarting.\n");
    od;
end;