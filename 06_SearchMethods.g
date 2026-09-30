
#https://math.stackexchange.com/questions/4381402/how-do-i-force-gap-give-me-a-true-random-number
Reset(GlobalMersenneTwister,CurrentDateTimeString());;

GreedyTriangleFamilyMultiGraph := function(multigraph)
    local triangles, active_edges, active_triangles, family,
          changed, e, possible, t, tri, counts, mincount;
    triangles := TrianglesMultiGraph(multigraph);
    active_edges := [1..Length(multigraph[2])];
    active_triangles := [1..Length(triangles)];
    family := [];
    while active_edges <> [] do
        repeat
            changed := false;
            for e in ShallowCopy(active_edges) do
                if e in active_edges then
                    possible := Filtered(
                        active_triangles,
                        t -> e in triangles[t]
                    );
                    if Length(possible) = 0 then
                        RemoveSet(active_edges,e);
                        changed := true;
                    elif Length(possible) = 1 then
                        t := possible[1];
                        tri := triangles[t];
                        Add(family,tri);
                        active_edges := Difference(active_edges,tri);
                        active_triangles := Filtered(
                            active_triangles,
                            t -> Intersection(triangles[t],tri) = []
                        );
                        changed := true;
                    fi;
                fi;
            od;
        until not changed;
        if active_edges <> [] then
            counts := List(
                active_edges,
                e -> Length(
                    Filtered(active_triangles,t -> e in triangles[t])
                )
            );
            mincount := Minimum(counts);
            e := active_edges[Position(counts,mincount)];
            possible := Filtered(
                active_triangles,
                t -> e in triangles[t]
            );
            t := possible[1];
            tri := triangles[t];
            Add(family,tri);
            active_edges := Difference(active_edges,tri);
            active_triangles := Filtered(
                active_triangles,
                t -> Intersection(triangles[t],tri) = []
            );
        fi;
    od;
    return family;
end;

GreedyScoreMultiGraph := function(multigraph)
    return
        Sum(List(GreedyTriangleFamilyMultiGraph(multigraph),t -> Length(t)));
end;

BestNeighbourMultiGraph := function(multigraph,tau,perms,checked)
    local best, best_score, move, new_tau, score;
    best := fail;
    best_score := -1;
    for move in perms do
        new_tau := tau*move;
        if not new_tau in checked then
            score := GreedyScoreMultiGraph(
                SingleActionMultiGraph(multigraph,new_tau)
            );
            if score > best_score then
                best := new_tau;
                best_score := score;
            fi;
        fi;
    od;
    return [best,best_score];
end;

BestNeighbourMultiGraphParallel := function(multigraph,tau,perms,checked,number_kernels)
    local number_perms, cuts, packages, worker,
          results, best_score, best;
    number_perms := Length(perms);
    cuts := List(
        [0..number_kernels],
        i -> Int(i*number_perms/number_kernels)
    );
    packages := List(
        [1..number_kernels],
        i -> perms{[cuts[i]+1..cuts[i+1]]}
    );
    worker := function(i)
        return BestNeighbourMultiGraph(
            multigraph,
            tau,
            packages[i],
            checked
        );
    end;
    results := ParListByFork(
        [1..number_kernels],
        worker,
        rec(NumberJobs := number_kernels)
    );
    results := Filtered(results,x -> x[1] <> fail);
    if results = [] then
        return [fail,-1];
    fi;
    best_score := Maximum(List(results,x -> x[2]));
    best := Filtered(results,x -> x[2] = best_score);
    return Random(best);
end;

MultiGraphSearcher := function(multi_graph, number_kernels)
    local Sym, perms, perfect_score, limit,
          tau, graph, family, score, checked,
          best_tau, best_score, very_best_score,
          steps_since_best, total_steps, start_time;
    Sym := SymmetricGroup(multi_graph[1]);
    perms := List(
        Combinations([1..multi_graph[1]],2),
        c -> (c[1],c[2])
    );
    perfect_score := Length(multi_graph[2]);
    limit := 500;
    start_time := CurrentDateTimeString();
    repeat
        tau := Random(Sym);
        checked := [];
        graph := SingleActionMultiGraph(multi_graph,tau);
        family := GreedyTriangleFamilyMultiGraph(graph);
        score := Sum(List(family,t -> Length(t)));
        very_best_score := score;
        steps_since_best := 0;
        total_steps := 0;
        Print("New start with score ",score,".\n");
        while score < perfect_score and steps_since_best <= limit do
            AddSet(checked,tau);
            best_tau := BestNeighbourMultiGraphParallel(
                multi_graph,
                tau,
                perms,
                checked,
                number_kernels
            );
            best_score := best_tau[2];
            if best_tau[1] = fail then
                break;
            fi;
            tau := best_tau[1];
            score := best_score;
            total_steps := total_steps+1;
            steps_since_best := steps_since_best+1;
            if score > very_best_score then
                very_best_score := score;
                steps_since_best := 0;
                Print(
                    "New best score ",score,
                    " after ",total_steps," steps.\n"
                );
            fi;
        od;
        if score = perfect_score then
            graph := SingleActionMultiGraph(multi_graph,tau);
            family := GreedyTriangleFamilyMultiGraph(graph);
        else
            Print(
                "Restarting. Best score was ",
                very_best_score,".\n"
            );
        fi;
    until score = perfect_score;
    Print("Start time: ",start_time,".\n");
    Print("End time: ",CurrentDateTimeString(),".\n");
    return [tau,graph,family];
end;
