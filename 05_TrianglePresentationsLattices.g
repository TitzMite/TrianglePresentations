CheapNonIsomorphismCheck := function(tp1,tp2)
    local n1,n2,diag1,diag2,double1,double2;
    n1 := Maximum(Flat(tp1));
    n2 := Maximum(Flat(tp2));
    if n1 <> n2 then
        return true;
    fi;
    diag1 := Number(tp1,t -> t[1]=t[2] and t[2]=t[3]);
    diag2 := Number(tp2,t -> t[1]=t[2] and t[2]=t[3]);
    if diag1 <> diag2 then
        return true;
    fi;
    double1 := Number(tp1,t -> t[1]=t[2] and t[2]<>t[3]);
    double2 := Number(tp2,t -> t[1]=t[2] and t[2]<>t[3]);
    if double1 <> double2 then
        return true;
    fi;
    return fail;
end;

TypePreservingIsomorphismCheck := function(tp1,tp2)
    local n,MakeGraph,InducedPerm,gamma1,gamma2,iso,sigma,aut1,enum,a;
    if CheapNonIsomorphismCheck(tp1,tp2)=true then
        return false;
    fi;
    n := Maximum(Flat(tp1));
    MakeGraph := function(tp,n)
        local adj,t,a,b,c,gamma;
        adj := List([1..3*n],i -> List([1..3*n],j -> false));
        for t in tp do
            a := t[1]; b := t[2]; c := t[3];
            adj[a][n+b] := true; adj[n+b][a] := true;
            adj[n+b][2*n+c] := true; adj[2*n+c][n+b] := true;
            adj[2*n+c][a] := true; adj[a][2*n+c] := true;
        od;
        gamma := Graph(Group(()),[1..3*n],OnPoints,function(x,y) return adj[x][y]; end,true);
        return rec(graph:=gamma,colourClasses:=[[1..n],[n+1..2*n],[2*n+1..3*n]]);
    end;
    InducedPerm := function(p,n)
        local images,i,sigma;
        images := List([1..n],i -> i^p);
        if not ForAll(images,x -> x>=1 and x<=n) then
            return fail;
        fi;
        sigma := PermList(images);
        for i in [1..n] do
            if (n+i)^p<>n+i^sigma or (2*n+i)^p<>2*n+i^sigma then
                return fail;
            fi;
        od;
        return sigma;
    end;
    gamma1 := MakeGraph(tp1,n);
    gamma2 := MakeGraph(tp2,n);
    iso := GraphIsomorphism(gamma1,gamma2);
    if iso=fail then
        return false;
    fi;
    sigma := InducedPerm(iso,n);
    if sigma<>fail and OnSetsTuples(tp1,sigma)=tp2 then
        return true;
    fi;
    aut1 := AutomorphismGroup(gamma1);
    enum := Enumerator(aut1);
    for a in enum do
        sigma := InducedPerm(a*iso,n);
        if sigma<>fail and OnSetsTuples(tp1,sigma)=tp2 then
            return true;
        fi;
    od;
    return false;
end;

ClustersTrianglePresentationsTypePreserving := function(triangle_presentations)
    local clusters,i,c;
    clusters := [];
    for i in [1..Length(triangle_presentations)] do
        c := PositionProperty(clusters,cluster ->
            CheapNonIsomorphismCheck(triangle_presentations[i],triangle_presentations[cluster[1]])=fail);
        if c=fail then
            Add(clusters,[i]);
        else
            Add(clusters[c],i);
        fi;
    od;
    return clusters;
end;

RepresentativesTrianglePresentationsTypePreserving := function(triangle_presentations)
    local clusters,representatives,cluster,reps,i,j;
    clusters := ClustersTrianglePresentationsTypePreserving(triangle_presentations);
    representatives := [];
    for cluster in clusters do
        reps := [];
        for i in cluster do
            if ForAll(reps,j -> not TypePreservingIsomorphismCheck(triangle_presentations[i],triangle_presentations[j])) then
                Add(reps,i);
            fi;
        od;
        Append(representatives,reps);
    od;
    return triangle_presentations{representatives};
end;

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
        return fail;
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
    innerColours := Filtered([[1], Difference([1..Length(B1)], [1])],C -> Length(C) > 0);
    A1 := AutGroupGraph(rec(graph := delta,colourClasses := innerColours));
    # Translate the A1-orbits back to gamma's vertex numbering.
    colours := List(Orbits(A1, [1..Length(B1)]),O -> Set(O, i -> B1[i]));
    # All vertices outside B_1 receive one additional colour.
    outer := Difference([1..gamma.order], B1);
    if Length(outer) > 0 then
        Add(colours, outer);
    fi;
    # Keep the known centre stabilizer associated with the graph.
    return AutGroupGraph(rec(graph := NewGroupGraph(H, gamma),colourClasses := colours));
end;