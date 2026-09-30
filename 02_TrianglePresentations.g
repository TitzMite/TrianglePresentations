#The triangle families are of the form like the output of the functions in the Score files
IsValidTriangleFamily := function(geo, triangle_family)
    local t, valid, a, arrows;
    ###
    valid := true;
    ###
    #We check if every triangle is actually a triangle in geo
    for t in triangle_family do
        for a in t do
            if not a[1] in geo[2][a[2]] then
                valid := false;
            fi;
        od;
    od;
    ###
    #We check that no arrow is contained in two different triangles
    arrows := [];
    for t in triangle_family do
        for a in t do
            Add(arrows, a);
        od;
    od;
    if not IsDuplicateFree(arrows) then
        valid := false;
    fi;
    ###
    return valid;
end;

PresentationFormsTriangle := function(t)
    local x, y, z;
    if Size(t) = 3 then
        x := t[1][1];
        y := t[2][1];
        z := t[3][1];
        return [ [x,y,z], [y,z,x], [z,x,y] ];
    elif Size(t) = 1 then
        x := t[1][1];
        return [[x,x,x]];
    else
        return fail;
    fi;
end;

TrianglePresentationFamily := function(geo, triangle_family)
    local triangle_presentation, t;
    if not IsPerfectGeo(geo) then
        Print("The geometry is not perfect...\n");
        triangle_presentation := fail;
    elif not IsValidTriangleFamily(geo,triangle_family) then
        Print("The triangle family is not valid...\n");
        triangle_presentation := fail;
    elif not OptimalScore(geo) = ScoreTriangleFamily(triangle_family) then
        Print("The triangle family is not perfect...\n");
        triangle_presentation := fail;
    else
        triangle_presentation := [];
        for t in triangle_family do
            Append(triangle_presentation, PresentationFormsTriangle(t));
        od;
        triangle_presentation := SSortedList(triangle_presentation);
    fi;
    return triangle_presentation;
end;

#returns all optimal triangle presentations for a geo
TrianglePresentationsGeo := function(geo)
    local is_perfect, triangle_families, triangle_presentations;
    if not IsPerfectGeo(geo) then
        triangle_presentations := fail;
    else
        triangle_families := OptimalDisjointTriangleFamilies(geo);
        triangle_presentations := List(triangle_families, f->TrianglePresentationFamily(geo, f));
    fi;
    return triangle_presentations;
end;

###
###
###next functions are to analyze further
###
###

IsValidTrianglePresentation := function(triangle_presentation)
    local closed, valid, gamma, delta;
    valid := true;
    if not IsDuplicateFree(triangle_presentation) then
        Print("The set of triangles is not duplicate-free...\n");
        valid := false;
        #we now check if the triangle_presentation is closed under cyclic permutation
    elif not ForAll(triangle_presentation, t->[t[2],t[3],t[1]] in triangle_presentation) then
        Print("The set of triangles is not closed under cyclic permutation...\n");
        valid := false;
        return valid;
    fi;
    return valid;
end;

#the link
GrapeGraphTrianglePresentation := function(triangle_presentation)
    local edges, points, n, gamma;
    if not IsValidTrianglePresentation(triangle_presentation) then
        Print("The triangle presentation is not valid...\n");
        gamma := fail;
    else
        points := List(triangle_presentation, t->t[1]);
        points := SSortedList(points);
        n := Size(points);
        edges := List(triangle_presentation, t-> [Position(points, t[1]), Position(points, t[2])]);
        edges := List(edges, a->[a[1], a[2]+n]);
        edges := Concatenation(edges, List(edges, a->[a[2],a[1]]));
        gamma := EdgeOrbitsGraph(Group(()), edges, 2*n);
    fi;
    return gamma;
end;

IsTrianglePresentationForGraph := function(tp, gamma)
    return not GraphIsomorphism(GrapeGraphTrianglePresentation(tp), gamma) = fail;
end;

# this function can be used to check for general isomorphisms
# rather than only type-preservings
# we include it here, but have no use for it in our project
RevertedTrianglePresentation := function(tp)
    return SSortedList(List(tp), t-> [t[3],t[2],t[1]]);
end;

TypeAutomorphismGroupTrianglePresentation := function(tp)
    local gamma,aut,aut_tp;
    gamma := GrapeGraphTrianglePresentation(tp);
    aut := AutomorphismGroup(gamma);
    aut_tp := Filtered(aut,a -> OnSetsTuples(tp,a) = tp);
    aut_tp := List(aut_tp,p -> RestrictedPerm(p,[1..Maximum(List(tp,t -> t[1]))]));
    aut_tp := Group(aut_tp);
    aut_tp := Group(SmallGeneratingSet(aut_tp));
    StructureDescription(aut_tp);
    return aut_tp;
end;

############

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