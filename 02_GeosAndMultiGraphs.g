
GeoToMultiGraph := function(geo)
    local edges, i, j;
    edges := [];
    for j in [1..Length(geo[2])] do
        for i in geo[2][j] do
            Add(edges, [i,j]);
        od;
    od;
    Sort(edges);
    return [geo[1], edges];
end;


MultiGraphToGeo := function(multigraph)
    local lines, e;
    lines := List([1..multigraph[1]], i -> []);
    for e in multigraph[2] do
        Add(lines[e[2]], e[1]);
    od;
    for e in lines do
        Sort(e);
    od;
    return [multigraph[1], lines];
end;

# Computes the point and line orbits under gen, ordered by their smallest elements.
# Chooses the smallest element of each orbit as its base lift.
# Constructs one quotient edge per incidence orbit, preserving parallel edges.
# Records each edge’s representative and its cyclic coordinate difference relative to the chosen base lifts.
# Returns [multigraph, point_lifts, line_lifts, edge_lifts, differences]; line lifts use indices in geo[2].
QuotientData := function(geo, gen)
    local n, order, point_orbits, line_orbits, point_lifts, line_lifts,
    line_exponents, edges, edge_lifts, differences, i, j, k, p, l;
    n := geo[1];
    order := Order(gen);
    point_orbits := List(Orbits(Group(gen), [1..n]), Set);
    Sort(point_orbits);
    line_orbits := List(
        Orbits(Group(gen), [n+1..2*n]),
        orbit -> Set(orbit, v -> v-n)
    );
    Sort(line_orbits);
    point_lifts := List(point_orbits, orbit -> orbit[1]);
    line_lifts := List(line_orbits, orbit -> orbit[1]);
    # Exponents relative to the chosen line lifts.
    line_exponents := [];
    for j in [1..Length(line_orbits)] do
        l := n + line_lifts[j];
        for k in [0..order-1] do
            line_exponents[l-n] := k;
            l := l^gen;
        od;
    od;
    edges := [];
    edge_lifts := [];
    differences := [];
    # One representative per incidence orbit.
    for i in [1..Length(point_orbits)] do
        p := point_lifts[i];
        for j in [1..Length(line_orbits)] do
            for l in line_orbits[j] do
                if p in geo[2][l] then
                    Add(edges, [i,j]);
                    Add(edge_lifts, [p,l]);
                    Add(differences, line_exponents[l]);
                fi;
            od;
        od;
    od;
    return [[Length(point_orbits), edges],
        point_lifts,
        line_lifts,
        edge_lifts,
        differences];
end;


# Builds one equation for each triangle in the quotient cover.
# Counts how often each vertex occurs as a target, giving the coefficients of its variable.
# Sets the right-hand side to minus the sum of the triangle's edge differences, modulo n.
# Treats a singleton [e] as (e,e,e), counting its target and difference three times.
# Returns [M,b], representing the system Ms=b mod n.
LiftEquationSystem := function(multigraph, cover, differences, n)
    local M, b, t, triangle, row, e;
    M := [];
    b := [];
    for t in cover do
        if Length(t) = 1 then
            triangle := [t[1], t[1], t[1]];
        else
            triangle := t;
        fi;
        row := List([1..multigraph[1]], j -> 0);
        for e in triangle do
            row[multigraph[2][e][2]] :=
                row[multigraph[2][e][2]] + 1;
        od;
        Add(M, row);
        Add(b, -Sum(List(triangle, e -> differences[e])) mod n);
    od;
    return [M, b];
end;
    
# Constructs the line enumeration determined by the base lifts and the solution.
# Reconstructs new_geo using this enumeration.
# Lifts each quotient triangle using the edge representatives, their differences, and the solution.
# Checks that the three lifted edges satisfy the cyclic compatibility conditions.
# Includes all translates under gen, removes duplicate cyclic rotations, and represents singleton triangles by one arrow.
# Returns [new_geo, lifted_cover], where lifted_cover is a perfect triangle cover of new_geo.
LiftCover := function(multigraph, gen, point_lifts, line_lifts, edge_lifts, cover, solution)
    local n, N, phi2, line_exponents, differences,
          new_lines, lifted_cover, translate,
          j, k, p, l, e, t, ord, h2, h3,
          delta, labels, triangle;
    n := Order(gen);
    N := multigraph[1] * n;
    # phi_{2,c}(gen^k.y_j) = gen^(k+s_j).x_j.
    phi2 := [];
    line_exponents := [];
    for j in [1..multigraph[1]] do
        for k in [0..n-1] do
            l := (N + line_lifts[j])^(gen^k) - N;
            phi2[l] := point_lifts[j]^(gen^(k + solution[j]));
            line_exponents[l] := k;
        od;
    od;
    differences := List(edge_lifts, e -> line_exponents[e[2]]);
    # Translate an original incidence by gen^k.
    translate := function(e, k)
        return [
            e[1]^(gen^k),
            (N + e[2])^(gen^k) - N
        ];
    end;
    # Construct the enumerated geometry F_c.
    new_lines := List([1..N], i -> []);
    for e in edge_lifts do
        for k in [0..n-1] do
            t := translate(e, k);
            Add(new_lines[phi2[t[2]]], t[1]);
        od;
    od;
    new_lines := List(new_lines, Set);
    lifted_cover := [];
    for t in cover do
        if Length(t) = 1 then
            ord := [t[1], t[1], t[1]];
        else
            ord := Minimum([t,[t[2], t[3], t[1]],[t[3], t[1], t[2]]]);
        fi;
        # Exponents of d_1 c_{t(omega_1)}
        # and d_1 c_{t(omega_1)} d_2 c_{t(omega_2)}.
        h2 := differences[ord[1]] + solution[multigraph[2][ord[1]][2]];
        h3 := h2 + differences[ord[2]] + solution[multigraph[2][ord[2]][2]];
        delta := [
            edge_lifts[ord[1]],
            translate(edge_lifts[ord[2]], h2),
            translate(edge_lifts[ord[3]], h3)
        ];
        # Apply (phi_1, phi_{2,c}), with phi_1 the identity.
        labels := List(delta, e -> [e[1], phi2[e[2]]]);
        if not (
            labels[1][2] = labels[2][1] and
            labels[2][2] = labels[3][1] and
            labels[3][2] = labels[1][1]
        ) then
            Error("The supplied solution does not close the lift");
        fi;
        # Include all A-translates.
        for k in [0..n-1] do
            triangle := List(
                labels,
                e -> [e[1]^(gen^k), e[2]^(gen^k)]
            );
            if triangle[1] = triangle[2] and triangle[2] = triangle[3] then
                triangle := [triangle[1]];
            else
                triangle := Minimum([
                    triangle,
                    [triangle[2], triangle[3], triangle[1]],
                    [triangle[3], triangle[1], triangle[2]]
                ]);
            fi;
            AddSet(lifted_cover, triangle);
        od;
    od;
    return [[N, new_lines], lifted_cover];
end;
