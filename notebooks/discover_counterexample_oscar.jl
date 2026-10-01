# Discover Plotnick's counterexample using Oscar's matrix group and group algebra.
# Tested with Oscar 1.8.2 and Hecke 0.39.22.
# Run from the repository root with:
# julia --project=notebooks notebooks/discover_counterexample_oscar.jl
# The Hecke subalgebra helper below is an internal API.
using Oscar

function reduced_order_units(O)
    M = maximal_order(O)
    # Keep units factored while computing their maximal-order coordinates.
    U, unitmap = unit_group_fac_elem(M; GRH=false)
    K, order_unitmap = unit_group_fac_elem(O; GRH=false)
    inclusion = hom(K, U, [unitmap \ order_unitmap(g) for g in gens(K)])
    @assert is_injective(inclusion)
    println("Order unit group: ", K)

    # Reduce exponents before evaluating units: unreduced units can be enormous.
    L = hnf(vcat(matrix(inclusion), Hecke.rels(U)))
    L = lll(L[1:ngens(U), 1:ngens(U)])
    return [Hecke.elem_in_algebra(O(evaluate(unitmap(U(L[i:i, :])))))
            for i in 1:nrows(L)]
end

function discover_counterexample_oscar()
    G = SL(2, 5)
    t = G(matrix(base_ring(G), [1 1; 0 1]))
    QG = group_algebra(QQ, G)
    n = order(Int, G)
    @assert n == 120 && order(t) == 5

    # Orbits under conjugation by t and inversion; these operations commute.
    orbits = [[one(G)]]
    seen = Set(first(orbits))
    for g in G
        g in seen && continue
        orbit = unique([t^j * h * t^(-j) for h in (g, inv(g)) for j in 0:4])
        push!(orbits, orbit)
        union!(seen, orbit)
    end
    B = [sum(QG(g) for g in orbit) for orbit in orbits]
    N = sum(B)
    @assert length(orbits) == 17

    # Realize QG/(N) inside QG using the complementary central idempotent p.
    p = one(QG) - inv(QQ(n))*N
    quotient_basis = [b*p for b in B[2:end]]
    A, embed = Hecke._subalgebra(QG, quotient_basis)
    @assert dim(A) == 16 && is_commutative(A)
    # _subalgebra changes basis, so use the original orbit lattice for the order.
    preimages = [has_preimage_with_preimage(embed, b) for b in quotient_basis]
    @assert all(first, preimages)
    O = order(A, last.(preimages); isbasis=true)
    println("Dimension: ", dim(A), "; discriminant: ", discriminant(O))

    units = reduced_order_units(O)
    # Our integral lifts have zero identity coefficient before projection by p.
    # Hence the augmentation modulo n is -n times that projected coefficient.
    augmentation_mod_n(x) = Int(mod(ZZ(-n*embed(x)[one(G)]), n))
    residues = augmentation_mod_n.(units)

    # Multiply (residue, unit) pairs, identifying pairs by residue only.
    representatives = Dict(Hecke.closure(
        collect(zip(residues, units)),
        (a, b) -> (mod(a[1]*b[1], n), a[2]*b[2]),
        (1, one(A)); eq=(a, b) -> a[1] == b[1]))
    println("Augmentation image: ", sort!(collect(keys(representatives))))
    @assert haskey(representatives, 49)

    x = representatives[49]
    u = embed(x) + QQ(49, n)*N
    v = embed(inv(x)) + QQ(49, n)*N
    @assert all(isone(denominator(c)) for w in (u, v) for c in coefficients(w))
    @assert sum(coefficients(u)) == sum(coefficients(v)) == 49
    @assert QG(t)*u == u*QG(t) && QG(t)*v == v*QG(t)
    @assert u*v == v*u == one(QG) + 20*N
    println("Verified integral lifts: augmentation 49, t-invariant, uv = vu = 1 + 20N.")
    println("Orbit coefficients of u: ", [ZZ(u[first(o)]) for o in orbits])
    println("Orbit coefficients of v: ", [ZZ(v[first(o)]) for o in orbits])
    return (; G, t, QG, N, O, orbits, u, v)
end

# Retain the actual group-algebra elements for inspection after include(...).
result = discover_counterexample_oscar()
