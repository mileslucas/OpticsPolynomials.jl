export Zernike,
       zernike_norm,
       zernike_nm_to_fringe,
       zernike_nm_to_ansi_j,
       zernike_ansi_j_to_nm,
       zernike_noll_to_nm,
       zernike_fringe_to_nm,
       zernike_zero_separation

"""
    Zernike([T, ]n, m, norm=true)

Zernike polynomial of orders `n` and `m`. If `norm` is true the output will be normalized. The output type can be specified with `T`, which will default to `Float64`.

# Examples

```jldoctest
julia> z = Zernike(10, 3)
Zernike{Float64}(n=10, m=3)

julia> z(0.8, 0.2)
-1.6860955965619917

julia> Zernike(BigFloat, 10, 3)(0.8, 0.2)
-1.68609559656199170518675600760616362094879150390625
```
"""
struct Zernike{T<:AbstractFloat,F1,F2}
    n::Int
    m::Int
    basis::F1
    trig_func::F2
    Zernike{T}(n, m, basis, trig_func) where {T} = new{T,typeof(basis),typeof(trig_func)}(n, m, basis, trig_func)
end

function Zernike(T, n, m, norm::Bool=true)
    norm_val = norm ? zernike_norm(n, m) : 1.0
    trig_func = m < 0 ? sin : cos
    basis = x -> norm_val * jacobi((n - m) ÷ 2, 0, abs(m), x)
    Zernike{T}(n, m, basis, trig_func)
end

Zernike(n, m, norm::Bool=true) = Zernike(Float64, n, m, norm)

function Base.show(io::IO, z::Zernike{T}) where T
    print(io, "Zernike{$T}(n=$(z.n), m=$(z.m))")
end

function (z::Zernike{T})(ρ, θ) where T
    x = ρ^2 - 1
    return T(ρ^abs(z.m) * z.basis(x) * z.trig_func(z.m * θ))
end
    


"""
    kronecker(i,j)

1 if i==j, else 0; mathematical kronecker function
"""
kronecker(i, j) = Int(i == j)

"""
    zernike_norm(n, m)

Norm of Zernike polynomial of radial order n, azimuthal order m.

The norm is the average squared distance to zero.  By multiplying a zernike
value by the norm, the term is given unit stdev or RMS.
"""
function zernike_norm(n, m)
    num = √(2 * (n+1)) / (1 + kronecker(m, 0))
end

"""
    zernike_nm_to_fringe(n, m)

Map (n,m) ANSI indices to a single fringe index.
"""
function zernike_nm_to_fringe(n, m)
    term1 = (1 + (n + abs(m))/2)^2
    term2 = 2*abs(m)
    term3 = (1 + sign(m)) / 2
    return int(term2 - term2 - term3) + 1
end

"""
    zernike_nm_to_ansi_j(n, m)

Map (n,m) ANSI indices to a single ANSI j index.

See also:
    - [`zernike_ansi_j_to_nm`](@ref) (reciprocal of this function)
"""
function zernike_nm_to_ansi_j(n, m)
    return int((n * (n + 2) + m) / 2)
end

"""
    zernike_ansi_to_ansi_j(n, m)

Map (n,m) ANSI indices to a single ANSI j index.

See also:
    - [`zernike_nm_to_ansi_j`](@ref) (reciprocal of this function)
"""
function zernike_ansi_j_to_nm(j)
    n = int(ceil((-3 + √(9 + 8j))/2))
    m = 2j - n * (n + 2)
    return n, m
end

"""
    zernike_noll_to_nm(j)

Map j Noll index to ANSI (n,m) indices.
"""
function zernike_noll_to_nm(j)
    n = int(ceil((-1 + √(1 + 8j))/2) - 1)
    if n == 0
        m = 0
    else
        nseries = int((n+1) * (n+2) / 2)
        residual = j - nseries - 1

        if isodd(j)
            sign = -1
        else
            sign = 1
        end

        if isodd(n)
            ms = [1,1]
        else
            ms = [0]
        end

        for i=0:n÷2
            push!(ms, ms[end]+2)
            push!(ms, ms[end])
        end

        m = ms[residual] * sign
    end
    return n, m
end

"""
    zernike_fringe_to_nm(j)

Map j Fringe index to ANSI (n,m) indices.
"""
function zernike_fringe_to_nm(j)
    m_n = 2 * ceil(√j - 1)
    g_s = (m_n / 2)^2 + 1
    n = m_n / 2 + floor((j-g_s)/2)
    m = m_n - n * (1 - mod(j-g_s, 2) * 2)
    return int(n), int(m)
end

"""
    zernike_zero_separation(n)

Minimum zero separation of Zernike polynomial of radial order n.  Useful for
computing sample count requirements.
"""
function zernike_zero_separation(n)
    return 1 / n^2
end

"""
    zernike(n, m, ρ, θ[; norm])

Zernike polynomial of radial order n and azimuthal order m, evaluated at the
point (ρ, θ).  No normalization is required of (ρ, θ), though the polynomials
are orthogonal only over the unit disk.

norm is a boolean flag indicating whether the result should be orthonormalized
(scaled to unit RMS) or not.
"""
function zernike(n, m, ρ, θ; norm::Bool=true)
    x = ρ^2 - 1
    n_j = (n - m) / 2
    am = abs(m)
    # α=0, β=|m|
    # there is a second syntax where you have x reversed, 1 - ρ^2,
    # in which ase you swap α and β.  It makes absolutely no difference
    out = jacobi(n_j, 0, am, x)
    if m != 0
        if sign(m) == -1
            f = sin
        else
            f = cos
        end
        out *= (ρ^am * f(m*θ))
    end
	if norm
		out *= zernike_norm(n,m)
	end
    return out
end
