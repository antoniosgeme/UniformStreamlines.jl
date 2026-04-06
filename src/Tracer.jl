"""
    box_intersection!(out, pin, pout, lower, upper)

Compute the point where the segment `pin → pout` crosses the axis-aligned box
`[lower, upper]` and write the result into `out`.
Assumes `pin` is inside the box and `pout` is outside.
"""
function box_intersection!(out::AbstractVector, pin::AbstractVector, pout::AbstractVector,
                           lower::AbstractVector, upper::AbstractVector)
    t = Inf
    for i in eachindex(pin)
        di = pout[i] - pin[i]
        if di != 0
            bound = di > 0 ? upper[i] : lower[i]
            t = min(t, (bound - pin[i]) / di)
        end
    end
    @. out = pin + t * (pout - pin)
    return out
end


"""
    trace!(verts, pos, x, x0, u, stepsize, maxvert, lower, upper, sgn) -> Int

Integrate a single streamline from `x0` into the pre-allocated buffer `verts`
(D × maxvert) using midpoint (RK2) steps. `pos` and `x` are D-length working
vectors. Returns the number of columns written.
"""
function trace!(verts::Matrix{Float64}, pos::Vector{Float64}, x::Vector{Float64},
                x0::AbstractVector, u::F,
                stepsize::Float64, maxvert::Int,
                lower::AbstractVector, upper::AbstractVector, sgn::Int) where F
    verts[:, 1] .= x0

    for i in 1:maxvert-1
        vi = view(verts, :, i)
        pos .= vi
        v = u(pos)
        any(isnan, v) && return i

        # RK2 midpoint
        @. x = vi + sgn * (stepsize / 2) * v
        if !all(j -> lower[j] <= x[j] <= upper[j], eachindex(x))
            box_intersection!(view(verts, :, i+1), vi, x, lower, upper)
            return i + 1
        end

        vm = u(x)
        any(isnan, vm) && return i

        @. x = vi + sgn * stepsize * vm
        if !all(j -> lower[j] <= x[j] <= upper[j], eachindex(x))
            box_intersection!(view(verts, :, i+1), vi, x, lower, upper)
            return i + 1
        end
        verts[:, i+1] .= x
    end
    return maxvert
end


"""
    stream(Val(D), lower, upper, u; kwargs...) -> Matrix{Float64}

Low-level core: compute evenly-spaced streamlines of the vector field `u` over
the axis-aligned box `[lower, upper]`. Works in any dimension D (D=2 or D=3 are
typical).

`u` must be a function `u(p::AbstractVector) -> AbstractVector` returning the
velocity at position `p`.

# Returns
A `D × N` matrix. Individual streamlines are concatenated column-by-column,
separated by columns of `NaN`.

# Keyword arguments
| Keyword            | Default    | Description                                                                 |
|:-------------------|:-----------|:----------------------------------------------------------------------------|
| `min_density`      | `4`        | Seeding grid density (domain divided into `10 × min_density` cells/axis). Higher → more seed candidates → denser coverage. |
| `max_density`      | `10`       | Collision grid density (domain divided into `10 × max_density` cells/axis). Higher → streamlines may pass closer together. |
| `allow_collisions` | `false`    | If `true`, streamlines pass through each other instead of being truncated.  |
| `seeds`            | `nothing`  | Explicit seed points; overrides density grids.                              |
| `min_length`       | `2`        | Discard streamlines with fewer than this many vertices.                     |
| `stepsize`         | adaptive   | Integration step size. By default set to `min(norm(domain) / (10 × max_density × 10), 0.05)`. |
"""
function stream(::Val{D}, lower::Vector{<:Real}, upper::Vector{<:Real}, u::F;
                min_density = 4,
                max_density = 10,
                allow_collisions::Bool = false,
                seeds::Union{Nothing, Tuple{Vararg{AbstractVector}}} = nothing,
                min_length::Int = 2,
                stepsize::Union{Nothing, Float64} = nothing
                ) where {D, F}

    num    = 10
    nstart = ceil(Int, num * min_density)
    nend   = ceil(Int, num * max_density)
    rng    = upper .- lower
    incstart   = rng ./ nstart
    irangecs   = nstart ./ rng .* (1 - eps())
    irangece   = nend   ./ rng .* (1 - eps())
    stepsize = isnothing(stepsize) ? min(norm(rng) / (nend * 10), 0.05) : stepsize
    maxvert = ceil(Int, 2 * norm(rng) / stepsize)

    # Pre-allocate working buffers (reused across all trace! calls)
    verts_f = Matrix{Float64}(undef, D, maxvert)
    verts_b = Matrix{Float64}(undef, D, maxvert)
    pos     = Vector{Float64}(undef, D)
    x_buf   = Vector{Float64}(undef, D)

    dims_start = ntuple(_ -> nstart, Val(D))
    dims_end   = ntuple(_ -> nend,   Val(D))
    startgrid  = falses(dims_start)
    endgrid    = falses(dims_end)

    seed_list = if isnothing(seeds)
        rc_list = collect(CartesianIndices(startgrid))
        shuffle!(rc_list)
        rc_list
    else
        seeds
    end

    # Trim a half-streamline: return the number of valid columns.
    function trim_count(v::AbstractMatrix, ncols::Int; allow_collisions = false)
        for i in 1:ncols
            any(isnan, view(v, :, i)) && return i - 1
        end

        t = CartesianIndex(ntuple(i -> floor(Int, (v[i, 1] - lower[i]) * irangece[i]) + 1, Val(D)))

        for j in 1:ncols
            s = CartesianIndex(ntuple(i -> clamp(floor(Int, (v[i, j] - lower[i]) * irangecs[i]) + 1, 1, nstart), Val(D)))
            startgrid[s] = true

            e = CartesianIndex(ntuple(i -> floor(Int, (v[i, j] - lower[i]) * irangece[i]) + 1, Val(D)))
            if !checkbounds(Bool, endgrid, e)
                return j
            end
            if !allow_collisions && endgrid[e] && e != t
                return j
            end
            endgrid[e] = true
            t = e
        end
        return ncols
    end

    B   = 100_000
    mat = Matrix{Float64}(undef, D, B)
    x0  = Vector{Float64}(undef, D)
    n   = 0

    for item in seed_list
        if item isa CartesianIndex
            for i in 1:D
                x0[i] = lower[i] + (item.I[i] - 0.5) * incstart[i]
            end
            ci = item
        else
            x0 .= item
            ci = CartesianIndex(ntuple(i -> clamp(floor(Int, (x0[i] - lower[i]) / incstart[i]) + 1, 1, nstart), Val(D)))
        end

        if isnothing(seeds)
            startgrid[ci] && continue
        end
        startgrid[ci] = true

        nf = trace!(verts_f, pos, x_buf, x0, u, stepsize, maxvert, lower, upper,  1)
        nb = trace!(verts_b, pos, x_buf, x0, u, stepsize, maxvert, lower, upper, -1)

        if !allow_collisions
            nb_t = trim_count(verts_b, nb)
            nf_t = trim_count(verts_f, nf)
            ncols = nb_t + max(nf_t - 1, 0)

            if ncols > min_length
                while n + ncols + 1 > size(mat, 2)
                    mat = hcat(mat, Matrix{Float64}(undef, D, B))
                end
                # Copy backward half reversed
                @views for j in 1:nb_t
                    mat[:, n + j] .= verts_b[:, nb_t - j + 1]
                end
                # Copy forward half, skip first (same as seed point)
                @views for j in 2:nf_t
                    mat[:, n + nb_t + j - 1] .= verts_f[:, j]
                end
                n += ncols
                mat[:, n+1] .= NaN
                n += 1
            end
        else
            # allow_collisions: build combined path into mat, then trim
            raw = nb + max(nf - 1, 0)
            while n + raw + 1 > size(mat, 2)
                mat = hcat(mat, Matrix{Float64}(undef, D, B))
            end
            # Write backward reversed
            @views for j in 1:nb
                mat[:, n + j] .= verts_b[:, nb - j + 1]
            end
            # Write forward, skip first
            @views for j in 2:nf
                mat[:, n + nb + j - 1] .= verts_f[:, j]
            end
            ncols = trim_count(view(mat, :, n+1:n+raw), raw; allow_collisions = true)

            if ncols > min_length
                n += ncols
                mat[:, n+1] .= NaN
                n += 1
            end
        end
    end
    return mat[:, 1:n]
end
