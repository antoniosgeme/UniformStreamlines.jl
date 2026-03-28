"""
    box_intersection(pin, pout, lower, upper)

Return the point where the segment `pin → pout` crosses the axis-aligned box
`[lower, upper]`. Assumes `pin` is inside the box and `pout` is outside.
Works for arbitrary dimension D.
"""
function box_intersection(pin::AbstractVector, pout::AbstractVector,
                          lower::AbstractVector, upper::AbstractVector)
    d = pout .- pin
    t = Inf
    for i in eachindex(pin)
        if d[i] != 0
            bound = d[i] > 0 ? upper[i] : lower[i]
            t = min(t, (bound - pin[i]) / d[i])
        end
    end
    return pin .+ t .* d
end


"""
    trace!(verts, pos, x, x0, u, stepsize, maxvert, lower, upper, sgn) -> Int

Integrate a single streamline from `x0` into the pre-allocated buffer `verts`
(D × maxvert) using midpoint (RK2) steps. `pos` and `x` are D-length working
vectors. Returns the number of columns written.
"""
function trace!(verts::Matrix{Float64}, pos::Vector{Float64}, x::Vector{Float64},
                x0::AbstractVector, u::Function,
                stepsize::Float64, maxvert::Int,
                lower::AbstractVector, upper::AbstractVector, sgn::Int)
    verts[:, 1] .= x0

    for i in 1:maxvert-1
        vi = view(verts, :, i)
        pos .= vi
        v = u(pos)
        any(isnan, v) && return i

        # RK2 midpoint
        @. x = vi + sgn * (stepsize / 2) * v
        if !all(j -> lower[j] <= x[j] <= upper[j], eachindex(x))
            verts[:, i+1] .= box_intersection(vi, x, lower, upper)
            return i + 1
        end

        vm = u(x)
        any(isnan, vm) && return i

        @. x = vi + sgn * stepsize * vm
        if !all(j -> lower[j] <= x[j] <= upper[j], eachindex(x))
            verts[:, i+1] .= box_intersection(vi, x, lower, upper)
            return i + 1
        end
        verts[:, i+1] .= x
    end
    return maxvert
end


"""
    evenstream(Val(D), lower, upper, u; kwargs...) -> Matrix{Float64}

Compute evenly-spaced streamlines of the vector field `u` over the axis-aligned
box `[lower, upper]`. Works in any dimension D (D=2 or D=3 are typical).

`u` must be a function `u(p::AbstractVector) -> AbstractVector` returning the
velocity at position `p`.

# Returns
A `D × N` matrix. Individual streamlines are concatenated column-by-column,
separated by columns of `NaN`.

# Keyword arguments
| Keyword        | Default | Description                                      |
|:---------------|:--------|:-------------------------------------------------|
| `min_density`  | `3`     | Coarse grid density for seeding start points     |
| `max_density`  | `10`    | Fine grid density for collision detection        |
| `allow_collisions`     | `false` | If `true`, lines pass through each other         |
| `seeds`        | `nothing` | Explicit seed points; overrides density grids  |
| `min_length`   | `2`     | Discard streamlines with fewer than this many points |
"""
function evenstream(::Val{D}, lower::Vector{<:Real}, upper::Vector{<:Real}, u;
                min_density = 3,
                max_density = 10,
                allow_collisions::Bool = false,
                seeds::Union{Nothing, Tuple{Vararg{AbstractVector}}} = nothing,
                min_length::Int = 2,
                stepsize::Union{Nothing, Float64} = nothing
                ) where D

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
    rc_list    = collect(CartesianIndices(startgrid))
    shuffle!(rc_list)

    start_points = if isnothing(seeds)
        [lower .+ (idx.I .- 0.5) .* incstart for idx in rc_list]
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
    n   = 0

    for x0 in start_points
        ci = CartesianIndex(ntuple(i -> clamp(floor(Int, (x0[i] - lower[i]) / incstart[i]) + 1, 1, nstart), Val(D)))
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
