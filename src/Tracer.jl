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
    trace(x0, u, stepsize, maxvert, lower, upper, sgn) -> Matrix{Float64}

Integrate a single streamline from `x0` using midpoint (RK2) steps of size
`stepsize`. Integration stops when the path exits `[lower, upper]` or reaches
`maxvert` points. `sgn = +1` traces forward, `sgn = -1` backward.

Returns a `D × K` matrix of positions (K ≤ maxvert).
Dimension D is inferred from `length(x0)`.
"""
function trace(x0::AbstractVector, u::Function,
               stepsize::Float64, maxvert::Int,
               lower::AbstractVector, upper::AbstractVector, sgn::Int)
    D = length(x0)
    verts = Matrix{Float64}(undef, D, maxvert)
    verts[:, 1] .= x0

    v = u(x0)
    x = similar(x0, Float64)

    for i in 1:maxvert-1
        v = u(verts[:, i])
        any(isnan, v) && return verts[:, 1:i]

        # RK2 midpoint
        x .= verts[:, i] .+ sgn .* (stepsize / 2) .* v
        if !all(lower .<= x .<= upper)
            verts[:, i+1] .= box_intersection(verts[:, i], x, lower, upper)
            return verts[:, 1:i+1]
        end

        vm = u(x)
        any(isnan, vm) && return verts[:, 1:i]

        x .= verts[:, i] .+ sgn .* stepsize .* vm
        if !all(lower .<= x .<= upper)
            verts[:, i+1] .= box_intersection(verts[:, i], x, lower, upper)
            return verts[:, 1:i+1]
        end
        verts[:, i+1] .= x
    end
    return verts
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
                min_length::Int = 2
                ) where D

    num    = 10
    nstart = ceil(Int, num * min_density)
    nend   = ceil(Int, num * max_density)
    rng    = upper .- lower
    incstart   = rng ./ nstart
    irangecs   = nstart ./ rng .* (1 - eps())
    irangece   = nend   ./ rng .* (1 - eps())
    stepsize = min(rng ./ (nend * 2)..., 0.1) 
    maxvert = 10000

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

    # Trim a half-streamline: cut on NaN, box exit, or density collision.
    function trim(v::AbstractMatrix; allow_collisions = false)
        for i in axes(v, 2)
            any(isnan, view(v, :, i)) && return v[:, 1:i-1]
        end

        t = CartesianIndex(ntuple(i -> floor(Int, (v[i, 1] - lower[i]) * irangece[i]) + 1, Val(D)))

        for j in axes(v, 2)
            s = CartesianIndex(ntuple(i -> clamp(floor(Int, (v[i, j] - lower[i]) * irangecs[i]) + 1, 1, nstart), Val(D)))
            startgrid[s] = true

            e = CartesianIndex(ntuple(i -> floor(Int, (v[i, j] - lower[i]) * irangece[i]) + 1, Val(D)))
            if !checkbounds(Bool, endgrid, e)
                return v[:, 1:j]
            end
            if !allow_collisions && endgrid[e] && e != t
                return v[:, 1:j]
            end
            endgrid[e] = true
            t = e
        end
        return v
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

        vf = trace(x0, u, stepsize, maxvert, lower, upper,  1)
        vb = trace(x0, u, stepsize, maxvert, lower, upper, -1)

        vo = if !allow_collisions
            hcat(trim(vb)[:, end:-1:1], trim(vf)[:, 2:end])
        else
            second = vb[:, end:-1:1]
            first  = vf[:, 2:end]
            trim(hcat(second, first); allow_collisions = true)
        end

        if size(vo, 2) > min_length
            ncols = size(vo, 2)
            while n + ncols + 1 > size(mat, 2)
                mat = hcat(mat, Matrix{Float64}(undef, D, B))
            end
            mat[:, n+1:n+ncols] .= vo
            n += ncols
            mat[:, n+1] .= NaN
            n += 1
        end
    end
    return mat[:, 1:n]
end
