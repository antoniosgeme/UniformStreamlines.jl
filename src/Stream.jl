
"""
    stream(axs::NTuple{D,AbstractVector}, fns::NTuple{D,Function};       kwargs...) -> StreamlineData{D}
    stream(axs::NTuple{D,AbstractVector}, arrs::NTuple{D,AbstractArray}; kwargs...) -> StreamlineData{D}
    stream(xs, ys,     ufn, vfn;          kwargs...) -> StreamlineData{2}
    stream(xs, ys,     U, V;              kwargs...) -> StreamlineData{2}
    stream(xs, ys, zs, ufn, vfn, wfn;    kwargs...) -> StreamlineData{3}
    stream(xs, ys, zs, U, V, W;          kwargs...) -> StreamlineData{3}

Compute evenly-spaced streamlines via the Jobard–Lefer algorithm.

**Input styles:**
- *Functions* — velocity components are evaluated on demand; no interpolation overhead.
- *Arrays* — pre-computed grids, linearly interpolated at integration points.
  Convention: `U[i, j, ...]` is the first velocity component at the point
  `(axs[1][i], axs[2][j], ...)`, created from the comprehension
  `U = [u(x, y) for x in xs, y in ys]`.

The 2-D/3-D flat forms are convenience wrappers around the N-D tuple form.

# Keyword arguments
- `min_density` — minimum separation between streamlines.
- `max_density` — maximum separation.
- `seeds` — A NTuple of vectors representing seed points.
- `min_length` — discard streamlines shorter than this many points.
- `unbroken` — Do not truncate streamlines when they intersect.

# Returns
A [`StreamlineData{D}`](@ref). Pass to [`colorize`](@ref) and [`streamarrows`](@ref).

# Examples
```julia
xs = LinRange(-2, 2, 200);  ys = LinRange(-2, 2, 200)

# From functions (flat 2-D form)
result = stream(xs, ys, (x,y) -> -y, (x,y) -> x)

# From pre-computed grids (flat 2-D form)
U = [-y for x in xs, y in ys];  V = [x for x in xs, y in ys]
result = stream(xs, ys, U, V)

# N-D tuple form (any dimension)
zs = ts = LinRange(-2, 2, 200)
result = stream((xs, ys, zs, ts), ((x,y,z,t) -> -y, (x,y,z,t) -> x, (x,y,z,t) -> z, (x,y,z,t) -> t))

colors = colorize(result, :speed)
arrows = streamarrows(result; every=15)
```
"""
function stream(axs::NTuple{D,AbstractVector}, fns::NTuple{D,Function}; kwargs...) where D
    lower = [minimum(ax) for ax in axs]
    upper = [maximum(ax) for ax in axs]
    field = p -> [f(p...) for f in fns]
    paths = evenstream(lower, upper, field; kwargs...)
    return StreamlineData{D}(paths, lower, upper, field)
end

function stream(axs::NTuple{D,AbstractVector}, arrs::NTuple{D,AbstractArray{<:Real}}; kwargs...) where D
    expected = Tuple(length(ax) for ax in axs)
    for (k, A) in enumerate(arrs)
        size(A) == expected || throw(ArgumentError(
            "Component $k has size $(size(A)); expected $expected"))
    end
    lower = [minimum(ax) for ax in axs]
    upper = [maximum(ax) for ax in axs]
    itps  = map(A -> linear_interp(axs, A), arrs)
    field = p -> [itp(p) for itp in itps]
    paths = evenstream(lower, upper, field; kwargs...)
    return StreamlineData{D}(paths, lower, upper, field)
end


# ──────────────────────────────────────────────────────────────────────────────
# stream — 2-D helpers
# ──────────────────────────────────────────────────────────────────────────────

function stream(xs::AbstractVector, ys::AbstractVector,
                ufn::Function, vfn::Function; kwargs...)
    lower = [Float64(minimum(xs)), Float64(minimum(ys))]
    upper = [Float64(maximum(xs)), Float64(maximum(ys))]
    field = p -> [ufn(p...), vfn(p...)]
    paths = evenstream(lower, upper, field; kwargs...)
    return StreamlineData{2}(paths, lower, upper, field)
end

# 2-D, grid data
function stream(xs::AbstractVector, ys::AbstractVector,
                U::AbstractMatrix{<:Real}, V::AbstractMatrix{<:Real}; kwargs...)
    @assert size(U) == (length(xs), length(ys)) "U must be size (length(xs), length(ys)) = ($(length(xs)), $(length(ys))); got $(size(U))"
    @assert size(V) == size(U) "V must be same size as U"
    lower  = [minimum(xs), minimum(ys)]
    upper  = [maximum(xs), maximum(ys)]
    u_itp  = linear_interp((xs, ys), U)
    v_itp  = linear_interp((xs, ys), V)
    field  = p -> [u_itp(p), v_itp(p)]
    paths  = evenstream(lower, upper, field; kwargs...)
    return StreamlineData{2}(paths, lower, upper, field)
end


# ──────────────────────────────────────────────────────────────────────────────
# stream — 3-D, helpers
# ──────────────────────────────────────────────────────────────────────────────

function stream(xs::AbstractVector, ys::AbstractVector, zs::AbstractVector,
                    ufn::Function, vfn::Function, wfn::Function; kwargs...)
    lower = [minimum(xs), minimum(ys), minimum(zs)]
    upper = [maximum(xs), maximum(ys), maximum(zs)]
    field = p -> [ufn(p...), vfn(p...), wfn(p...)]
    paths = evenstream(lower, upper, field; kwargs...)
    return StreamlineData{3}(paths, lower, upper, field)
end

# 3-D, grid data
function stream(xs::AbstractVector, ys::AbstractVector, zs::AbstractVector,
                    U::AbstractArray{<:Real,3}, V::AbstractArray{<:Real,3},
                    W::AbstractArray{<:Real,3}; kwargs...)
    expected = (length(xs), length(ys), length(zs))
    @assert size(U) == expected "U must be size $(expected); got $(size(U))"
    @assert size(V) == expected "V must be same size as U"
    @assert size(W) == expected "W must be same size as U"
    lower  = [minimum(xs), minimum(ys), minimum(zs)]
    upper  = [maximum(xs), maximum(ys), maximum(zs)]
    u_itp  = linear_interp((xs, ys, zs), U)
    v_itp  = linear_interp((xs, ys, zs), V)
    w_itp  = linear_interp((xs, ys, zs), W)
    field  = p -> [u_itp(p), v_itp(p), w_itp(p)]
    paths  = evenstream(lower, upper, field; kwargs...)
    return StreamlineData{3}(paths, lower, upper, field)
end


# ──────────────────────────────────────────────────────────────────────────────
# colorize
# ──────────────────────────────────────────────────────────────────────────────

"""
    colorize(data::StreamlineData, f) -> Vector{Float64}

Compute a per-point scalar value along the streamlines for color-mapping.
Returns a `Vector{Float64}` of length `size(data.paths, 2)` with `NaN` at
separator columns.

# The color function `f`

`f` must have the signature

    f(pos::AbstractVector, vel::AbstractVector) -> Real

where `pos` and `vel` are both length-D vectors.

**Built-in symbols** (shortcuts):

| Symbol      | Equivalent function                        |
|:------------|:-------------------------------------------|
| `:speed`    | `(p, v) -> norm(v)`                        |
| `:vx`       | `(p, v) -> v[1]`                           |
| `:vy`       | `(p, v) -> v[2]`                           |
| `:vz`       | `(p, v) -> v[3]`  (3-D only)               |
| `:x`        | `(p, v) -> p[1]`                           |
| `:y`        | `(p, v) -> p[2]`                           |
| `:z`        | `(p, v) -> p[3]`  (3-D only)               |

# Examples
```julia
colors = colorize(result, :speed)
colors = colorize(result, (p, v) -> p[1]^2 + p[2]^2)   # distance² from origin
colors = colorize(result, (p, v) -> v[1] / norm(v))     # cos(angle) with x-axis
```
"""
function colorize(data::StreamlineData, f)
    fn = resolve_color_fn(f)
    N  = size(data.paths, 2)
    c  = Vector{Float64}(undef, N)
    field = data.field
    for i in 1:N
        p = view(data.paths, :, i)
        if any(isnan, p)
            c[i] = NaN
        else
            vel   = field(p)
            c[i]  = Float64(fn(p, vel))
        end
    end
    return c
end

# Resolve symbol shortcuts to (pos, vel) -> Real functions.
function resolve_color_fn(f::Symbol)
    f === :speed && return (p, v) -> norm(v)
    f === :vx    && return (p, v) -> v[1]
    f === :vy    && return (p, v) -> v[2]
    f === :vz    && return (p, v) -> v[3]
    f === :x     && return (p, v) -> p[1]
    f === :y     && return (p, v) -> p[2]
    f === :z     && return (p, v) -> p[3]
    throw(ArgumentError("Unknown color symbol :$f. Valid options: :speed, :vx, :vy, :vz, :x, :y, :z"))
end
resolve_color_fn(f::Function) = f


# ──────────────────────────────────────────────────────────────────────────────
# streamarrows
# ──────────────────────────────────────────────────────────────────────────────

"""
    streamarrows(data::StreamlineData{D}; every=10, scale=1.0) -> ArrowData{D}

Extract arrow glyphs from a `StreamlineData` result, placed every `every`
steps along each streamline segment.

# Keyword arguments
- `every  :: Int   = 10`   — place an arrow every this many path vertices.
- `scale  :: Real  = 1.0`  — length factor applied to each unit-tangent arrow
  vector; pass this directly to your quiver/arrow plotting call.

# Returns
An [`ArrowData{D}`](@ref) with fields:
- `points`  — `D × N` base positions.
- `vectors` — `D × N` arrow vectors (unit tangent × `scale`).
- `speeds`  — `N`-vector of `‖velocity‖` at each arrow, for color-mapping.

# Example
```julia
arrows = streamarrows(result; every=15, scale=0.05)
# In Makie:
arrows!(ax, arrows.points[1,:], arrows.points[2,:],
            arrows.vectors[1,:], arrows.vectors[2,:],
            color = arrows.speeds)
```
"""
function streamarrows(data::StreamlineData{D}; every::Int=10, scale::Real=1.0) where D
    paths  = data.paths
    field  = data.field
    N      = size(paths, 2)

    pts    = Vector{Float64}[]   # will hold D-vectors
    vecs   = Vector{Float64}[]
    speeds = Float64[]

    seg_idx = 0   # position within the current segment
    for i in 1:N
        p = view(paths, :, i)
        if any(isnan, p)
            seg_idx = 0
            continue
        end
        seg_idx += 1

        # Need a left and right neighbour, both non-NaN.
        if seg_idx > 1 && i < N &&
           !any(isnan, view(paths, :, i-1)) &&
           !any(isnan, view(paths, :, i+1)) &&
           seg_idx % every == 0

            tangent = paths[:, i+1] .- paths[:, i-1]
            n       = norm(tangent)
            if n > 0
                push!(pts,    copy(p))
                push!(vecs,   tangent .* (scale / n))
                push!(speeds, norm(field(p)))
            end
        end
    end

    M = length(pts)
    if M == 0
        return ArrowData{D}(Matrix{Float64}(undef, D, 0),
                            Matrix{Float64}(undef, D, 0),
                            Float64[])
    end

    point_mat = Matrix{Float64}(undef, D, M)
    vec_mat   = Matrix{Float64}(undef, D, M)
    for j in 1:M
        point_mat[:, j] = pts[j]
        vec_mat[:,   j] = vecs[j]
    end

    return ArrowData{D}(point_mat, vec_mat, speeds)
end
