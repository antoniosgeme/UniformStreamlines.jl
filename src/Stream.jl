
"""
    evenstream(axs::NTuple{D,AbstractVector}, fns::NTuple{D,Function};       kwargs...) -> StreamlineData{D}
    evenstream(axs::NTuple{D,AbstractVector}, arrs::NTuple{D,AbstractArray}; kwargs...) -> StreamlineData{D}
    evenstream(xs, ys,     ufn, vfn;          kwargs...) -> StreamlineData{2}
    evenstream(xs, ys,     U, V;              kwargs...) -> StreamlineData{2}
    evenstream(xs, ys, zs, ufn, vfn, wfn;    kwargs...) -> StreamlineData{3}
    evenstream(xs, ys, zs, U, V, W;          kwargs...) -> StreamlineData{3}

Compute evenly-spaced streamlines via the Jobard–Lefer algorithm.

**Input styles:**
- *Functions* — velocity components are evaluated on demand; no interpolation overhead.
- *Arrays* — pre-computed grids, linearly interpolated at integration points.
  Convention: `U[i, j, ...]` is the first velocity component at the point
  `(axs[1][i], axs[2][j], ...)`, created from the comprehension
  `U = [u(x, y) for x in xs, y in ys]`.

The 2-D/3-D flat forms are convenience wrappers around the N-D tuple form.

# Keyword arguments

| Keyword            | Default    | Description                                                                 |
|:-------------------|:-----------|:----------------------------------------------------------------------------|
| `min_density`      | `3`        | Seeding grid density (domain divided into `10 × min_density` cells/axis). Higher → more seed candidates → denser coverage. |
| `max_density`      | `10`       | Collision grid density (domain divided into `10 × max_density` cells/axis). Higher → streamlines may pass closer together. |
| `seeds`            | `nothing`  | Explicit seed points (a tuple/vector of D-vectors); overrides density grids. |
| `min_length`       | `2`        | Discard streamlines with fewer than this many vertices.                     |
| `allow_collisions` | `false`    | If `true`, streamlines pass through each other instead of being truncated.  |
| `stepsize`         | adaptive   | Integration step size. By default set to `min(norm(domain) / (10 × max_density × 10), 0.05)`. |

# Returns
A [`StreamlineData{D}`](@ref). Pass to [`colorize`](@ref) and [`streamarrows`](@ref).

# Examples
```julia
xs = LinRange(-2, 2, 200);  ys = LinRange(-2, 2, 200)

# From functions (flat 2-D form)
result = evenstream(xs, ys, (x,y) -> -y, (x,y) -> x)

# From pre-computed grids (flat 2-D form)
U = [-y for x in xs, y in ys];  V = [x for x in xs, y in ys]
result = evenstream(xs, ys, U, V)

# N-D tuple form (any dimension)
zs = ts = LinRange(-2, 2, 200)
result = evenstream((xs, ys, zs, ts), ((x,y,z,t) -> -y, (x,y,z,t) -> x, (x,y,z,t) -> z, (x,y,z,t) -> t))

colors = colorize(result, :norm)
arrows = streamarrows(result; every=15)
```
"""
function evenstream(axs::NTuple{D,AbstractVector}, fns::NTuple{D,Function}; kwargs...) where D
    lower = Float64[minimum(ax) for ax in axs]
    upper = Float64[maximum(ax) for ax in axs]
    field = p -> (args = ntuple(j -> p[j], Val(D)); SVector(ntuple(i -> fns[i](args...), Val(D))))
    paths = stream(Val(D), lower, upper, field; kwargs...)
    return StreamlineData{D}(paths, lower, upper, field)
end

function evenstream(axs::NTuple{D,AbstractVector}, arrs::NTuple{D,AbstractArray{<:Real}}; kwargs...) where D
    expected = Tuple(length(ax) for ax in axs)
    for (k, A) in enumerate(arrs)
        size(A) == expected || throw(ArgumentError(
            "Component $k has size $(size(A)); expected $expected"))
    end
    lower = Float64[minimum(ax) for ax in axs]
    upper = Float64[maximum(ax) for ax in axs]
    itps  = map(A -> linear_interp(axs, A), arrs)
    field = p -> SVector(ntuple(i -> itps[i](p), Val(D)))
    paths = stream(Val(D), lower, upper, field; kwargs...)
    return StreamlineData{D}(paths, lower, upper, field)
end


# ──────────────────────────────────────────────────────────────────────────────
# evenstream — 2-D helpers
# ──────────────────────────────────────────────────────────────────────────────

function evenstream(xs::AbstractVector, ys::AbstractVector,
                    ufn::Function, vfn::Function; kwargs...)
    lower = [Float64(minimum(xs)), Float64(minimum(ys))]
    upper = [Float64(maximum(xs)), Float64(maximum(ys))]
    field = p -> (args = ntuple(j -> p[j], Val(2)); SVector(ufn(args...), vfn(args...)))
    paths = stream(Val(2), lower, upper, field; kwargs...)
    return StreamlineData{2}(paths, lower, upper, field)
end

# 2-D, grid data
function evenstream(xs::AbstractVector, ys::AbstractVector,
                    U::AbstractMatrix{<:Real}, V::AbstractMatrix{<:Real}; kwargs...)
    @assert size(U) == (length(xs), length(ys)) "U must be size (length(xs), length(ys)) = ($(length(xs)), $(length(ys))); got $(size(U))"
    @assert size(V) == size(U) "V must be same size as U"
    lower  = Float64[minimum(xs), minimum(ys)]
    upper  = Float64[maximum(xs), maximum(ys)]
    u_itp  = linear_interp((xs, ys), U)
    v_itp  = linear_interp((xs, ys), V)
    field  = p -> SVector(u_itp(p), v_itp(p))
    paths  = stream(Val(2), lower, upper, field; kwargs...)
    return StreamlineData{2}(paths, lower, upper, field)
end


# ──────────────────────────────────────────────────────────────────────────────
# evenstream — 3-D helpers
# ──────────────────────────────────────────────────────────────────────────────

function evenstream(xs::AbstractVector, ys::AbstractVector, zs::AbstractVector,
                    ufn::Function, vfn::Function, wfn::Function; kwargs...)
    lower = Float64[minimum(xs), minimum(ys), minimum(zs)]
    upper = Float64[maximum(xs), maximum(ys), maximum(zs)]
    field = p -> (args = ntuple(j -> p[j], Val(3)); SVector(ufn(args...), vfn(args...), wfn(args...)))
    paths = stream(Val(3), lower, upper, field; kwargs...)
    return StreamlineData{3}(paths, lower, upper, field)
end

# 3-D, grid data
function evenstream(xs::AbstractVector, ys::AbstractVector, zs::AbstractVector,
                    U::AbstractArray{<:Real,3}, V::AbstractArray{<:Real,3},
                    W::AbstractArray{<:Real,3}; kwargs...)
    expected = (length(xs), length(ys), length(zs))
    @assert size(U) == expected "U must be size $(expected); got $(size(U))"
    @assert size(V) == expected "V must be same size as U"
    @assert size(W) == expected "W must be same size as U"
    lower  = Float64[minimum(xs), minimum(ys), minimum(zs)]
    upper  = Float64[maximum(xs), maximum(ys), maximum(zs)]
    u_itp  = linear_interp((xs, ys, zs), U)
    v_itp  = linear_interp((xs, ys, zs), V)
    w_itp  = linear_interp((xs, ys, zs), W)
    field  = p -> [u_itp(p), v_itp(p), w_itp(p)]
    paths  = stream(Val(3), lower, upper, field; kwargs...)
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

| Symbol          | Equivalent function                        |
|:----------------|:-------------------------------------------|
| `:norm`/`:speed`| `(p, v) -> norm(v)`                        |
| `:vx` / `:u`    | `(p, v) -> v[1]`                           |
| `:vy` / `:v`    | `(p, v) -> v[2]`                           |
| `:vz` / `:w`    | `(p, v) -> v[3]`  (3-D only)               |
| `:x`             | `(p, v) -> p[1]`                           |
| `:y`             | `(p, v) -> p[2]`                           |
| `:z`             | `(p, v) -> p[3]`  (3-D only)               |

# Examples
```julia
colors = colorize(result, :norm)
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
    f === :norm  && return (p, v) -> norm(v)
    f === :speed && return (p, v) -> norm(v)
    f === :vx    && return (p, v) -> v[1]
    f === :u     && return (p, v) -> v[1]
    f === :vy    && return (p, v) -> v[2]
    f === :v     && return (p, v) -> v[2]
    f === :vz    && return (p, v) -> v[3]
    f === :w     && return (p, v) -> v[3]
    f === :x     && return (p, v) -> p[1]
    f === :y     && return (p, v) -> p[2]
    f === :z     && return (p, v) -> p[3]
    throw(ArgumentError("Unknown color symbol :$f. Valid options: :norm, :speed, :vx, :u, :vy, :v, :vz, :w, :x, :y, :z"))
end
resolve_color_fn(f::Function) = f


# ──────────────────────────────────────────────────────────────────────────────
# streamarrows
# ──────────────────────────────────────────────────────────────────────────────

"""
    streamarrows(data::StreamlineData{D}; every=10, spacing=nothing, scale=1.0) -> ArrowData{D}

Extract arrow glyphs from a `StreamlineData` result.

Arrows can be placed in two modes:

1. **Vertex mode** (default) — an arrow every `every` path vertices.
   Fast, but non-uniform when vertex spacing varies.
2. **Arc-length mode** — an arrow every `spacing` units of arc length.
   Set `spacing` to a positive number to enable uniform placement.
   When `spacing` is given, `every` is ignored.

# Keyword arguments
- `every    :: Int   = 10`      — (vertex mode) place an arrow every this many path vertices.
- `spacing  :: Real  = nothing` — (arc-length mode) place an arrow every this many units of
  arc length along each streamline. Overrides `every` when set.
- `scale    :: Real  = 1.0`     — length factor applied to each unit-tangent arrow
  vector; pass this directly to your quiver/arrow plotting call.

# Returns
An [`ArrowData{D}`](@ref) with fields:
- `points`  — `D × N` base positions.
- `vectors` — `D × N` arrow vectors (unit tangent × `scale`).
- `speeds`  — `N`-vector of `‖velocity‖` at each arrow, for color-mapping.
- `indices` — column indices into the original `StreamlineData.paths` matrix
  (nearest vertex to each arrow position).

# Examples
```julia
# Vertex-based (legacy)
arrows = streamarrows(result; every=15, scale=0.05)

# Uniform arc-length spacing
arrows = streamarrows(result; spacing=0.3, scale=0.05)
```
"""
function streamarrows(data::StreamlineData{D}; every::Int=10, spacing::Union{Nothing,Real}=nothing, scale::Real=1.0) where D
    if spacing !== nothing
        return _streamarrows_arclength(data, Float64(spacing), Float64(scale))
    else
        return _streamarrows_vertex(data, every, Float64(scale))
    end
end

# ── Vertex-based placement (original algorithm) ──────────────────────────────

function _streamarrows_vertex(data::StreamlineData{D}, every::Int, scale::Float64) where D
    paths  = data.paths
    field  = data.field
    N      = size(paths, 2)

    pts    = Vector{Float64}[]
    vecs   = Vector{Float64}[]
    speeds = Float64[]
    idxs   = Int[]

    seg_idx = 0
    for i in 1:N
        p = paths[:, i]
        if any(isnan, p)
            seg_idx = 0
            continue
        end
        seg_idx += 1

        if seg_idx > 1 && i < N &&
           !any(isnan, paths[:, i-1]) && 
           !any(isnan, paths[:, i+1]) &&
           seg_idx % every == 0

            tangent = paths[:, i+1] .- paths[:, i-1]
            n = norm(tangent)
            push!(pts, p)
            push!(vecs, tangent .* (scale / n))
            push!(speeds, norm(field(p)))
            push!(idxs, i)
        end
    end

    return _pack_arrowdata(Val(D), pts, vecs, speeds, idxs)
end

# ── Arc-length-based placement ───────────────────────────────────────────────

function _streamarrows_arclength(data::StreamlineData{D}, spacing::Float64, scale::Float64) where D
    paths = data.paths
    field = data.field
    N     = size(paths, 2)

    pts    = Vector{Float64}[]
    vecs   = Vector{Float64}[]
    speeds = Float64[]
    idxs   = Int[]

    # Collect segments (runs of non-NaN columns)
    seg_start = 0
    for i in 1:N+1
        is_nan = i > N || any(isnan, @view(paths[:, i]))
        if is_nan
            if seg_start > 0
                _arclength_segment!(pts, vecs, speeds, idxs,
                                    paths, field, seg_start, i - 1, spacing, scale)
                seg_start = 0
            end
        else
            if seg_start == 0
                seg_start = i
            end
        end
    end

    return _pack_arrowdata(Val(D), pts, vecs, speeds, idxs)
end

function _arclength_segment!(pts, vecs, speeds, idxs,
                             paths::Matrix{Float64}, field,
                             i0::Int, i1::Int, spacing::Float64, scale::Float64)
    n_seg = i1 - i0 + 1
    n_seg < 2 && return  # need at least 2 points for a tangent

    # Compute cumulative arc length
    cumlen = Vector{Float64}(undef, n_seg)
    cumlen[1] = 0.0
    for k in 2:n_seg
        cumlen[k] = cumlen[k-1] + norm(@view(paths[:, i0 + k - 1]) .- @view(paths[:, i0 + k - 2]))
    end

    total = cumlen[end]
    total <= 0 && return

    # Center arrows: offset so they're symmetric about the midpoint
    n_arrows = floor(Int, total / spacing)
    n_arrows < 1 && return
    offset = (total - n_arrows * spacing) / 2 + spacing / 2

    for a in 0:n_arrows-1
        s_target = offset + a * spacing

        # Binary search for the segment containing s_target
        lo = searchsortedlast(cumlen, s_target)
        lo = clamp(lo, 1, n_seg - 1)
        hi = lo + 1

        # Interpolation fraction within edge [lo, hi]
        ds = cumlen[hi] - cumlen[lo]
        t = ds > 0 ? (s_target - cumlen[lo]) / ds : 0.0

        # Interpolated position
        p = (1 - t) .* @view(paths[:, i0 + lo - 1]) .+ t .* @view(paths[:, i0 + hi - 1])

        # Tangent from finite difference of neighbors
        # Use lo/hi for tangent (central difference when possible)
        tl = max(lo - 1, 1)
        th = min(hi + 1, n_seg)
        tangent = @view(paths[:, i0 + th - 1]) .- @view(paths[:, i0 + tl - 1])
        tn = norm(tangent)
        tn <= 0 && continue

        push!(pts, collect(p))
        push!(vecs, collect(tangent .* (scale / tn)))
        push!(speeds, norm(field(p)))
        # Nearest vertex index in the original paths matrix
        push!(idxs, t < 0.5 ? i0 + lo - 1 : i0 + hi - 1)
    end
end

# ── Common packing helper ────────────────────────────────────────────────────

function _pack_arrowdata(::Val{D}, pts, vecs, speeds, idxs) where D
    M = length(pts)
    if M == 0
        return ArrowData{D}(Matrix{Float64}(undef, D, 0),
                            Matrix{Float64}(undef, D, 0),
                            Float64[],
                            Int[])
    end
    point_mat = Matrix{Float64}(undef, D, M)
    vec_mat   = Matrix{Float64}(undef, D, M)
    for j in 1:M
        point_mat[:, j] = pts[j]
        vec_mat[:,   j] = vecs[j]
    end
    return ArrowData{D}(point_mat, vec_mat, speeds, idxs)
end
