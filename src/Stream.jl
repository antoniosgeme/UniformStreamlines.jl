
"""
    evenstream(axs::NTuple{D,AbstractVector}, fns::NTuple{D,Function};       kwargs...) -> StreamlineData{D}
    evenstream(axs::NTuple{D,AbstractVector}, fn::Function;                  kwargs...) -> StreamlineData{D}
    evenstream(axs::NTuple{D,AbstractVector}, arrs::NTuple{D,AbstractArray}; kwargs...) -> StreamlineData{D}
    evenstream(xs, ys,     ufn, vfn;          kwargs...) -> StreamlineData{2}
    evenstream(xs, ys,     fn;                kwargs...) -> StreamlineData{2}
    evenstream(xs, ys,     U, V;              kwargs...) -> StreamlineData{2}
    evenstream(xs, ys, zs, ufn, vfn, wfn;    kwargs...) -> StreamlineData{3}
    evenstream(xs, ys, zs, fn;               kwargs...) -> StreamlineData{3}
    evenstream(xs, ys, zs, U, V, W;          kwargs...) -> StreamlineData{3}

Compute evenly-spaced streamlines via the Jobard–Lefer algorithm.

**Input styles:**
- *Separate component functions* — one function per velocity component, each called as
  `ufn(x, y)` / `ufn(x, y, z)`. Zero allocations; the fastest option when the
  components are independent expressions.
- *Single vector-valued function* — one function `fn(x)` that returns all D velocity
  components together. Convenient when the components share computation (e.g. a struct
  carrying parameters). `x` is passed as an `SVector{D,Float64}` so the function
  receives a concrete, stack-allocated input.

  **Return-type performance note:** `evenstream` converts the return value of `fn` to an
  `SVector` internally, but it cannot eliminate an allocation that happens *inside* `fn`.
  Choose the return type carefully:

  | Return type | Example | Allocations per step |
  |:------------|:--------|:---------------------|
  | `Tuple` | `x -> (x[2], -x[1])` | **0** — recommended |
  | `SVector` | `x -> SA[x[2], -x[1]]` | **0** — recommended |
  | `Vector` | `x -> [x[2], -x[1]]` | 1 per call — ~100× slower |

- *Arrays* — pre-computed grids, linearly interpolated at integration points.
  Convention: `U[i, j, ...]` is the first velocity component at the point
  `(axs[1][i], axs[2][j], ...)`, created from the comprehension
  `U = [u(x, y) for x in xs, y in ys]`.

The 2-D/3-D flat forms are convenience wrappers around the N-D tuple form.

# Keyword arguments

| Keyword            | Default    | Description                                                                 |
|:-------------------|:-----------|:----------------------------------------------------------------------------|
| `min_density`      | `4`        | Seeding grid density (domain divided into `10 × min_density` cells/axis). Higher → more seed candidates → denser coverage. |
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

# Separate component functions (flat 2-D form) — zero allocations
result = evenstream(xs, ys, (x,y) -> -y, (x,y) -> x)

# Single vector-valued function — return Tuple or SVector for best performance
result = evenstream(xs, ys, x -> (x[2], -x[1]))         # Tuple  — zero allocations ✓
result = evenstream(xs, ys, x -> SA[x[2], -x[1]])        # SVector — zero allocations ✓
result = evenstream(xs, ys, x -> [x[2], -x[1]])          # Vector — works, but allocates ✗

# Useful when the field carries shared state:
struct MyField; params end
(F::MyField)(x) = (x[2] * F.params, -x[1])              # returns Tuple — fast
result = evenstream(xs, ys, MyField(1.5))

# From pre-computed grids (flat 2-D form)
U = [-y for x in xs, y in ys];  V = [x for x in xs, y in ys]
result = evenstream(xs, ys, U, V)

# N-D tuple form (any dimension)
zs = ts = LinRange(-2, 2, 200)
result = evenstream((xs, ys, zs, ts), ((x,y,z,t) -> -y, (x,y,z,t) -> x, (x,y,z,t) -> z, (x,y,z,t) -> t))

# N-D single function
result = evenstream((xs, ys), x -> (x[2], -x[1]))

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

function evenstream(axs::NTuple{D,AbstractVector}, fn::Function; kwargs...) where D
    lower = Float64[minimum(ax) for ax in axs]
    upper = Float64[maximum(ax) for ax in axs]
    sample = fn(SVector{D,Float64}(ntuple(i -> lower[i], Val(D))))
    length(sample) == D || throw(ArgumentError(
        "fn must return a $D-element vector; got length $(length(sample))"))
    field = p -> (r = fn(SVector{D,Float64}(ntuple(j -> p[j], Val(D)))); SVector{D,Float64}(ntuple(i -> @inbounds(r[i]), Val(D))))
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

function evenstream(xs::AbstractVector, ys::AbstractVector,
                    fn::Function; kwargs...)
    lower = [Float64(minimum(xs)), Float64(minimum(ys))]
    upper = [Float64(maximum(xs)), Float64(maximum(ys))]
    sample = fn(SVector{2,Float64}(lower[1], lower[2]))
    length(sample) == 2 || throw(ArgumentError(
        "fn must return a 2-element vector; got length $(length(sample))"))
    field = p -> (r = fn(SVector{2,Float64}(p[1], p[2])); SVector{2,Float64}(r[1], r[2]))
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

function evenstream(xs::AbstractVector, ys::AbstractVector, zs::AbstractVector,
                    fn::Function; kwargs...)
    lower = Float64[minimum(xs), minimum(ys), minimum(zs)]
    upper = Float64[maximum(xs), maximum(ys), maximum(zs)]
    sample = fn(SVector{3,Float64}(lower[1], lower[2], lower[3]))
    length(sample) == 3 || throw(ArgumentError(
        "fn must return a 3-element vector; got length $(length(sample))"))
    field = p -> (r = fn(SVector{3,Float64}(p[1], p[2], p[3])); SVector{3,Float64}(r[1], r[2], r[3]))
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
    colorize(data::StreamlineData, f) -> Vector

Compute a per-point value along the streamlines for use as a `color` argument
to `streamlines` / `streamlines!`.

The return type matches the return type of `f`:
- `f` returns a `Real`     → `Vector{Float64}` (use with `colormap`)
- `f` returns a `Colorant` → `Vector{<:Colorant}` (used directly as color, no colormap needed)

# The color function `f`

`f` has the signature

    f(pos::AbstractVector, vel::AbstractVector) -> Union{Real, Colorant}

where `pos` is the point position and `vel = field(pos)` is the velocity vector,
both of length D.

**Built-in symbols** (shortcuts for scalar coloring):

| Symbol           | Equivalent function           |
|:-----------------|:------------------------------|
| `:norm`/`:speed` | `(p, v) -> norm(v)`           |
| `:vx` / `:u`     | `(p, v) -> v[1]`              |
| `:vy` / `:v`     | `(p, v) -> v[2]`              |
| `:vz` / `:w`     | `(p, v) -> v[3]`  (3-D only)  |
| `:x`             | `(p, v) -> p[1]`              |
| `:y`             | `(p, v) -> p[2]`              |
| `:z`             | `(p, v) -> p[3]`  (3-D only)  |

# Examples
```julia
# Scalar coloring — pair with colormap
c = colorize(str, :norm)
c = colorize(str, (p, v) -> p[1]^2 + p[2]^2)   # distance² from origin
c = colorize(str, (p, v) -> v[1] / norm(v))     # cos(angle) with x-axis
streamlines!(ax, str; color=c, colormap=:viridis)

# Direct color — no colormap needed
c = colorize(str, (p, v) -> RGBAf(v[1], v[2], 0, 1))   # color by velocity direction
streamlines!(ax, str; color=c)
```
"""
function colorize(data::StreamlineData, f)
    fn    = resolve_color_fn(f)
    field = data.field
    c = [fn(p,field(p)) for p in eachcol(data.paths)]
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
