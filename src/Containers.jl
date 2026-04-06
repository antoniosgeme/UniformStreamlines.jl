"""
    StreamlineData{D}

Holds the output of `evenstream`. `D` is the spatial dimension.

# Fields
- `paths  :: Matrix{Float64}` — `D × N` matrix. Streamlines are stored
  column-by-column and separated by columns of `NaN`. Pass `paths[1,:]` and
  `paths[2,:]` (and `paths[3,:]` in 3-D) directly to a plotting library.
- `lower  :: Vector{Float64}` — lower corner of the domain box.
- `upper  :: Vector{Float64}` — upper corner of the domain box.
- `field  :: Function`        — the interpolated velocity function
  `field(p::AbstractVector) -> Vector{Float64}`. Useful for querying
  velocities at arbitrary points, and required internally by `colorize` and
  `streamarrows`.
"""
struct StreamlineData{D}
    paths  :: Matrix{Float64}
    lower  :: Vector{Float64}
    upper  :: Vector{Float64}
    field  :: Function
end

function Base.show(io::IO, data::StreamlineData{D}) where D
    ncols = size(data.paths, 2)
    nlines = count(i -> any(isnan, @view(data.paths[:, i])), 1:ncols)
    # If the last column isn't NaN, there's one more segment
    if ncols > 0 && !any(isnan, @view(data.paths[:, ncols]))
        nlines += 1
    end
    npts = ncols - (nlines > 0 ? nlines - 1 : 0)  # subtract NaN separators
    # Format domain
    bounds = join(["[$(data.lower[i]), $(data.upper[i])]" for i in 1:D], " × ")
    print(io, "StreamlineData{$D}: $nlines streamlines, $npts points, domain $bounds")
end


"""
    ArrowData{D}

Holds arrow glyphs computed by `streamarrows`.

# Fields
- `points  :: Matrix{Float64}` — `D × N` matrix of arrow base positions.
- `vectors :: Matrix{Float64}` — `D × N` matrix of arrow direction vectors
  (unit tangent scaled by `scale`).
- `speeds  :: Vector{Float64}` — speed `‖field(p)‖` at each arrow, suitable
  for color-mapping.
- `indices :: Vector{Int}`     — column indices into the original
  `StreamlineData.paths` matrix, so per-vertex color arrays can be subsampled.
"""
struct ArrowData{D}
    points  :: Matrix{Float64}
    vectors :: Matrix{Float64}
    speeds  :: Vector{Float64}
    indices :: Vector{Int}
end
