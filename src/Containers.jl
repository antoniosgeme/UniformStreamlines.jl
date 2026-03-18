"""
    StreamlineData{D}

Holds the output of `stream`. `D` is the spatial dimension.

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


"""
    ArrowData{D}

Holds arrow glyphs computed by `streamarrows`.

# Fields
- `points  :: Matrix{Float64}` — `D × N` matrix of arrow base positions.
- `vectors :: Matrix{Float64}` — `D × N` matrix of arrow direction vectors
  (unit tangent scaled by `scale`).
- `speeds  :: Vector{Float64}` — speed `‖field(p)‖` at each arrow, suitable
  for color-mapping.
"""
struct ArrowData{D}
    points  :: Matrix{Float64}
    vectors :: Matrix{Float64}
    speeds  :: Vector{Float64}
end
