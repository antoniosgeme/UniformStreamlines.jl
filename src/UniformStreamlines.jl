module UniformStreamlines

using FastInterpolations
using LinearAlgebra: norm
using Random: shuffle!
using StaticArrays

include("Containers.jl")
include("Tracer.jl")  
include("Stream.jl")  

export stream                                         # low-level core
export evenstream                                     # high-level entry point
export colorize, streamarrows                         # post-processing
export StreamlineData, ArrowData                      # result types

"""
    streamlines(str::StreamlineData; kwargs...)
    streamlines!(ax, str::StreamlineData; kwargs...)

Plot evenly-spaced streamlines from precomputed [`StreamlineData`](@ref).

`streamlines` creates a new figure; `streamlines!` adds to an existing axis or plot.
Both work with **Makie** (via `MakieExt`) and **Plots.jl** (via `PlotsExt`).

# Common keyword arguments

| Keyword | Default | Description |
|:--------|:--------|:------------|
| `with_arrows` | `false` | Draw directional arrowheads along streamlines |
| `arrows_every` | `10` | Place an arrowhead every N path vertices |

# Makie-specific keyword arguments

| Keyword | Default | Description |
|:--------|:--------|:------------|
| `color` | `:blue` | A color, or a per-point `Vector` from [`colorize`](@ref) |
| `colormap` | — | Colormap when `color` is numeric (e.g. `:viridis`, `:plasma`) |
| `colorrange` | automatic | `(min, max)` for the color mapping |
| `linewidth` | inherited | Streamline width |
| `markersize` | inherited (2-D) | Arrowhead size. In 2-D this is in pixels; in 3-D it is in data units (default auto-scaled to ~2 % of domain diagonal) |

# Plots.jl-specific keyword arguments

| Keyword | Default | Description |
|:--------|:--------|:------------|
| `line_z` | — | Per-point color values from [`colorize`](@ref) |
| `arrow_scale` | `1.0` | Multiplicative scale factor for arrowhead size |
| `color` | `:blue` | Series color / colormap when using `line_z` |

# Examples

```julia
using UniformStreamlines

xs = LinRange(-2, 2, 200)
ys = LinRange(-2, 2, 200)
str = stream(xs, ys, (x, y) -> -y, (x, y) -> x)

# Makie
using CairoMakie
streamlines(str; color=colorize(str, :speed), colormap=:viridis,
            with_arrows=true, arrows_every=30)

# Plots.jl
using Plots
c = colorize(str, :speed)
streamlines(str; line_z=c, color=:viridis,
            with_arrows=true, arrows_every=30)
```
"""
function streamlines end
function streamlines! end
export streamlines, streamlines!

end
