# UniformStreamlines.jl

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://antoniosgeme.github.io/UniformStreamlines.jl/stable)
[![Development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://antoniosgeme.github.io/UniformStreamlines.jl/dev)
[![Test workflow status](https://github.com/antoniosgeme/UniformStreamlines.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/antoniosgeme/UniformStreamlines.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/antoniosgeme/UniformStreamlines.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/antoniosgeme/UniformStreamlines.jl)
[![Docs workflow Status](https://github.com/antoniosgeme/UniformStreamlines.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/antoniosgeme/UniformStreamlines.jl/actions/workflows/Docs.yml?query=branch%3Amain)
[![BestieTemplate](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JuliaBesties/BestieTemplate.jl/main/docs/src/assets/badge.json)](https://github.com/JuliaBesties/BestieTemplate.jl)

Evenly-spaced streamlines for 2-D, 3-D, and N-D vector fields in Julia, using the Jobard–Lefer algorithm. Works with function-defined or grid-defined velocity fields, with built-in support for Plots.jl and Makie.jl.

## Installation

```julia
using Pkg
Pkg.add("UniformStreamlines")
```

## Quick Start

```julia
using UniformStreamlines

xs = LinRange(-2, 2, 200)
ys = LinRange(-2, 2, 200)

# From functions
str = stream(xs, ys, (x, y) -> -y, (x, y) -> x)

# From matrices
U = [-y for x in xs, y in ys]
V = [ x for x in xs, y in ys]
str = stream(xs, ys, U, V)
```

Plot with Plots.jl:

```julia
using Plots
streamlines(str)
```

Or with Makie:

```julia
using CairoMakie
streamlines(str)
```

## Features

### Density Control

```julia
# Sparse streamlines
str = stream(xs, ys, (x, y) -> -1 - x^2 + y, (x, y) -> 1 + x - y^2;
             min_density=2, max_density=4)

# Dense streamlines
str = stream(xs, ys, (x, y) -> -1 - x^2 + y, (x, y) -> 1 + x - y^2;
             min_density=5, max_density=15)
```

### Color-Mapping

```julia
str = stream(xs, ys, (x, y) -> sin(π*x) * cos(π*y), (x, y) -> 0.2y)
c = colorize(str, :speed)

# Plots.jl
streamlines(str; line_z=c, color=:viridis)
```

Built-in color symbols: `:speed`, `:vx`, `:vy`, `:vz`, `:x`, `:y`, `:z`, or pass any `(pos, vel) -> scalar` function.

### Directional Arrows

```julia
streamlines(str; with_arrows=true, arrows_every=20)
```

### NaN Masking

Return `NaN` to mask regions — streamlines stop at boundaries:

```julia
u(x, y) = (x+1)^2 + y^2 < 1 ? NaN : x + y
v(x, y) = (x+1)^2 + y^2 < 1 ? NaN : x - y
str = stream(xs, ys, u, v)
```

### Seed Points

```julia
str = stream(xs, ys, (x, y) -> x + y, (x, y) -> x - y;
             seeds=([-1.0, 0.0, 1.0], [0.0, 0.0, 0.0]))
```

### Unbroken Streamlines

```julia
str = stream(xs, ys, (x, y) -> -y / (x^2 + y^2 + 0.1),
                     (x, y) ->  x / (x^2 + y^2 + 0.1);
             allow_collisions=true)
```

### 3-D

```julia
xs = LinRange(-2, 2, 50); ys = LinRange(-2, 2, 50); zs = LinRange(-2, 2, 50)
str3 = stream(xs, ys, zs, (x,y,z) -> -y, (x,y,z) -> x, (x,y,z) -> 0.3z)
```

### N-D (Tuple Form)

```julia
axs = ntuple(_ -> LinRange(-2, 2, 50), 4)
fns = ((x,y,z,t) -> -y, (x,y,z,t) -> x, (x,y,z,t) -> z, (x,y,z,t) -> -t)
str4 = stream(axs, fns)
```

## API

| Function | Description |
|:---------|:------------|
| `stream` | Compute evenly-spaced streamlines |
| `colorize` | Compute per-point color values |
| `streamarrows` | Extract arrow glyphs for visualization |
| `streamlines` / `streamlines!` | Plot recipe (Plots.jl or Makie) |

See the [documentation](https://antoniosgeme.github.io/UniformStreamlines.jl/stable) for full details.

