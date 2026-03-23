# UniformStreamlines.jl

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://antoniosgeme.github.io/UniformStreamlines.jl/stable)
[![Development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://antoniosgeme.github.io/UniformStreamlines.jl/dev)
[![Test workflow status](https://github.com/antoniosgeme/UniformStreamlines.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/antoniosgeme/UniformStreamlines.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/antoniosgeme/UniformStreamlines.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/antoniosgeme/UniformStreamlines.jl)
[![Docs workflow Status](https://github.com/antoniosgeme/UniformStreamlines.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/antoniosgeme/UniformStreamlines.jl/actions/workflows/Docs.yml?query=branch%3Amain)
[![BestieTemplate](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JuliaBesties/BestieTemplate.jl/main/docs/src/assets/badge.json)](https://github.com/JuliaBesties/BestieTemplate.jl)


<p align="center">
  <img src="docs/src/assets/logo.png" alt="ex1" width="400">
</p>


Evenly-spaced streamlines for 2-D, 3-D, and N-D vector fields in Julia, using the Jobard–Lefer algorithm. Works with function-defined or grid-defined velocity fields, with built-in support for Plots.jl and Makie.jl.

## Installation

```julia
using Pkg
Pkg.add("https://github.com/antoniosgeme/UniformStreamlines.jl.git")
```

## Quick Start

```julia
using UniformStreamlines

xs = LinRange(-2, 2, 200)
ys = LinRange(-2, 2, 200)

# From functions
u(x,y) = -y
v(x,y) = x
str = stream(xs, ys, u, v)

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

<p align="center">
  <img src="docs/src/assets/quickstart.png" alt="Quick Start" width="500">
</p>

## Features

### Density Control

```julia
u(x,y) = -1 - x^2 + y
v(x,y) = 1 + x - y^2
# Sparse streamlines
str = stream(xs, ys, u, v ;
             min_density=2, max_density=4)

# Dense streamlines
str = stream(xs, ys, u, v; 
             min_density=5, max_density=15)
```

![Density Control](docs/src/assets/density_control.png)

### Color-Mapping

```julia
str = stream(xs, ys, (x, y) -> sin(π*x) * cos(π*y), (x, y) -> 0.2y)
c = colorize(str, :norm)

# Makie
streamlines(str; color=colors, color=:viridis)
```

Built-in color symbols: `:norm`, `:vx`, `:vy`, `:vz`, `:x`, `:y`, `:z`, or pass any `(pos, vel) -> scalar` function.

<p align="center">
  <img src="docs/src/assets/coloring.png" alt="Coloring by Speed" width="500">
</p>

### NaN Masking

Return `NaN` to mask regions — streamlines stop at boundaries:

```julia
u(x, y) = (x+1)^2 + y^2 < 1 ? NaN : x + y
v(x, y) = (x+1)^2 + y^2 < 1 ? NaN : x - y

str = stream(xs, ys, u, v)
c = colorize(str,:norm)
streamlines(str; color=c, colormap=:coolwarm, linewidth=2)

```

<p align="center">
  <img src="docs/src/assets/nan_masking.png" alt="NaN Masking" width="500">
</p>

### Seed Points

```julia
xs = -3:3
ys = -3:3
seed_x = rand(10) .* 6 .- 3  # 10 random x-coordinates in [-3, 3]
seed_y = rand(10) .* 6 .- 3  # 10 random y-coordinates in [-3, 3]

# create an Ntuple of seed 2D vectors [x,y] from the separate x and y vectors
seeds = ntuple(i -> [seed_x[i], seed_y[i]], length(seed_x))

str = stream(xs, ys, (x, y) -> x + y, (x, y) -> x - y;
             seeds=seeds)
c_seed = colorize(str, (p, v) -> atan(v[2], v[1]))
streamlines(str; color=c_seed, colormap=:hsv, linewidth=2)
scatter!(seed_x, seed_y; markersize=12, color=:black, label="seeds")


```

<p align="center">
  <img src="docs/src/assets/seeds.png" alt="Seed Points" width="500">
</p>

### Unbroken Streamlines

```julia
str_normal = stream(xs6, ys6,u_vdp, v_vdp)  
str_unbroken = stream(xs6, ys6,u_vdp, v_vdp;
                      allow_collisions=true)

fig = Figure(size=(1000, 450));
c_normal = colorize(str_normal, :norm)
c_unbroken = colorize(str_unbroken, :norm)
ax7a = Axis(fig[1, 1]; aspect=DataAspect(), title="Default (truncated)")
streamlines!(ax7a, str_normal; color=c_normal, colormap=:turbo, linewidth=2)
tightlimits!(ax7a)
ax7b = Axis(fig[1, 2]; aspect=DataAspect(), title="allow_collisions = true")
streamlines!(ax7b, str_unbroken; color=c_unbroken, colormap=:turbo, linewidth=2)
tightlimits!(ax7b)
```

<p align="center">
  <img src="docs/src/assets/unbroken.png" alt="Unbroken Streamlines" width="800">
</p>



### 3-D

```julia
xs = LinRange(-2, 2, 50); ys = LinRange(-2, 2, 50); zs = LinRange(-2, 2, 50)
str3 = stream(xs, ys, zs, (x,y,z) -> -y, (x,y,z) -> x, (x,y,z) -> 0.3z)
```

ABC flow with arrows and speed coloring:

```julia
A, B, C = 1.0, √2, √3
str3 = stream(xs, ys, zs,
              (x,y,z) -> A*sin(z) + C*cos(y),
              (x,y,z) -> B*sin(x) + A*cos(z),
              (x,y,z) -> C*sin(y) + B*cos(x);
              min_density=2, max_density=4)
c3 = colorize(str3, :norm)
streamlines(str3; color=c3, colormap=:magma,
            with_arrows=true, arrows_every=25, markersize=0.12)
```

<p align="center">
  <img src="docs/src/assets/3d_arrows.png" alt="3-D ABC Flow" width="500">
</p>

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

