# Generate sample plots for README and documentation.
# Run from the project root:
#   julia --project assets/generate_plots.jl

using UniformStreamlines
using CairoMakie

const ASSETS = joinpath(@__DIR__, "..", "assets")
mkpath(ASSETS)

# ── 1. Quick Start: simple vortex ────────────────────────────────────────────

xs = LinRange(-2, 2, 200)
ys = LinRange(-2, 2, 200)
str = stream(xs, ys, (x, y) -> -y, (x, y) -> x)

fig = Figure(size=(600, 500));
ax = Axis(fig[1, 1]; aspect=DataAspect(), title="Rigid-Body Rotation")
streamlines!(ax, str; color=:teal, linewidth=2)
tightlimits!(ax)
save(joinpath(ASSETS, "quickstart.png"), fig; px_per_unit=2)
println("  ✓ quickstart.png")

# ── 2. Density Control ───────────────────────────────────────────────────────

xs2 = LinRange(-3, 3, 200)
ys2 = LinRange(-3, 3, 200)
u_vdp(x, y) = -1 - x^2 + y
v_vdp(x, y) = 1 + x - y^2

str_sparse = stream(xs2, ys2, u_vdp, v_vdp; min_density=2, max_density=4)
str_dense  = stream(xs2, ys2, u_vdp, v_vdp; min_density=5, max_density=15)

fig2 = Figure(size=(1000, 450));
ax1 = Axis(fig2[1, 1]; aspect=DataAspect(), title="Sparse (min=2, max=4)")
streamlines!(ax1, str_sparse; color=:blue, linewidth=2)
tightlimits!(ax1)
ax2 = Axis(fig2[1, 2]; aspect=DataAspect(), title="Dense (min=5, max=15)")
streamlines!(ax2, str_dense; color=:blue)
tightlimits!(ax2)
save(joinpath(ASSETS, "density_control.png"), fig2; px_per_unit=2)
println("  ✓ density_control.png")

# ── 3. Coloring by speed ─────────────────────────────────────────────────────

xs3 = LinRange(-2, 2, 200)
ys3 = LinRange(-2, 2, 200)
str_wave = stream(xs3, ys3,
                  (x, y) -> sin(π * x) * cos(π * y),
                  (x, y) -> 0.2y)
c = colorize(str_wave, :speed)

fig3 = Figure(size=(650, 500));
ax3 = Axis(fig3[1, 1]; aspect=DataAspect(), title="Coloring by Speed")
streamlines!(ax3, str_wave; color=c, colormap=:inferno, linewidth=2)
tightlimits!(ax3)
Colorbar(fig3[1, 2]; colormap=:inferno, label="Speed",
         limits=extrema(filter(!isnan, c)))
save(joinpath(ASSETS, "coloring.png"), fig3; px_per_unit=2)
println("  ✓ coloring.png")

# ── 4. Arrows ────────────────────────────────────────────────────────────────

str_saddle = stream(xs, ys, (x, y) -> x + y, (x, y) -> x - y)
c_saddle = colorize(str_saddle, :speed)

fig4 = Figure(size=(600, 500));
ax4 = Axis(fig4[1, 1]; aspect=DataAspect(), title="Arrows — Saddle Field")
streamlines!(ax4, str_saddle; color=c_saddle, colormap=:plasma,
             with_arrows=true, arrows_every=30)
tightlimits!(ax4)
save(joinpath(ASSETS, "arrows.png"), fig4; px_per_unit=2)
println("  ✓ arrows.png")

# ── 4b. Arrow Size Control (markersize) ──────────────────────────────────────

fig4b = Figure(size=(1000, 450));
ax4b_a = Axis(fig4b[1, 1]; aspect=DataAspect(), title="markersize = 8")
streamlines!(ax4b_a, str_saddle; color=c_saddle, colormap=:plasma,
             with_arrows=true, arrows_every=30, markersize=8)
tightlimits!(ax4b_a)
ax4b_b = Axis(fig4b[1, 2]; aspect=DataAspect(), title="markersize = 20")
streamlines!(ax4b_b, str_saddle; color=c_saddle, colormap=:plasma,
             with_arrows=true, arrows_every=30, markersize=20)
tightlimits!(ax4b_b)
save(joinpath(ASSETS, "arrow_sizes.png"), fig4b; px_per_unit=2)
println("  ✓ arrow_sizes.png")

# ── 5. NaN Masking ───────────────────────────────────────────────────────────

xs5 = LinRange(-3, 3, 300)
ys5 = LinRange(-3, 3, 300)
u_mask(x, y) = (x + 1)^2 + y^2 < 1 ? NaN : x + y
v_mask(x, y) = (x + 1)^2 + y^2 < 1 ? NaN : x - y

str_mask = stream(xs5, ys5, u_mask, v_mask)

θ = LinRange(0, 2π, 100)
cx = -1 .+ cos.(θ)
cy = sin.(θ)

fig5 = Figure(size=(600, 500));
ax5 = Axis(fig5[1, 1]; aspect=DataAspect(), title="NaN Masking — Circular Obstacle")
c_mask = colorize(str_mask, :speed)
streamlines!(ax5, str_mask; color=c_mask, colormap=:coolwarm, linewidth=2)
poly!(ax5, Point2f.(cx, cy); color=(:gray, 0.3), strokecolor=:black, strokewidth=2)
tightlimits!(ax5)
save(joinpath(ASSETS, "nan_masking.png"), fig5; px_per_unit=2)
println("  ✓ nan_masking.png")

# ── 6. Seed Points ───────────────────────────────────────────────────────────

xs6 = LinRange(-3, 3, 200)
ys6 = LinRange(-3, 3, 200)
seed_x = rand(10) .* 6 .- 3  # 10 random x-coordinates in [-3, 3]
seed_y = rand(10) .* 6 .- 3  # 10 random y-coordinates in [-3, 3]
# create an Ntuple of seed 2D vectors [x,y] from the separate x and y vectors
seeds = ntuple(i -> [seed_x[i], seed_y[i]], length(seed_x))

u_vdp(x, y) = -1 - x^2 + y
v_vdp(x, y) = 1 + x - y^2

str_seed = stream(xs6, ys6,u_vdp, v_vdp;
                  seeds=seeds)

fig6 = Figure(size=(600, 500));
ax6 = Axis(fig6[1, 1]; aspect=DataAspect(), title="Seed Points — Vortex");
c_seed = colorize(str_seed, (p, v) -> atan(v[2], v[1]))
streamlines!(ax6, str_seed; color=c_seed, colormap=:hsv, linewidth=2)
scatter!(ax6, seed_x, seed_y; markersize=12, color=:black, label="seeds")
axislegend(ax6; position=:rt)
tightlimits!(ax6)
save(joinpath(ASSETS, "seeds.png"), fig6; px_per_unit=2)
println("  ✓ seeds.png")

# ── 7. Unbroken Streamlines (allow_collisions) ──────────────────────────────

str_normal = stream(xs6, ys6,u_vdp, v_vdp)  # default: allow_collisions=false, so streamlines truncate at collisions
str_unbroken = stream(xs6, ys6,u_vdp, v_vdp;
                      allow_collisions=true)

fig7 = Figure(size=(1000, 450));
c_normal = colorize(str_normal, :speed)
c_unbroken = colorize(str_unbroken, :speed)
ax7a = Axis(fig7[1, 1]; aspect=DataAspect(), title="Default (truncated)")
streamlines!(ax7a, str_normal; color=c_normal, colormap=:turbo, linewidth=2)
tightlimits!(ax7a)
ax7b = Axis(fig7[1, 2]; aspect=DataAspect(), title="allow_collisions = true")
streamlines!(ax7b, str_unbroken; color=c_unbroken, colormap=:turbo, linewidth=2)
tightlimits!(ax7b)
save(joinpath(ASSETS, "unbroken.png"), fig7; px_per_unit=2)
println("  ✓ unbroken.png")

# ── 8. 3-D Streamlines with Arrows ──────────────────────────────────────────

xs8 = LinRange(-2, 2, 40)
ys8 = LinRange(-2, 2, 40)
zs8 = LinRange(-2, 2, 40)

# ABC (Arnold–Beltrami–Childress) flow — a classic chaotic 3-D field
A, B, C = 1.0, √2, √3
str3 = stream(xs8, ys8, zs8,
              (x, y, z) -> A * sin(z) + C * cos(y),
              (x, y, z) -> B * sin(x) + A * cos(z),
              (x, y, z) -> C * sin(y) + B * cos(x);
              min_density=1, max_density=3)
c3 = colorize(str3, :speed)

fig8 = Figure(size=(700, 600));
ax8 = Axis3(fig8[1, 1]; title="3-D ABC Flow with Arrows",
            xlabel="x", ylabel="y", zlabel="z")
streamlines!(ax8, str3; color=c3, colormap=:magma, linewidth=2,
             with_arrows=true, arrows_every=25, markersize=0.12)
save(joinpath(ASSETS, "3d_arrows.png"), fig8; px_per_unit=2)
println("  ✓ 3d_arrows.png")

println("\nAll plots saved to $(ASSETS)/")

