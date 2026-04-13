using UniformStreamlines
using CairoMakie

struct FitzhughNagumo{T}
    ϵ::T
    s::T
    γ::T
    β::T
end

P = FitzhughNagumo(0.1, 0.0, 1.5, 0.8)

f(x, P::FitzhughNagumo) = Point2f(
    (x[1]-x[2]-x[1]^3+P.s)/P.ϵ,
    P.γ*x[1]-x[2] + P.β
)

f(x) = f(x, P)


str = evenstream(-1.5:1.5, -1.5:1.5, f, min_density = 3, max_density = 7)



fig = Figure(resolution = (1000, 800));

# Top-left: Makie streamplot with colormap
ax11 = Axis(fig[1, 1], title = "Makie streamplot");
streamplot!(ax11, f, -1.5..1.5, -1.5..1.5, colormap = :magma);

# Top-right: Makie streamplot with explicit RGBA color function
ax12 = Axis(fig[1, 2], title = "Makie streamplot (RGBAf color)");
streamplot!(ax12, f, -1.5..1.5, -1.5..1.5, color = p -> RGBAf(p..., 0.0, 1.0));

c = colorize(str, :norm)
colors = colorize(str, (p, v) -> RGBAf(v[1], v[2], 0.0, 1.0));

# Bottom-left: UniformStreamlines with numeric color values + colormap
ax21 = Axis(fig[2, 1], title = "UniformStreamlines");
streamlines!(ax21, str; color = c, colormap = :magma, with_arrows = true, arrows_spacing = 0.4);

# Bottom-right: UniformStreamlines with explicit RGBA colors
ax22 = Axis(fig[2, 2], title = "UniformStreamlines (RGBAf color)");
streamlines!(ax22, str; color = colors, with_arrows = true, arrows_spacing = 0.4);

fig
