using CairoMakie, UniformStreamlines

dots = [
    (0.0, 1, RGBf(0.22, 0.596, 0.149)),   # green
    (cosd(30), -sind(30), RGBf(0.584, 0.345, 0.698)),   # purple
    (cosd(150), -sind(150), RGBf(0.796, 0.235, 0.2)),   # red
]

xc  = [d[1] for d in dots]
yc  = [d[2] for d in dots]
colors = [d[3] for d in dots]

x = y = collect(-3:0.01:3)
X = [i for j in y, i in x]
Y = [j for j in y, i in x]
Z = @. X + im * Y

V_inf = 1.0
a     = 3/4
α     = -π/2
zc    = xc .+ im .* yc
w     = V_inf * exp(-im*α) .* ones(size(Z))
Γ     = zeros(length(zc))

for i in eachindex(zc)
    z0 = zc[i]
    inside = abs.(Z .- z0) .< a
    Z[inside] .= NaN
    w .= w .- V_inf * a^2 * exp(im*α) ./ (Z .- z0).^2
    w .= w .- im*Γ[i] ./ (2π*(Z .- z0))
end

U = collect(real.(w)')
V = collect(-imag.(w)')

# Compute streamlines
str = evenstream(x, y, U, V; min_density=4, max_density=8)

# Color by pressure: Cp = 1 - (u² + v²)
c = colorize(str, (p, v) -> 1 - (v[1]^2 + v[2]^2))


# Build figure
fig = Figure(size=(800, 800), backgroundcolor=:transparent);
ax = Axis(fig[1, 1]; aspect=DataAspect(),
          limits=(extrema(x), extrema(y)))
hidedecorations!(ax)
hidespines!(ax)

# Streamlines with pressure coloring and arrows
streamlines!(ax, str; color=c, colorrange=(-1, 1), colormap=:jet, linewidth=3,
             with_arrows=true, markersize=20)

# Draw the Julia-colored cylinders
t = LinRange(0, 2π, 100)
for i in 1:3
    poly!(ax, Point2f.(a * cos.(t) .+ xc[i], a * sin.(t) .+ yc[i]);
          color=colors[i], strokecolor=:black, strokewidth=1)
end

tightlimits!(ax)
save(joinpath(@__DIR__, "..", "docs", "src", "assets", "logo.png"), fig; px_per_unit=2)

