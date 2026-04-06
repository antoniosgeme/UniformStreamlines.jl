@testsnippet MakieHelpers begin
    using Random
    using LinearAlgebra: norm
end


# ──────────────────────────────────────────────────────────────────────────────
# Makie 2-D
# ──────────────────────────────────────────────────────────────────────────────

@testitem "Makie 2D streamlines — no arrows" tags=[:ext] setup=[MakieHelpers] begin
    using UniformStreamlines, CairoMakie
    Random.seed!(100)

    xs = collect(LinRange(-1, 1, 61))
    ys = collect(LinRange(-1, 1, 61))
    str = evenstream(xs, ys, (x, y) -> -y, (x, y) -> x; min_density=0.5, max_density=1.0)

    fig = streamlines(str)
    @test fig isa Makie.FigureAxisPlot
end

@testitem "Makie 2D streamlines — with arrows and colormap" tags=[:ext] setup=[MakieHelpers] begin
    using UniformStreamlines, CairoMakie
    Random.seed!(101)

    xs = collect(LinRange(-1, 1, 61))
    ys = collect(LinRange(-1, 1, 61))
    str = evenstream(xs, ys, (x, y) -> -y, (x, y) -> x; min_density=0.5, max_density=1.0)
    c = colorize(str, :norm)

    fig = streamlines(str; color=c, colormap=:viridis, with_arrows=true, arrows_every=5)
    @test fig isa Makie.FigureAxisPlot
end

@testitem "Makie 2D streamlines! — bang variant" tags=[:ext] setup=[MakieHelpers] begin
    using UniformStreamlines, CairoMakie
    Random.seed!(102)

    xs = collect(LinRange(-1, 1, 61))
    ys = collect(LinRange(-1, 1, 61))
    str = evenstream(xs, ys, (x, y) -> -y, (x, y) -> x; min_density=0.5, max_density=1.0)

    fig = Figure()
    ax = Axis(fig[1, 1])
    p = streamlines!(ax, str; color=:red, with_arrows=true, arrows_every=10)
    @test p isa Makie.Plot
end


# ──────────────────────────────────────────────────────────────────────────────
# Makie 3-D
# ──────────────────────────────────────────────────────────────────────────────

@testitem "Makie 3D streamlines — with arrows" tags=[:ext] setup=[MakieHelpers] begin
    using UniformStreamlines, CairoMakie
    Random.seed!(103)

    xs = collect(LinRange(-1, 1, 21))
    ys = collect(LinRange(-1, 1, 21))
    zs = collect(LinRange(-1, 1, 17))
    str = evenstream(xs, ys, zs, (x, y, z) -> 1.0, (x, y, z) -> 0.0, (x, y, z) -> z;
        min_density=0.4, max_density=0.8)

    fig = streamlines(str; with_arrows=true, arrows_every=4)
    @test fig isa Makie.FigureAxisPlot
end

@testitem "Makie 3D streamlines — no arrows" tags=[:ext] setup=[MakieHelpers] begin
    using UniformStreamlines, CairoMakie
    Random.seed!(106)

    xs = collect(LinRange(-1, 1, 21))
    ys = collect(LinRange(-1, 1, 21))
    zs = collect(LinRange(-1, 1, 17))
    str = evenstream(xs, ys, zs, (x, y, z) -> 1.0, (x, y, z) -> 0.0, (x, y, z) -> z;
        min_density=0.4, max_density=0.8)

    fig = streamlines(str)
    @test fig isa Makie.FigureAxisPlot
end


# ──────────────────────────────────────────────────────────────────────────────
# convert_arguments
# ──────────────────────────────────────────────────────────────────────────────

@testitem "Makie convert_arguments — StreamlineData 2D" tags=[:ext] setup=[MakieHelpers] begin
    using UniformStreamlines, CairoMakie
    Random.seed!(104)

    xs = collect(LinRange(-1, 1, 61))
    ys = collect(LinRange(-1, 1, 61))
    str = evenstream(xs, ys, (x, y) -> -y, (x, y) -> x; min_density=0.5, max_density=1.0)

    result = Makie.convert_arguments(PointBased(), str)
    @test !isempty(result)
end

@testitem "Makie convert_arguments — StreamlineData 3D" tags=[:ext] setup=[MakieHelpers] begin
    using UniformStreamlines, CairoMakie
    Random.seed!(107)

    xs = collect(LinRange(-1, 1, 21))
    ys = collect(LinRange(-1, 1, 21))
    zs = collect(LinRange(-1, 1, 17))
    str = evenstream(xs, ys, zs, (x, y, z) -> 1.0, (x, y, z) -> 0.0, (x, y, z) -> z;
        min_density=0.4, max_density=0.8)

    result = Makie.convert_arguments(PointBased(), str)
    @test !isempty(result)
end


# ──────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ──────────────────────────────────────────────────────────────────────────────

@testitem "Makie helpers — arrow_color, colorrange_from" tags=[:ext] setup=[MakieHelpers] begin
    using UniformStreamlines, CairoMakie
    Random.seed!(105)

    xs = collect(LinRange(-1, 1, 61))
    ys = collect(LinRange(-1, 1, 61))
    str = evenstream(xs, ys, (x, y) -> -y, (x, y) -> x; min_density=0.5, max_density=1.0)

    ME = Base.get_extension(UniformStreamlines, :MakieExt)

    # arrow_color with scalar
    arr = streamarrows(str; every=5)
    @test ME.arrow_color(:blue, arr) === :blue

    # arrow_color with vector
    c = colorize(str, :norm)
    ac = ME.arrow_color(c, arr)
    @test length(ac) == size(arr.points, 2)

    # colorrange_from with vector
    cr = ME.colorrange_from(c)
    @test cr isa Tuple{Float64, Float64}
    @test cr[1] <= cr[2]

    # colorrange_from with scalar
    @test ME.colorrange_from(:blue) === nothing

    # colorrange_from with all-NaN vector
    @test ME.colorrange_from([NaN, NaN]) == (0.0, 1.0)
end

@testitem "rotation_z_to edge cases" tags=[:ext] setup=[MakieHelpers] begin
    using UniformStreamlines, CairoMakie

    ME = Base.get_extension(UniformStreamlines, :MakieExt)

    # Aligned with +z → identity quaternion
    q1 = ME.rotation_z_to(Vec3f(0, 0, 1))
    @test q1 isa Makie.Quaternion

    # Aligned with -z → 180° rotation
    q2 = ME.rotation_z_to(Vec3f(0, 0, -1))
    @test q2 isa Makie.Quaternion

    # General direction
    q3 = ME.rotation_z_to(Vec3f(1, 1, 0))
    @test q3 isa Makie.Quaternion
end
