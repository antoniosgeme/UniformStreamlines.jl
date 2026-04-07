@testsnippet PlotsHelpers begin
    using Random
end


# ──────────────────────────────────────────────────────────────────────────────
# Plots.jl recipe (via RecipesBase)
# ──────────────────────────────────────────────────────────────────────────────

@testitem "Plots recipe — lines only" tags=[:ext] setup=[PlotsHelpers] begin
    using UniformStreamlines, RecipesBase
    Random.seed!(110)

    xs = collect(LinRange(-1, 1, 61))
    ys = collect(LinRange(-1, 1, 61))
    str = evenstream(xs, ys, (x, y) -> -y, (x, y) -> x; min_density=0.5, max_density=1.0)

    PE = Base.get_extension(UniformStreamlines, :PlotsExt)
    sl = PE.Streamlines((str,))
    recipes = RecipesBase.apply_recipe(Dict{Symbol, Any}(), sl)
    @test !isempty(recipes)
end

@testitem "Plots recipe — with arrows" tags=[:ext] setup=[PlotsHelpers] begin
    using UniformStreamlines, RecipesBase
    Random.seed!(111)

    xs = collect(LinRange(-1, 1, 61))
    ys = collect(LinRange(-1, 1, 61))
    str = evenstream(xs, ys, (x, y) -> -y, (x, y) -> x; min_density=0.5, max_density=1.0)

    PE = Base.get_extension(UniformStreamlines, :PlotsExt)
    sl = PE.Streamlines((str,))
    recipes = RecipesBase.apply_recipe(
        Dict{Symbol, Any}(:with_arrows => true, :arrows_every => 5, :markersize => 0.5),
        sl)
    @test length(recipes) >= 2  # line series + shape series
end

@testitem "Plots recipe — with line_z color mapping and arrows" tags=[:ext] setup=[PlotsHelpers] begin
    using UniformStreamlines, RecipesBase
    Random.seed!(112)

    xs = collect(LinRange(-1, 1, 61))
    ys = collect(LinRange(-1, 1, 61))
    str = evenstream(xs, ys, (x, y) -> -y, (x, y) -> x; min_density=0.5, max_density=1.0)
    c = colorize(str, :norm)

    PE = Base.get_extension(UniformStreamlines, :PlotsExt)
    sl = PE.Streamlines((str,))
    recipes = RecipesBase.apply_recipe(
        Dict{Symbol, Any}(
            :line_z => c,
            :with_arrows => true,
            :arrows_every => 10,
            :markersize => 1.0,
        ),
        sl)
    @test length(recipes) >= 2
end

@testitem "Plots recipe — default spacing arrows" tags=[:ext] setup=[PlotsHelpers] begin
    using UniformStreamlines, RecipesBase
    Random.seed!(113)

    xs = collect(LinRange(-1, 1, 61))
    ys = collect(LinRange(-1, 1, 61))
    str = evenstream(xs, ys, (x, y) -> -y, (x, y) -> x; min_density=0.5, max_density=1.0)

    PE = Base.get_extension(UniformStreamlines, :PlotsExt)
    sl = PE.Streamlines((str,))
    recipes = RecipesBase.apply_recipe(
        Dict{Symbol, Any}(:with_arrows => true, :markersize => 0.5),
        sl)
    @test length(recipes) >= 2
end
