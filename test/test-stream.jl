@testsnippet StreamHelpers begin
    using Random
    using LinearAlgebra: norm

    function split_streams(paths::AbstractMatrix{<:Real})
        segments = Matrix{Float64}[]
        start_col = 1
        ncols = size(paths, 2)
        for i in 1:ncols
            if any(isnan, @view(paths[:, i]))
                if i > start_col
                    push!(segments, Matrix{Float64}(paths[:, start_col:i-1]))
                end
                start_col = i + 1
            end
        end
        if start_col <= ncols
            push!(segments, Matrix{Float64}(paths[:, start_col:ncols]))
        end
        return segments
    end
end


# ──────────────────────────────────────────────────────────────────────────────
# 2-D, flat form
# ──────────────────────────────────────────────────────────────────────────────

@testitem "2D stream — function input" tags=[:unit] setup=[StreamHelpers] begin
    using UniformStreamlines
    Random.seed!(1)

    xs = collect(LinRange(-1, 1, 61))
    ys = collect(LinRange(-1, 1, 61))
    ufn(x, y) = 1.0
    vfn(x, y) = 0.0

    data = stream(xs, ys, ufn, vfn; min_density=0.5, max_density=1.0)
    @test data isa StreamlineData{2}
    @test data.lower == [-1.0, -1.0]
    @test data.upper == [1.0, 1.0]
    @test size(data.paths, 1) == 2
    @test size(data.paths, 2) > 0
    @test any(isnan, data.paths)

    # Horizontal field → no vertical drift
    segments = split_streams(data.paths)
    @test !isempty(segments)
    for seg in segments
        if size(seg, 2) > 1
            @test maximum(abs.(diff(seg[2, :]))) <= 1e-8
        end
    end
end

@testitem "2D stream — grid input and shape checks" tags=[:unit] setup=[StreamHelpers] begin
    using UniformStreamlines
    Random.seed!(2)

    xs = collect(LinRange(-1, 1, 41))
    ys = collect(LinRange(-1, 1, 39))
    U = [-y for x in xs, y in ys]
    V = [x for x in xs, y in ys]

    data = stream(xs, ys, U, V; min_density=0.5, max_density=1.0)
    @test data isa StreamlineData{2}
    @test size(data.paths, 1) == 2

    # Interpolated field should approximate the analytic field
    p = [0.25, -0.4]
    vel = data.field(p)
    @test isapprox(vel[1], -p[2]; atol=0.1)
    @test isapprox(vel[2], p[1]; atol=0.1)

    # Mismatched array sizes should throw
    badU = zeros(length(xs), length(ys) + 1)
    @test_throws AssertionError stream(xs, ys, badU, V)
    @test_throws AssertionError stream(xs, ys, U, badU)
end


# ──────────────────────────────────────────────────────────────────────────────
# N-D tuple form
# ──────────────────────────────────────────────────────────────────────────────

@testitem "N-D stream — tuple function input" tags=[:unit] setup=[StreamHelpers] begin
    using UniformStreamlines
    Random.seed!(10)

    xs = collect(LinRange(-1, 1, 51))
    ys = collect(LinRange(-1, 1, 51))
    ufn(x, y) = -y
    vfn(x, y) = x

    data = stream((xs, ys), (ufn, vfn); min_density=0.5, max_density=1.0)
    @test data isa StreamlineData{2}
    @test data.lower == [-1.0, -1.0]
    @test data.upper == [1.0, 1.0]
    @test size(data.paths, 1) == 2
    @test size(data.paths, 2) > 0
end

@testitem "N-D stream — tuple grid input" tags=[:unit] setup=[StreamHelpers] begin
    using UniformStreamlines
    Random.seed!(11)

    xs = collect(LinRange(-1, 1, 41))
    ys = collect(LinRange(-1, 1, 39))
    U = [-y for x in xs, y in ys]
    V = [x for x in xs, y in ys]

    data = stream((xs, ys), (U, V); min_density=0.5, max_density=1.0)
    @test data isa StreamlineData{2}
    @test size(data.paths, 1) == 2

    p = [0.25, -0.4]
    vel = data.field(p)
    @test isapprox(vel[1], -p[2]; atol=0.1)
    @test isapprox(vel[2], p[1]; atol=0.1)

    # N-D grid form throws ArgumentError for mismatched sizes
    badU = zeros(length(xs), length(ys) + 1)
    @test_throws ArgumentError stream((xs, ys), (badU, V))
    @test_throws ArgumentError stream((xs, ys), (U, badU))
end


# ──────────────────────────────────────────────────────────────────────────────
# Seeds
# ──────────────────────────────────────────────────────────────────────────────

@testitem "Explicit seeds appear in paths" tags=[:unit] setup=[StreamHelpers] begin
    using UniformStreamlines
    Random.seed!(3)

    xs = collect(LinRange(-1, 1, 51))
    ys = collect(LinRange(-1, 1, 51))
    ufn(x, y) = -y
    vfn(x, y) = x
    seeds = ([0.65, 0.0],)

    data = stream(xs, ys, ufn, vfn; seeds=seeds, min_density=0.5, max_density=1.0)
    valid_cols = [j for j in 1:size(data.paths, 2) if !any(isnan, @view(data.paths[:, j]))]
    @test !isempty(valid_cols)

    for s in seeds
        found = any(norm(data.paths[:, j] .- s) <= 1e-10 for j in valid_cols)
        @test found
    end
end

@testitem "Seeds work with N-D tuple form" tags=[:unit] setup=[StreamHelpers] begin
    using UniformStreamlines
    Random.seed!(30)

    xs = collect(LinRange(-1, 1, 51))
    ys = collect(LinRange(-1, 1, 51))
    ufn(x, y) = -y
    vfn(x, y) = x
    seeds = ([0.65, 0.0],)

    data = stream((xs, ys), (ufn, vfn); seeds=seeds, min_density=0.5, max_density=1.0)
    valid_cols = [j for j in 1:size(data.paths, 2) if !any(isnan, @view(data.paths[:, j]))]
    @test !isempty(valid_cols)

    for s in seeds
        found = any(norm(data.paths[:, j] .- s) <= 1e-10 for j in valid_cols)
        @test found
    end
end

@testitem "Seeds work with 3D" tags=[:unit] setup=[StreamHelpers] begin
    using UniformStreamlines
    Random.seed!(31)

    xs = collect(LinRange(-1, 1, 21))
    ys = collect(LinRange(-1, 1, 21))
    zs = collect(LinRange(-1, 1, 17))
    ufn(x, y, z) = 1.0
    vfn(x, y, z) = 0.0
    wfn(x, y, z) = 0.0
    seeds = ([0.0, 0.0, 0.0],)

    data = stream(xs, ys, zs, ufn, vfn, wfn; seeds=seeds, min_density=0.4, max_density=0.8)
    valid_cols = [j for j in 1:size(data.paths, 2) if !any(isnan, @view(data.paths[:, j]))]
    @test !isempty(valid_cols)

    for s in seeds
        found = any(norm(data.paths[:, j] .- s) <= 1e-10 for j in valid_cols)
        @test found
    end
end


# ──────────────────────────────────────────────────────────────────────────────
# colorize
# ──────────────────────────────────────────────────────────────────────────────

@testitem "colorize output contract" tags=[:unit] setup=[StreamHelpers] begin
    using UniformStreamlines
    Random.seed!(4)

    xs = collect(LinRange(-1, 1, 61))
    ys = collect(LinRange(-1, 1, 61))
    ufn(x, y) = -y
    vfn(x, y) = x

    data = stream(xs, ys, ufn, vfn; min_density=0.5, max_density=1.0)

    # :norm
    speed = colorize(data, :norm)
    @test length(speed) == size(data.paths, 2)
    nan_mask = [any(isnan, @view(data.paths[:, i])) for i in 1:size(data.paths, 2)]
    @test all(isnan, speed[nan_mask])
    @test all(isfinite, speed[.!nan_mask])
    @test all(>=(0.0), speed[.!nan_mask])

    # :x should match the first coordinate
    xvals = colorize(data, :x)
    valid = .!nan_mask
    @test all(isapprox.(xvals[valid], data.paths[1, valid]; atol=1e-10))

    # Custom function
    custom = colorize(data, (p, vel) -> p[1] + 2 * p[2] + vel[1])
    @test length(custom) == size(data.paths, 2)

    # Invalid symbol
    @test_throws ArgumentError colorize(data, :magnitude)
end


# ──────────────────────────────────────────────────────────────────────────────
# streamarrows
# ──────────────────────────────────────────────────────────────────────────────

@testitem "streamarrows output contract" tags=[:unit] setup=[StreamHelpers] begin
    using UniformStreamlines
    using LinearAlgebra: norm
    Random.seed!(5)

    xs = collect(LinRange(-1, 1, 61))
    ys = collect(LinRange(-1, 1, 61))
    ufn(x, y) = -y
    vfn(x, y) = x

    data = stream(xs, ys, ufn, vfn; min_density=0.8, max_density=1.2)

    arrows = streamarrows(data; every=5, scale=0.2)
    @test arrows isa ArrowData{2}
    @test size(arrows.points, 1) == 2
    @test size(arrows.vectors, 1) == 2
    @test size(arrows.points, 2) == size(arrows.vectors, 2) == length(arrows.speeds)
    @test all(isfinite, arrows.speeds)
    @test all(arrows.speeds .>= 0)

    # Arrow vectors should all have magnitude == scale
    if size(arrows.vectors, 2) > 0
        mags = [norm(@view(arrows.vectors[:, j])) for j in 1:size(arrows.vectors, 2)]
        @test all(isapprox.(mags, 0.2; atol=1e-10, rtol=1e-8))
    end

    # Huge `every` → no arrows
    no_arrows = streamarrows(data; every=10^9, scale=0.2)
    @test size(no_arrows.points) == (2, 0)
    @test size(no_arrows.vectors) == (2, 0)
    @test isempty(no_arrows.speeds)
end


# ──────────────────────────────────────────────────────────────────────────────
# 3-D
# ──────────────────────────────────────────────────────────────────────────────

@testitem "3D stream and helpers" tags=[:unit] setup=[StreamHelpers] begin
    using UniformStreamlines
    using LinearAlgebra: norm
    Random.seed!(6)

    xs = collect(LinRange(-1, 1, 25))
    ys = collect(LinRange(-1, 1, 21))
    zs = collect(LinRange(-1, 1, 17))

    ufn(x, y, z) = 1.0
    vfn(x, y, z) = 0.0
    wfn(x, y, z) = z

    # 3-D flat form, functions
    data = stream(xs, ys, zs, ufn, vfn, wfn; min_density=0.4, max_density=0.8)
    @test data isa StreamlineData{3}
    @test size(data.paths, 1) == 3

    # colorize :z
    zvals = colorize(data, :z)
    @test length(zvals) == size(data.paths, 2)

    # streamarrows 3D
    arrows = streamarrows(data; every=4, scale=0.1)
    @test arrows isa ArrowData{3}
    @test size(arrows.points, 1) == 3
    @test size(arrows.vectors, 1) == 3
    @test size(arrows.points, 2) == size(arrows.vectors, 2) == length(arrows.speeds)

    # 3-D flat form, grids
    U = [1.0 for x in xs, y in ys, z in zs]
    V = [0.0 for x in xs, y in ys, z in zs]
    W = [z  for x in xs, y in ys, z in zs]
    data_grid = stream(xs, ys, zs, U, V, W; min_density=0.4, max_density=0.8)
    @test data_grid isa StreamlineData{3}

    # Mismatched grid sizes
    badW = zeros(length(xs), length(ys), length(zs) + 1)
    @test_throws AssertionError stream(xs, ys, zs, U, V, badW)
end

@testitem "3D stream — N-D tuple form" tags=[:unit] setup=[StreamHelpers] begin
    using UniformStreamlines
    Random.seed!(60)

    xs = collect(LinRange(-1, 1, 25))
    ys = collect(LinRange(-1, 1, 21))
    zs = collect(LinRange(-1, 1, 17))

    ufn(x, y, z) = 1.0
    vfn(x, y, z) = 0.0
    wfn(x, y, z) = z

    data = stream((xs, ys, zs), (ufn, vfn, wfn); min_density=0.4, max_density=0.8)
    @test data isa StreamlineData{3}
    @test size(data.paths, 1) == 3
    @test size(data.paths, 2) > 0

    # Grid tuple form
    U = [1.0 for x in xs, y in ys, z in zs]
    V = [0.0 for x in xs, y in ys, z in zs]
    W = [z  for x in xs, y in ys, z in zs]
    data_grid = stream((xs, ys, zs), (U, V, W); min_density=0.4, max_density=0.8)
    @test data_grid isa StreamlineData{3}

    # N-D grid form throws ArgumentError for mismatched sizes
    badW = zeros(length(xs), length(ys), length(zs) + 1)
    @test_throws ArgumentError stream((xs, ys, zs), (U, V, badW))
end
