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

    data = evenstream(xs, ys, ufn, vfn; min_density=0.5, max_density=1.0)
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

    data = evenstream(xs, ys, U, V; min_density=0.5, max_density=1.0)
    @test data isa StreamlineData{2}
    @test size(data.paths, 1) == 2

    # Interpolated field should approximate the analytic field
    p = [0.25, -0.4]
    vel = data.field(p)
    @test isapprox(vel[1], -p[2]; atol=0.1)
    @test isapprox(vel[2], p[1]; atol=0.1)

    # Mismatched array sizes should throw
    badU = zeros(length(xs), length(ys) + 1)
    @test_throws AssertionError evenstream(xs, ys, badU, V)
    @test_throws AssertionError evenstream(xs, ys, U, badU)
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

    data = evenstream((xs, ys), (ufn, vfn); min_density=0.5, max_density=1.0)
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

    data = evenstream((xs, ys), (U, V); min_density=0.5, max_density=1.0)
    @test data isa StreamlineData{2}
    @test size(data.paths, 1) == 2

    p = [0.25, -0.4]
    vel = data.field(p)
    @test isapprox(vel[1], -p[2]; atol=0.1)
    @test isapprox(vel[2], p[1]; atol=0.1)

    # N-D grid form throws ArgumentError for mismatched sizes
    badU = zeros(length(xs), length(ys) + 1)
    @test_throws ArgumentError evenstream((xs, ys), (badU, V))
    @test_throws ArgumentError evenstream((xs, ys), (U, badU))
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

    data = evenstream(xs, ys, ufn, vfn; seeds=seeds, min_density=0.5, max_density=1.0)
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

    data = evenstream((xs, ys), (ufn, vfn); seeds=seeds, min_density=0.5, max_density=1.0)
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

    data = evenstream(xs, ys, zs, ufn, vfn, wfn; seeds=seeds, min_density=0.4, max_density=0.8)
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

    data = evenstream(xs, ys, ufn, vfn; min_density=0.5, max_density=1.0)

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

    data = evenstream(xs, ys, ufn, vfn; min_density=0.8, max_density=1.2)

    arrows = streamarrows(data; every=5, scale=0.2)
    @test arrows isa UniformStreamlines.ArrowData{2}
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
    data = evenstream(xs, ys, zs, ufn, vfn, wfn; min_density=0.4, max_density=0.8)
    @test data isa StreamlineData{3}
    @test size(data.paths, 1) == 3

    # colorize :z
    zvals = colorize(data, :z)
    @test length(zvals) == size(data.paths, 2)

    # streamarrows 3D
    arrows = streamarrows(data; every=4, scale=0.1)
    @test arrows isa UniformStreamlines.ArrowData{3}
    @test size(arrows.points, 1) == 3
    @test size(arrows.vectors, 1) == 3
    @test size(arrows.points, 2) == size(arrows.vectors, 2) == length(arrows.speeds)

    # 3-D flat form, grids
    U = [1.0 for x in xs, y in ys, z in zs]
    V = [0.0 for x in xs, y in ys, z in zs]
    W = [z  for x in xs, y in ys, z in zs]
    data_grid = evenstream(xs, ys, zs, U, V, W; min_density=0.4, max_density=0.8)
    @test data_grid isa StreamlineData{3}

    # Mismatched grid sizes
    badW = zeros(length(xs), length(ys), length(zs) + 1)
    @test_throws AssertionError evenstream(xs, ys, zs, U, V, badW)
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

    data = evenstream((xs, ys, zs), (ufn, vfn, wfn); min_density=0.4, max_density=0.8)
    @test data isa StreamlineData{3}
    @test size(data.paths, 1) == 3
    @test size(data.paths, 2) > 0

    # Grid tuple form
    U = [1.0 for x in xs, y in ys, z in zs]
    V = [0.0 for x in xs, y in ys, z in zs]
    W = [z  for x in xs, y in ys, z in zs]
    data_grid = evenstream((xs, ys, zs), (U, V, W); min_density=0.4, max_density=0.8)
    @test data_grid isa StreamlineData{3}

    # N-D grid form throws ArgumentError for mismatched sizes
    badW = zeros(length(xs), length(ys), length(zs) + 1)
    @test_throws ArgumentError evenstream((xs, ys, zs), (U, V, badW))
end


# ──────────────────────────────────────────────────────────────────────────────
# Coverage-focused cases (Tracer internals)
# ──────────────────────────────────────────────────────────────────────────────

@testitem "allow_collisions=true exercises collision branch" tags=[:unit] setup=[StreamHelpers] begin
    using UniformStreamlines
    using LinearAlgebra: norm
    Random.seed!(70)

    xs = collect(LinRange(-1, 1, 61))
    ys = collect(LinRange(-1, 1, 61))
    ufn(x, y) = -y
    vfn(x, y) = x

    yes_col = evenstream(xs, ys, ufn, vfn; min_density=0.8, max_density=1.2, allow_collisions=true)

    @test yes_col isa StreamlineData{2}
    @test size(yes_col.paths, 1) == 2
    @test size(yes_col.paths, 2) > 0
    @test any(isnan, yes_col.paths)

    segs_yes = split_streams(yes_col.paths)
    @test !isempty(segs_yes)

    # All non-NaN points must remain inside the domain.
    valid_cols = [j for j in 1:size(yes_col.paths, 2) if !any(isnan, @view(yes_col.paths[:, j]))]
    @test !isempty(valid_cols)
    @test all(-1.0 - 1e-10 <= yes_col.paths[1, j] <= 1.0 + 1e-10 for j in valid_cols)
    @test all(-1.0 - 1e-10 <= yes_col.paths[2, j] <= 1.0 + 1e-10 for j in valid_cols)

    # Velocity should stay finite on valid points.
    probe = yes_col.paths[:, valid_cols[clamp(div(length(valid_cols), 2), 1, length(valid_cols))]]
    @test isfinite(norm(yes_col.field(probe)))
end

@testitem "explicit seeds with allow_collisions path" tags=[:unit] setup=[StreamHelpers] begin
    using UniformStreamlines
    Random.seed!(71)

    xs = collect(LinRange(-1, 1, 41))
    ys = collect(LinRange(-1, 1, 41))
    ufn(x, y) = -y
    vfn(x, y) = x
    seeds = ([0.9, 0.0], [0.0, 0.9], [-0.9, 0.0])

    data = evenstream(xs, ys, ufn, vfn;
        seeds=seeds,
        allow_collisions=true,
        min_density=0.6,
        max_density=1.0,
    )

    @test data isa StreamlineData{2}
    @test size(data.paths, 1) == 2
    @test size(data.paths, 2) > 0

    valid_cols = [j for j in 1:size(data.paths, 2) if !any(isnan, @view(data.paths[:, j]))]
    @test !isempty(valid_cols)

    for s in seeds
        found = any(norm(data.paths[:, j] .- s) <= 1e-10 for j in valid_cols)
        @test found
    end
end


# ──────────────────────────────────────────────────────────────────────────────
# Base.show for StreamlineData
# ──────────────────────────────────────────────────────────────────────────────

@testitem "StreamlineData show method — 2D" tags=[:unit] setup=[StreamHelpers] begin
    using UniformStreamlines
    Random.seed!(80)

    xs = collect(LinRange(-1, 1, 61))
    ys = collect(LinRange(-1, 1, 61))
    ufn(x, y) = -y
    vfn(x, y) = x

    data = evenstream(xs, ys, ufn, vfn; min_density=0.5, max_density=1.0)
    buf = IOBuffer()
    show(buf, data)
    s = String(take!(buf))
    @test occursin("StreamlineData{2}:", s)
    @test occursin("streamlines", s)
    @test occursin("points", s)
    @test occursin("domain", s)
    @test occursin("[-1.0, 1.0]", s)
end

@testitem "StreamlineData show method — 3D" tags=[:unit] setup=[StreamHelpers] begin
    using UniformStreamlines
    Random.seed!(81)

    xs = collect(LinRange(-1, 1, 21))
    ys = collect(LinRange(-1, 1, 21))
    zs = collect(LinRange(-1, 1, 17))
    ufn(x, y, z) = 1.0
    vfn(x, y, z) = 0.0
    wfn(x, y, z) = 0.0

    data = evenstream(xs, ys, zs, ufn, vfn, wfn; min_density=0.4, max_density=0.8)
    buf = IOBuffer()
    show(buf, data)
    s = String(take!(buf))
    @test occursin("StreamlineData{3}:", s)
    @test occursin("streamlines", s)
    @test occursin("domain", s)
end

@testitem "StreamlineData show — empty paths" tags=[:unit] setup=[StreamHelpers] begin
    using UniformStreamlines

    # Construct a StreamlineData with an empty paths matrix
    paths = Matrix{Float64}(undef, 2, 0)
    data = UniformStreamlines.StreamlineData{2}(paths, [-1.0, -1.0], [1.0, 1.0], p -> [0.0, 0.0])
    buf = IOBuffer()
    show(buf, data)
    s = String(take!(buf))
    @test occursin("StreamlineData{2}:", s)
    @test occursin("0 streamlines", s)

    # Paths that don't end with NaN separator (last column is a real point)
    paths2 = [0.0 0.5 1.0; 0.0 0.5 1.0]   # 3 points, no NaN
    data2 = UniformStreamlines.StreamlineData{2}(paths2, [-1.0, -1.0], [1.0, 1.0], p -> [0.0, 0.0])
    buf2 = IOBuffer()
    show(buf2, data2)
    s2 = String(take!(buf2))
    @test occursin("StreamlineData{2}:", s2)
    @test occursin("1 streamlines", s2)
    @test occursin("3 points", s2)
end


# ──────────────────────────────────────────────────────────────────────────────
# colorize — individual symbol coverage
# ──────────────────────────────────────────────────────────────────────────────

@testitem "colorize all symbol shortcuts — 2D" tags=[:unit] setup=[StreamHelpers] begin
    using UniformStreamlines
    using LinearAlgebra: norm
    Random.seed!(82)

    xs = collect(LinRange(-1, 1, 61))
    ys = collect(LinRange(-1, 1, 61))
    ufn(x, y) = -y
    vfn(x, y) = x
    data = evenstream(xs, ys, ufn, vfn; min_density=0.5, max_density=1.0)
    nan_mask = [any(isnan, @view(data.paths[:, i])) for i in 1:size(data.paths, 2)]
    valid = .!nan_mask

    # :speed (alias for :norm)
    sp = colorize(data, :speed)
    @test length(sp) == size(data.paths, 2)
    @test all(isfinite, sp[valid])
    @test all(>=(0.0), sp[valid])

    # :vx / :u — should equal first velocity component
    vx_vals = colorize(data, :vx)
    u_vals  = colorize(data, :u)
    @test all(isapprox.(vx_vals[valid], u_vals[valid]; atol=1e-10))

    # :vy / :v — should equal second velocity component
    vy_vals = colorize(data, :vy)
    v_vals  = colorize(data, :v)
    @test all(isapprox.(vy_vals[valid], v_vals[valid]; atol=1e-10))

    # :y — should match second coordinate
    y_vals = colorize(data, :y)
    @test all(isapprox.(y_vals[valid], data.paths[2, valid]; atol=1e-10))
end

@testitem "colorize 3D-only symbol shortcuts" tags=[:unit] setup=[StreamHelpers] begin
    using UniformStreamlines
    Random.seed!(83)

    xs = collect(LinRange(-1, 1, 21))
    ys = collect(LinRange(-1, 1, 21))
    zs = collect(LinRange(-1, 1, 17))
    ufn(x, y, z) = z
    vfn(x, y, z) = -x
    wfn(x, y, z) = y

    data = evenstream(xs, ys, zs, ufn, vfn, wfn; min_density=0.4, max_density=0.8)
    nan_mask = [any(isnan, @view(data.paths[:, i])) for i in 1:size(data.paths, 2)]
    valid = .!nan_mask

    # :vz / :w — should equal third velocity component
    vz_vals = colorize(data, :vz)
    w_vals  = colorize(data, :w)
    @test all(isapprox.(vz_vals[valid], w_vals[valid]; atol=1e-10))

    # :z — should match third coordinate
    z_vals = colorize(data, :z)
    @test all(isapprox.(z_vals[valid], data.paths[3, valid]; atol=1e-10))
end


# ──────────────────────────────────────────────────────────────────────────────
# Edge cases in stream
# ──────────────────────────────────────────────────────────────────────────────

@testitem "explicit stepsize kwarg" tags=[:unit] setup=[StreamHelpers] begin
    using UniformStreamlines
    Random.seed!(84)

    xs = collect(LinRange(-1, 1, 41))
    ys = collect(LinRange(-1, 1, 41))
    ufn(x, y) = 1.0
    vfn(x, y) = 0.0

    data = evenstream(xs, ys, ufn, vfn; min_density=0.5, max_density=1.0, stepsize=0.01)
    @test data isa StreamlineData{2}
    @test size(data.paths, 2) > 0
end

@testitem "min_length filters short streamlines" tags=[:unit] setup=[StreamHelpers] begin
    using UniformStreamlines
    Random.seed!(85)

    xs = collect(LinRange(-1, 1, 41))
    ys = collect(LinRange(-1, 1, 41))
    ufn(x, y) = -y
    vfn(x, y) = x

    # Very high min_length should filter out most/all streamlines
    data_strict = evenstream(xs, ys, ufn, vfn; min_density=0.3, max_density=0.5, min_length=100_000)
    @test size(data_strict.paths, 2) == 0 || size(data_strict.paths, 2) < 10

    # Normal min_length should keep streamlines
    data_normal = evenstream(xs, ys, ufn, vfn; min_density=0.3, max_density=0.5, min_length=2)
    @test size(data_normal.paths, 2) > size(data_strict.paths, 2)
end

@testitem "zero-velocity field produces minimal output" tags=[:unit] setup=[StreamHelpers] begin
    using UniformStreamlines
    Random.seed!(86)

    xs = collect(LinRange(-1, 1, 21))
    ys = collect(LinRange(-1, 1, 21))

    # A field that returns NaN everywhere (simulating zero velocity / stagnation)
    data = evenstream(xs, ys, (x,y) -> 0.0, (x,y) -> 0.0; min_density=0.3, max_density=0.5)
    # With zero velocity, RK2 midpoint just stays at the seed → short streamlines
    @test data isa StreamlineData{2}
end

@testitem "seeds with 3D grid data" tags=[:unit] setup=[StreamHelpers] begin
    using UniformStreamlines
    Random.seed!(87)

    xs = collect(LinRange(-1, 1, 21))
    ys = collect(LinRange(-1, 1, 21))
    zs = collect(LinRange(-1, 1, 17))
    U = [1.0 for x in xs, y in ys, z in zs]
    V = [0.0 for x in xs, y in ys, z in zs]
    W = [0.0 for x in xs, y in ys, z in zs]
    seeds = ([0.0, 0.0, 0.0],)

    data = evenstream(xs, ys, zs, U, V, W; seeds=seeds, min_density=0.4, max_density=0.8)
    @test data isa StreamlineData{3}
    @test size(data.paths, 2) > 0
end

@testitem "streamarrows with 3D and default keywords" tags=[:unit] setup=[StreamHelpers] begin
    using UniformStreamlines
    using LinearAlgebra: norm
    Random.seed!(88)

    xs = collect(LinRange(-1, 1, 25))
    ys = collect(LinRange(-1, 1, 21))
    zs = collect(LinRange(-1, 1, 17))
    ufn(x, y, z) = 1.0
    vfn(x, y, z) = 0.0
    wfn(x, y, z) = z

    data = evenstream(xs, ys, zs, ufn, vfn, wfn; min_density=0.4, max_density=0.8)
    # Use default every and scale
    arrows = streamarrows(data)
    @test arrows isa UniformStreamlines.ArrowData{3}
    @test size(arrows.points, 2) == length(arrows.speeds)
    @test length(arrows.indices) == length(arrows.speeds)
end

@testitem "high density to trigger buffer expansion" tags=[:unit] setup=[StreamHelpers] begin
    using UniformStreamlines
    Random.seed!(89)

    # Use a domain with many seeds and tiny stepsize to generate lots of points
    xs = collect(LinRange(-5, 5, 101))
    ys = collect(LinRange(-5, 5, 101))
    ufn(x, y) = -y
    vfn(x, y) = x

    # High density + small stepsize → many vertices, potentially exceeding initial buffer
    data = evenstream(xs, ys, ufn, vfn;
        min_density=8, max_density=15, stepsize=0.001)
    @test data isa StreamlineData{2}
    @test size(data.paths, 2) > 0
end

@testitem "allow_collisions with high density" tags=[:unit] setup=[StreamHelpers] begin
    using UniformStreamlines
    Random.seed!(90)

    xs = collect(LinRange(-5, 5, 101))
    ys = collect(LinRange(-5, 5, 101))
    ufn(x, y) = -y
    vfn(x, y) = x

    data = evenstream(xs, ys, ufn, vfn;
        min_density=8, max_density=15, stepsize=0.001, allow_collisions=true)
    @test data isa StreamlineData{2}
    @test size(data.paths, 2) > 0
end


# ──────────────────────────────────────────────────────────────────────────────
# Small-magnitude fields (issue: streamlines appearing as dots)
# ──────────────────────────────────────────────────────────────────────────────

@testitem "small magnitude field produces continuous streamlines" tags=[:unit] setup=[StreamHelpers] begin
    using UniformStreamlines
    using LinearAlgebra: norm
    Random.seed!(91)

    xs = collect(LinRange(-2, 2, 51))
    ys = collect(LinRange(-2, 2, 51))

    # Scale a normal rotation field by 1e-7 to simulate tiny-magnitude fields
    scale = 1e-7
    ufn(x, y) = -y * scale
    vfn(x, y) = x * scale

    data = evenstream(xs, ys, ufn, vfn; min_density=0.5, max_density=1.0)
    segments = split_streams(data.paths)

    # With normalization, streamlines should be multi-vertex (not dots)
    long_segs = filter(s -> size(s, 2) >= 5, segments)
    @test !isempty(long_segs)

    # Each long segment should span a significant fraction of the domain
    for seg in long_segs
        dx = maximum(seg[1, :]) - minimum(seg[1, :])
        dy = maximum(seg[2, :]) - minimum(seg[2, :])
        span = hypot(dx, dy)
        @test span > 0.1   # should cover at least some physical distance
    end
end

@testitem "endgrid out-of-bounds triggers trim" tags=[:unit] setup=[StreamHelpers] begin
    using UniformStreamlines
    Random.seed!(91)

    # Use seeds near domain boundary with a field pointing outward
    # to trigger the endgrid out-of-bounds check
    xs = collect(LinRange(-1, 1, 21))
    ys = collect(LinRange(-1, 1, 21))
    # Field pointing radially outward — streamlines quickly leave domain
    ufn(x, y) = x * 10.0
    vfn(x, y) = y * 10.0
    seeds = ([0.9, 0.9], [-0.9, -0.9], [0.9, -0.9], [-0.9, 0.9])

    data = evenstream(xs, ys, ufn, vfn; seeds=seeds, min_density=0.5, max_density=1.0)
    @test data isa StreamlineData{2}
end

@testitem "seeds with N-D tuple grid form" tags=[:unit] setup=[StreamHelpers] begin
    using UniformStreamlines
    Random.seed!(92)

    xs = collect(LinRange(-1, 1, 41))
    ys = collect(LinRange(-1, 1, 39))
    U = [-y for x in xs, y in ys]
    V = [x  for x in xs, y in ys]
    seeds = ([0.5, 0.0], [-0.5, 0.0])

    data = evenstream((xs, ys), (U, V); seeds=seeds, min_density=0.5, max_density=1.0)
    @test data isa StreamlineData{2}
    @test size(data.paths, 2) > 0
end


# ──────────────────────────────────────────────────────────────────────────────
# streamarrows — arc-length (uniform) mode
# ──────────────────────────────────────────────────────────────────────────────

@testitem "streamarrows spacing — uniform arc-length placement 2D" tags=[:unit] setup=[StreamHelpers] begin
    using UniformStreamlines
    using LinearAlgebra: norm
    Random.seed!(93)

    xs = collect(LinRange(-2, 2, 101))
    ys = collect(LinRange(-2, 2, 101))
    # Circular field → closed-ish streamlines with known radius
    data = evenstream(xs, ys, (x, y) -> -y, (x, y) -> x; min_density=0.8, max_density=1.5)

    arrows = streamarrows(data; spacing=0.5, scale=0.1)
    @test arrows isa UniformStreamlines.ArrowData{2}
    @test size(arrows.points, 1) == 2
    @test size(arrows.vectors, 1) == 2
    @test size(arrows.points, 2) == size(arrows.vectors, 2) == length(arrows.speeds)
    @test length(arrows.indices) == length(arrows.speeds)
    @test all(isfinite, arrows.speeds)
    @test all(arrows.speeds .>= 0)

    # All arrow vectors should have magnitude == scale
    if size(arrows.vectors, 2) > 0
        mags = [norm(@view(arrows.vectors[:, j])) for j in 1:size(arrows.vectors, 2)]
        @test all(isapprox.(mags, 0.1; atol=1e-10))
    end
end

@testitem "streamarrows spacing — 3D" tags=[:unit] setup=[StreamHelpers] begin
    using UniformStreamlines
    using LinearAlgebra: norm
    Random.seed!(94)

    xs = collect(LinRange(-1, 1, 25))
    ys = collect(LinRange(-1, 1, 21))
    zs = collect(LinRange(-1, 1, 17))
    data = evenstream(xs, ys, zs, (x,y,z)->1.0, (x,y,z)->0.0, (x,y,z)->z;
        min_density=0.4, max_density=0.8)

    arrows = streamarrows(data; spacing=0.4, scale=0.05)
    @test arrows isa UniformStreamlines.ArrowData{3}
    @test size(arrows.points, 1) == 3
    @test size(arrows.points, 2) > 0
    @test all(isfinite, arrows.speeds)
end

@testitem "streamarrows spacing — huge spacing gives no arrows" tags=[:unit] setup=[StreamHelpers] begin
    using UniformStreamlines
    Random.seed!(95)

    xs = collect(LinRange(-1, 1, 41))
    ys = collect(LinRange(-1, 1, 41))
    data = evenstream(xs, ys, (x,y)->-y, (x,y)->x; min_density=0.5, max_density=1.0)

    arrows = streamarrows(data; spacing=1000.0, scale=0.1)
    @test size(arrows.points, 2) == 0
    @test isempty(arrows.speeds)
end

@testitem "streamarrows spacing — arrows are uniformly spaced" tags=[:unit] setup=[StreamHelpers] begin
    using UniformStreamlines
    using LinearAlgebra: norm
    Random.seed!(96)

    xs = collect(LinRange(-2, 2, 201))
    ys = collect(LinRange(-2, 2, 201))
    # Uniform horizontal flow → straight streamlines with constant spacing
    data = evenstream(xs, ys, (x,y) -> 1.0, (x,y) -> 0.0; min_density=0.5, max_density=1.0)

    spacing_val = 0.3
    arrows = streamarrows(data; spacing=spacing_val, scale=0.05)

    # Group arrow points by their y-coordinate (each streamline is a horizontal line)
    if size(arrows.points, 2) > 1
        # For this uniform flow, all arrows on the same streamline share the same y
        ys_arr = arrows.points[2, :]
        # Find groups by rounding y
        unique_y = unique(round.(ys_arr; digits=6))
        for yval in unique_y
            mask = findall(y -> isapprox(y, yval; atol=1e-5), ys_arr)
            length(mask) < 2 && continue
            xs_group = sort(arrows.points[1, mask])
            diffs = diff(xs_group)
            # All inter-arrow distances should be approximately equal to spacing_val
            for d in diffs
                @test isapprox(d, spacing_val; atol=0.02)
            end
        end
    end
end
