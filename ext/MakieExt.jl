module MakieExt

using UniformStreamlines
using LinearAlgebra
import UniformStreamlines: streamlines, streamlines!
using UniformStreamlines: ArrowData
using Makie

Makie.convert_arguments(P::PointBased, str::StreamlineData{2}) =
    convert_arguments(P, str.paths[1, :], str.paths[2, :])

Makie.convert_arguments(P::PointBased, str::StreamlineData{3}) =
    convert_arguments(P, str.paths[1, :], str.paths[2, :], str.paths[3, :])

Makie.convert_arguments(::Makie.ArrowLike, arr::ArrowData{2}) =
    convert_arguments(
        PointBased(),
        vec(arr.points[1, :]),
        vec(arr.points[2, :]),
        vec(arr.vectors[1, :]),
        vec(arr.vectors[2, :]),
    )

Makie.convert_arguments(::Makie.ArrowLike, arr::ArrowData{3}) =
    convert_arguments(
        PointBased(),
        vec(arr.points[1, :]),
        vec(arr.points[2, :]),
        vec(arr.points[3, :]),
        vec(arr.vectors[1, :]),
        vec(arr.vectors[2, :]),
        vec(arr.vectors[3, :]),
    )

@recipe Streamlines (object,) begin
    with_arrows = false
    arrows_spacing = Makie.automatic
    arrows_every = nothing
    markersize = @inherit markersize
    color = :blue
    linewidth = @inherit linewidth
    linestyle = @inherit linestyle
    linecap = @inherit linecap
    joinstyle = @inherit joinstyle
    miter_limit = @inherit miter_limit
    Makie.mixin_colormap_attributes()...
    Makie.mixin_generic_plot_attributes()...
end

Makie.args_preferred_axis(::Type{<:Streamlines}, obj::StreamlineData{2}) = Axis
Makie.args_preferred_axis(::Type{<:Streamlines}, obj::StreamlineData{3}) = Axis3

function arrow_color(color, arr::ArrowData)
    return color isa AbstractVector ? color[arr.indices] : color
end

function colorrange_from(color)
    if color isa AbstractVector{<:Real}
        vals = filter(!isnan, color)
        return isempty(vals) ? (0.0, 1.0) : (Float64(minimum(vals)), Float64(maximum(vals)))
    else
        return nothing
    end
end

function resolved_colorrange(plot)
    cr = plot[:colorrange][]
    if cr === Makie.automatic || cr === nothing
        return colorrange_from(plot[:color][])
    else
        return cr
    end
end

function plot_arrowheads!(plot, arr::ArrowData{2};
    color = :blue,
    colorrange = nothing,
    colormap = nothing,
    markersize = 12
)
    pts = Point2f.(arr.points[1, :], arr.points[2, :])
    vecs = Vec2f.(arr.vectors[1, :], arr.vectors[2, :])

    rots = atan.(getindex.(vecs, 2), getindex.(vecs, 1)) .- π/2

    kw = Dict{Symbol, Any}(
        :marker => :utriangle,
        :rotation => rots,
        :markersize => markersize,
        :color => color,
    )

    if colorrange !== nothing
        kw[:colorrange] = colorrange
    end
    if colormap !== nothing
        kw[:colormap] = colormap
    end

    scatter!(plot, pts; kw...)
end

function plot_arrowheads3d!(plot, arr::ArrowData{3};
    color = :blue,
    colorrange = nothing,
    colormap = nothing,
    markersize = 0.08
)
    pts = Point3f.(arr.points[1, :], arr.points[2, :], arr.points[3, :])
    vecs = Vec3f.(arr.vectors[1, :], arr.vectors[2, :], arr.vectors[3, :])
    vecs = normalize.(vecs)

    cone_marker = Makie.Tessellation(
        Makie.Cone(Makie.Point3f(0), Makie.Point3f(0, 0, 1), 0.5f0),
        16,
    )

    rots = rotation_z_to.(vecs)

    kw3 = Dict{Symbol, Any}(
        :marker => cone_marker,
        :rotation => rots,
        :markersize => markersize,
        :color => color,
    )

    if colorrange !== nothing
        kw3[:colorrange] = colorrange
    end
    if colormap !== nothing
        kw3[:colormap] = colormap
    end

    meshscatter!(plot, pts; kw3...)
end

function Makie.plot!(plot::Streamlines{<:Tuple{StreamlineData{2}}})
    str = plot[:object][]

    lines!(plot, plot.attributes, str)

    if plot[:with_arrows][]
        ae = plot[:arrows_every][]
        if ae !== nothing
            arr = streamarrows(str; every=ae)
        else
            as_val = plot[:arrows_spacing][]
            sp = as_val === Makie.automatic ? norm(str.upper .- str.lower) / 20 : Float64(as_val)
            arr = streamarrows(str; spacing=sp)
        end
        ac = arrow_color(plot[:color][], arr)
        cr = resolved_colorrange(plot)
        cm = plot[:colormap][]
        plot_arrowheads!(
            plot,
            arr;
            color = ac,
            colorrange = cr,
            colormap = cm,
            markersize = plot[:markersize][],
        )
    end

    return plot
end

function Makie.plot!(plot::Streamlines{<:Tuple{StreamlineData{3}}})
    str = plot[:object][]

    lines!(plot, plot.attributes, str)

    if plot[:with_arrows][]
        ae = plot[:arrows_every][]
        if ae !== nothing
            arr = streamarrows(str; every=ae)
        else
            as_val = plot[:arrows_spacing][]
            sp = as_val === Makie.automatic ? norm(str.upper .- str.lower) / 20 : Float64(as_val)
            arr = streamarrows(str; spacing=sp)
        end
        ac = arrow_color(plot[:color][], arr)
        cr = resolved_colorrange(plot)
        cm = plot[:colormap][]
        # For 3D, meshscatter markersize is in data units.
        # Scale to ~2% of the domain diagonal when using the inherited theme default.
        ms = plot[:markersize][]
        if ms isa Number && ms >= 1  # likely inherited pixel value from theme
            diag = norm(str.upper .- str.lower)
            ms = 0.02 * diag
        end
        plot_arrowheads3d!(
            plot,
            arr;
            color = ac,
            colorrange = cr,
            colormap = cm,
            markersize = ms,
        )
    end

    return plot
end

function rotation_z_to(v::Vec3f)
    z = Vec3f(0, 0, 1)
    vn = normalize(v)

    c = clamp(dot(z, vn), -1f0, 1f0)

    if isapprox(c, 1f0; atol = 1f-6)
        return Makie.Quaternionf(1, 0, 0, 0)
    end

    if isapprox(c, -1f0; atol = 1f-6)
        return Makie.qrotation(Vec3f(1, 0, 0), π)
    end

    axis = normalize(cross(z, vn))
    angle = acos(c)

    return Makie.qrotation(axis, angle)
end

end