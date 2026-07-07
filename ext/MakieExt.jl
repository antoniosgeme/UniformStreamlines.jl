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
        return Makie.automatic
    end
end

function compute_arrow_data(str::StreamlineData, ae, as_val)
    if ae !== nothing
        return streamarrows(str; every = ae)
    else
        sp = as_val === Makie.automatic ? norm(str.upper .- str.lower) / 20 : Float64(as_val)
        return streamarrows(str; spacing = sp)
    end
end

function lifted_colorrange(plot)
    return @lift begin
        cr = $(plot[:colorrange])
        if cr === Makie.automatic || cr === nothing
            colorrange_from($(plot[:color]))
        else
            cr
        end
    end
end

function Makie.plot!(plot::Streamlines{<:Tuple{StreamlineData{2}}})
    lines!(plot, plot.attributes, plot[:object])

    if plot[:with_arrows][]
        arr = @lift compute_arrow_data(
            $(plot[:object]),
            $(plot[:arrows_every]),
            $(plot[:arrows_spacing]),
        )
        pts = @lift begin
            a = $arr
            Point2f.(a.points[1, :], a.points[2, :])
        end
        rots = @lift begin
            a = $arr
            vecs = Vec2f.(a.vectors[1, :], a.vectors[2, :])
            atan.(getindex.(vecs, 2), getindex.(vecs, 1)) .- π / 2
        end
        ac = @lift arrow_color($(plot[:color]), $arr)
        scatter!(
            plot,
            pts;
            marker = :utriangle,
            rotation = rots,
            markersize = plot[:markersize],
            color = ac,
            colorrange = lifted_colorrange(plot),
            colormap = plot[:colormap],
        )
    end

    return plot
end

function Makie.plot!(plot::Streamlines{<:Tuple{StreamlineData{3}}})
    lines!(plot, plot.attributes, plot[:object])

    if plot[:with_arrows][]
        arr = @lift compute_arrow_data(
            $(plot[:object]),
            $(plot[:arrows_every]),
            $(plot[:arrows_spacing]),
        )
        pts = @lift begin
            a = $arr
            Point3f.(a.points[1, :], a.points[2, :], a.points[3, :])
        end
        rots = @lift begin
            a = $arr
            vecs = normalize.(Vec3f.(a.vectors[1, :], a.vectors[2, :], a.vectors[3, :]))
            rotation_z_to.(vecs)
        end
        ac = @lift arrow_color($(plot[:color]), $arr)
        cone_marker = Makie.Tessellation(
            Makie.Cone(Makie.Point3f(0), Makie.Point3f(0, 0, 1), 0.5f0),
            16,
        )
        ms = @lift begin
            str = $(plot[:object])
            ms_val = $(plot[:markersize])
            if ms_val isa Number && ms_val >= 1
                return 0.02 * norm(str.upper .- str.lower)
            else
                return ms_val
            end
        end
        meshscatter!(
            plot,
            pts;
            marker = cone_marker,
            rotation = rots,
            markersize = ms,
            color = ac,
            colorrange = lifted_colorrange(plot),
            colormap = plot[:colormap],
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