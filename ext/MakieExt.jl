module MakieExt

using UniformStreamlines
using LinearAlgebra
import UniformStreamlines: streamlines, streamlines!
using Makie

Makie.convert_arguments(P::PointBased,str::StreamlineData{2}) = convert_arguments(P,str.paths[1,:], str.paths[2,:])
Makie.convert_arguments(P::PointBased,str::StreamlineData{3}) = convert_arguments(P,str.paths[1,:], str.paths[2,:], str.paths[3,:])

Makie.convert_arguments(::Makie.ArrowLike, arr::ArrowData{2}) = convert_arguments(PointBased(), 
                                                                                vec(arr.points[1,:]), 
                                                                                vec(arr.points[2,:]), 
                                                                                vec(arr.vectors[1,:]), 
                                                                                vec(arr.vectors[2,:]))

Makie.convert_arguments(::Makie.ArrowLike, arr::ArrowData{3}) = convert_arguments(PointBased(), 
                                                                                vec(arr.points[1,:]), 
                                                                                vec(arr.points[2,:]), 
                                                                                vec(arr.points[3,:]), 
                                                                                vec(arr.vectors[1,:]), 
                                                                                vec(arr.vectors[2,:]), 
                                                                                vec(arr.vectors[3,:]))


@recipe Streamlines (object,) begin
    with_arrows = false
    arrows_every = 10
    markersize = @inherit markersize
    color = :blue
    linewidth = @inherit linewidth
    linestyle = @inherit linestyle
    linecap   = @inherit linecap
    joinstyle = @inherit joinstyle
    miter_limit = @inherit miter_limit
    Makie.mixin_colormap_attributes()...
    Makie.mixin_generic_plot_attributes()...
end

Makie.args_preferred_axis(::Type{<:Streamlines}, obj::StreamlineData{2}) = Axis
Makie.args_preferred_axis(::Type{<:Streamlines}, obj::StreamlineData{3}) = Axis3

function Makie.plot!(plot::Streamlines{<:Tuple{StreamlineData{2}}})
    str = plot[:object][]  

    lines!(plot,plot.attributes,str)
    if plot[:with_arrows][]
        arr = streamarrows(str; every=plot[:arrows_every][])
        ac = arrow_color(plot[:color][], arr)
        plot_arrowheads!(plot, arr; color=ac, markersize=plot[:markersize][])
    end 

    return plot
end

function arrow_color(color, arr::ArrowData)
    color isa AbstractVector ? color[arr.indices] : color
end

function plot_arrowheads!(plot, arr::ArrowData{2}; color=:blue, markersize=12)
    pts = Point2f.(arr.points[1, :], arr.points[2, :])
    vecs = Vec2f.(arr.vectors[1, :], arr.vectors[2, :])

    rots = atan.(getindex.(vecs, 2), getindex.(vecs, 1)) .- π/2

    scatter!(
        plot,
        plot.attributes,
        pts;
        marker = :utriangle,
        rotation = rots,
        markersize = markersize,
        color = color
    )
end


function plot_arrowheads3d!(plot, arr::ArrowData{3}; markersize=0.08)
    pts = Point3f.(arr.points[1, :], arr.points[2, :], arr.points[3, :])
    vecs = Vec3f.(arr.vectors[1, :], arr.vectors[2, :], arr.vectors[3, :])
    vecs = normalize.(vecs)

    cone_marker = Makie.Tessellation(
        Makie.Cone(Makie.Point3f(0), Makie.Point3f(0, 0, 1), 0.5f0),
        16
    )

    rots = rotation_z_to.(vecs)

    meshscatter!(
        plot,
        plot.attributes,
        pts;
        marker = cone_marker,
        rotation = rots,
        markersize = markersize,
        color = arr.speeds,
    )
end


function Makie.plot!(plot::Streamlines{<:Tuple{StreamlineData{3}}})
    str = plot[:object][]
    arr = streamarrows(str; every=10, scale=0.1)

    lines!(plot, plot.attributes, str)
    plot_arrowheads3d!(plot, arr)

    return plot
end

using LinearAlgebra
using Makie

function rotation_z_to(v::Vec3f)
    z = Vec3f(0, 0, 1)
    vn = normalize(v)

    c = clamp(dot(z, vn), -1f0, 1f0)

    # same direction
    if isapprox(c, 1f0; atol=1f-6)
        return Makie.Quaternionf(1, 0, 0, 0)
    end

    # opposite direction
    if isapprox(c, -1f0; atol=1f-6)
        return Makie.qrotation(Vec3f(1, 0, 0), π)
    end

    axis = normalize(cross(z, vn))
    angle = acos(c)

    return Makie.qrotation(axis, angle)
end

end 