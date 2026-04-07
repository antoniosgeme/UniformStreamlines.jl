module PlotsExt

using UniformStreamlines
using RecipesBase
import UniformStreamlines: streamlines, streamlines!
using UniformStreamlines: StreamlineData, ArrowData, streamarrows, colorize
using LinearAlgebra: norm

# ── 2-D Streamlines ──────────────────────────────────────────────────────────

@userplot Streamlines

@recipe function f(sl::Streamlines)
    str = sl.args[1]::StreamlineData{2}

    with_arrows     = pop!(plotattributes, :with_arrows, false)
    arrows_spacing  = pop!(plotattributes, :arrows_spacing, :auto)
    arrows_every    = pop!(plotattributes, :arrows_every, nothing)
    markersize      = pop!(plotattributes, :markersize, 1.0)

    has_map = haskey(plotattributes, :line_z) && plotattributes[:line_z] !== nothing
    col     = get(plotattributes, :seriescolor, :blue)

    xs = vec(str.paths[1, :])
    ys = vec(str.paths[2, :])

    # Lines first (drawn below)
    @series begin
        seriestype := :path
        label      := ""
        if has_map
            line_z := plotattributes[:line_z]
        else
            seriescolor := col
        end
        xs, ys
    end

    # Arrowheads on top (filled triangles, no shaft)
    if with_arrows
        if arrows_every !== nothing
            arr = streamarrows(str; every=arrows_every, scale=markersize)
        else
            sp = arrows_spacing === :auto ? norm(str.upper .- str.lower) / 20 : Float64(arrows_spacing)
            arr = streamarrows(str; spacing=sp, scale=markersize)
        end
        if size(arr.points, 2) > 0
            dx = vec(arr.vectors[1, :])
            dy = vec(arr.vectors[2, :])
            ax = vec(arr.points[1, :])
            ay = vec(arr.points[2, :])

            n = length(ax)
            tri_x = Float64[]
            tri_y = Float64[]
            for i in 1:n
                θ = atan(dy[i], dx[i])
                s, c = sincos(θ)
                len   = hypot(dx[i], dy[i])
                fwd   = len * 0.06   # forward length of arrowhead
                hw    = len * 0.03  # half-width of triangle base

                # tip (forward from point)
                tx = ax[i] + fwd * c
                ty = ay[i] + fwd * s
                # base-left (perpendicular)
                blx = ax[i] - hw * s
                bly = ay[i] + hw * c
                # base-right
                brx = ax[i] + hw * s
                bry = ay[i] - hw * c

                append!(tri_x, (blx, tx, brx, blx, NaN))
                append!(tri_y, (bly, ty, bry, bly, NaN))
            end

            arrow_c = if has_map
                c_all = plotattributes[:line_z]
                # one color value per triangle vertex + close + NaN = 5 entries
                vcat([fill(c_all[idx], 5) for idx in arr.indices]...)
            else
                nothing
            end

            @series begin
                seriestype := :shape
                linewidth  := 0
                label      := ""
                primary    := false
                if has_map
                    fill_z    := arrow_c
                    linecolor := :match
                else
                    fillcolor := col
                    linecolor := col
                end
                tri_x, tri_y
            end
        end
    end
end

end
