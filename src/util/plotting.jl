const ϵrel = Base.rtoldefault(Float64)  # machine epsilon
const N_PLOT_SAMPLING_POINTS = 101

# Define drawshape() and drawshape!() functions.
@recipe(DrawShape) do scene
    Attributes()
    Theme()
end

# Implement drawshape!() for 2D shapes.
function Makie.plot!(ds::DrawShape{<:Tuple{Shape{2}}})
    shp = ds[1]
    contour!(ds, shp, levels=SVector(0.0); ds.attributes...)  # Makie.convert_arguments() below handles this new signature

    return ds
end

# Implement drawshape!() for 3D shapes.
function Makie.plot!(ds::DrawShape{<:Tuple{Shape{3},Tuple{Symbol,Real}}})
    shp = ds[1]
    cs = ds[2]
    contour!(ds, shp, cs, levels=SVector(0.0); ds.attributes...)  # Makie.convert_arguments() below handles this new signature

    return ds
end

# Define the new signature of contour!() used in drawshape!() for 2D shapes.
function Makie.convert_arguments(P::SurfaceLike, shp::Shape{2})
    lower, upper = bounds(shp)
    ∆ = upper - lower

    nw = 1; xs = range(lower[nw] - ϵrel*∆[nw], upper[nw] + ϵrel*∆[nw], length=N_PLOT_SAMPLING_POINTS)
    nw = 2; ys = range(lower[nw] - ϵrel*∆[nw], upper[nw] + ϵrel*∆[nw], length=N_PLOT_SAMPLING_POINTS)

    lvs = [level(@SVector([x,y]), shp) for x = xs, y = ys]

    return convert_arguments(P, xs, ys, lvs)
end

# Define the new signature of contour!() used in drawshape!() for 3D shapes.
function Makie.convert_arguments(P::SurfaceLike, shp::Shape{3},
                                 cs::Tuple{Symbol,Real})  # (:x or :y or :z, intercept): cross section spec
    ax, cept = cs  # axis normal to cross section, intercept

    ax==:x || ax==:y || ax==:z || @error "cs[1] = $(cs[1]) should be :x or :y or :z."
    nw = (ax==:x) + 2(ax==:y) + 3(ax==:z)  # nw = 1, 2, 3 for ax = :x, :y, :z
    nu, nv = mod1(nw+1,3), mod1(nw+2,3)

    # Set the unit vectors along the u-, v-, w-axes.
    û = SVector(ntuple(identity,Val(3))) .== nu
    v̂ = SVector(ntuple(identity,Val(3))) .== nv
    ŵ = SVector(ntuple(identity,Val(3))) .== nw

    lower, upper = bounds(shp)
    ∆ = upper - lower

    us = range(lower[nu] - ϵrel*∆[nu], upper[nu] + ϵrel*∆[nu], length=N_PLOT_SAMPLING_POINTS)
    vs = range(lower[nv] - ϵrel*∆[nv], upper[nv] + ϵrel*∆[nv], length=N_PLOT_SAMPLING_POINTS)

    lvs = [level(u*û + v*v̂ + cept*ŵ, shp) for u = us, v = vs]

    return convert_arguments(P, us, vs, lvs)
end
