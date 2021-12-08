const EPS_REL = τᵣ₀  # machine epsilon

# Define drawshape() and drawshape!() functions.
# hres and vres are the number of sampling points in the corresponding Cartesion derctions.
# They should be at least 2 because range(start, end, length) with start≠end requires length
# ≥ 2.  The default values of vres and hres are intentionally invalid.
@recipe(DrawShape) do scene
    Attributes(res=101, vres=1, hres=1)
end

# Implement drawshape!() for 2D shapes.
function Makie.plot!(ds::DrawShape{<:Tuple{Shape2}})
    res₀ = ds.res.val  # fields found by examining typeof(ds) and fieldnames(typeof(ds)), etc
    hres = ds.hres.val
    vres = ds.vres.val

    # If hres and vres are set to valid values, use them; otherwise use res₀.
    hres≤1 && (hres = res₀)
    vres≤1 && (vres = res₀)

    hres>1 || @error "hres = $hres should be at least 2."
    vres>1 || @error "vres = $vres should be at least 2."

    shp = ds[1]
    res = (hres, vres)

    # Makie.convert_arguments() defined below handles this new signature.
    contour!(ds, shp, res, levels=SVec(0.0); ds.attributes...)

    return ds
end

# Define the new signature of contour!() used in drawshape!() for 2D shapes.
function Makie.convert_arguments(P::SurfaceLike, shp::Shape2, res::Tuple2{Integer})
    lower, upper = bounds(shp)
    ∆ = upper - lower

    ϵrel = EPS_REL
    nw = 1; xs = range(lower[nw] - ϵrel*∆[nw], upper[nw] + ϵrel*∆[nw], length=res[nw])
    nw = 2; ys = range(lower[nw] - ϵrel*∆[nw], upper[nw] + ϵrel*∆[nw], length=res[nw])

    lvs = [level(SVec(x,y), shp) for x = xs, y = ys]

    return convert_arguments(P, xs, ys, lvs)
end
