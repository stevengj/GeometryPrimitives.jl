export Cylinder

const Cylinder = Prism{Ball{2,4}}

# Below, if we called Cylinder(c, ...) in the function body, it would call the inner
# constructor Prism{Ball{2,4,Nothing}}(c, ...) because Cylinder = Prism{Ball{2,4,Nothing}},
# which is not what we want.
# To call the outer constructor of Prism, we should call Prism(c, ...) instead of Cylinder(c, ...).
Cylinder(c::SReal{3},
         r::Real,
         h::Real=Inf,
         a::SReal{3}=SVec(0.0,0.0,1.0)) =
    (â = normalize(a); Prism(c, Ball(SVec(0.0,0.0),r), h, [orthoaxes(â)... â]))

Cylinder(c::AbsVecReal,  # center of cylinder
         r::Real,  # radius of base
         h::Real=Inf,  # height of cylinder
         a::AbsVecReal=[0.0,0.0,1.0]) =  # axis direction of cylinder
    Cylinder(SVec{3}(c), r, h, SVec{3}(a))

# Return the bounds of the center cut with respect to the prism center.
function bounds_ctrcut(s::Cylinder)
    ax = s.p'  # prism axes: columns are not only unit vectors, but also orthogonal
    r = s.b.r
    el = Ellipsoid(SVec(0.0,0.0,0.0), SVec(r,r,0.0), ax)  # center is set at origin to return bounds with respect to prism center

    return bounds(el)
end
