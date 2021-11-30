export PolygonalPrism

const PolygonalPrism{K,K2} = Prism{Polygon{K,K2}}

# Below, if we called PolygonalPrism(c, ...) in the function body, it would call the inner
# constructor Prism{Polygon{K,K2,Nothing}}(c, ...) because PolygonalPrism = Prism{Polygon{K,K2,Nothing}},
# which is not what we want.
# To call the outer constructor of Prism, we should call Prism(c, ...) instead of PolygonalPrism(c, ...).
PolygonalPrism(c::SReal{3},
               v::SMat{2,K,<:Real},  # 2D coordinates of base vertices in projected prism coordinates
               h::Real=Inf,
               a::SReal{3}=SVec(0.0,0.0,1.0)
               ) where {K} =
    (â = normalize(a); Prism(c, Polygon(v), h, [orthoaxes(â)... â]))

PolygonalPrism(c::AbsVecReal,  # center of prism
               v::AbsMatReal,  # vertices of base polygon
               h::Real=Inf,  # height of prism
               a::AbsVecReal=[0.0,0.0,1.0]) =  # axis direction of prism
    (K = size(v,1); PolygonalPrism(SVec{3}(c), SMat{2,K}(v), h, SVec{3}(a)))

# Return the bounds of the center cut with respect to the prism center.
function bounds_ctrcut(s::PolygonalPrism{K}) where {K}
    p = s.p'  # projection matrix to prism coordinates: rows are not only unit vectors, but also orthogonal
    v = [s.b.v; @SMatrix(zeros(1,K))]  # SMat{3,K}: 3D vectices in prism axis coordinates
    w = p * v  # SMat{3,K}: vertices in external coordinates

    return minimum(w, dims=Val(2))[:,1], maximum(w, dims=Val(2))[:,1]
end
