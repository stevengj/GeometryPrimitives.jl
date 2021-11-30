export CrossSection

mutable struct CrossSection{S<:Shape3} <: Shape2
    shp::S  # 3D shape
    p::S²Float{3,9}  # projection matrix to cross-sectional coordinates; must be orthonormal; 3rd entry of projection is along normal axis
    c::Float  # intercept on axis normal to cross section
end

CrossSection(shp::S,
             n::SReal{3},
             c::Real
             ) where {S<:Shape3} =
    (n̂ = normalize(n); CrossSection{S}(shp, [orthoaxes(n̂)... n̂]', c))

CrossSection(shp::Shape3, n::AbsVecReal, c::Real) = CrossSection(shp, SVec{3}(n), c)

function (shp::Shape3)(ax::Symbol, c::Real)
    ax==:x || ax==:y || ax==:z || @error "ax = $(ax) should be :x or :y or :z."

    ind_n̂ = (ax==:x) + 2(ax==:y) + 3(ax==:z)  # ind_n̂ = 1, 2, 3 for ax = :x, :y, :z
    n̂ = SVec(ntuple(identity,Val(3))) .== ind_n̂

    return CrossSection(shp, n̂, c)
end

coord3d(x::SReal{2}, s::CrossSection) = (y = SFloat{3}(x.data..., s.c); s.p' * y)

level(x::SReal{2}, s::CrossSection) = level(coord3d(x,s), s.shp)

function surfpt_nearby(x::SReal{2}, s::CrossSection)
    pt, nout = surfpt_nearby(coord3d(x,s), s.shp)

    uv = SVec(1,2)
    pt2 = (s.p * pt)[uv]
    nout2 = normalize((s.p * nout)[uv])

    return pt2, nout2
end

translate(s::CrossSection, ∆::SReal{2}) = CrossSection(translate(s.shp, coord3d(∆)), s.p, s.c)

# This does not create the tightest bounds.  For the tightest bounds, this function should
# be implemented for CrossSection{S} for each S<:Shape3.
function bounds(s::CrossSection)
    bₙ, bₚ = bounds(s.shp)
    pbₙ = s.p * bₙ
    pbₚ = s.p * bₚ

    return min.(pbₙ,pbₚ), max.(pbₙ,pbₚ)
end
