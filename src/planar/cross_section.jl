export CrossSection

mutable struct CrossSection{S<:Shape3} <: Shape2
    shp::S  # 3D shape
    p::SMatrix{3,3,Float64,9}  # projection matrix to cross-sectional coordinates; must be orthonormal; 3rd entry of projection is along normal axis
    c::Float64  # intercept on axis normal to cross section
end

CrossSection(shp::S,
             n::SVector{3,<:Real},
             c::Real
             ) where {S<:Shape3} =
    (n̂ = normalize(n); CrossSection{S}(shp, [orthoaxes(n̂)... n̂]', c))

CrossSection(shp::Shape3, n::AbstractVector{<:Real}, c::Real) = CrossSection(shp, SVector{3}(n), c)

function (shp::Shape3)(ax::Symbol, c::Real)
    ax==:x || ax==:y || ax==:z || @error "ax = $(ax) should be :x or :y or :z."

    ind_n̂ = (ax==:x) + 2(ax==:y) + 3(ax==:z)  # ind_n̂ = 1, 2, 3 for ax = :x, :y, :z
    n̂ = SVector(ntuple(identity,Val(3))) .== ind_n̂

    return CrossSection(shp, n̂, c)
end

Base.:(==)(s1::CrossSection, s2::CrossSection) = s1.shp==s2.shp && s1.p==s2.p  && s1.c==s2.c
Base.isapprox(s1::CrossSection, s2::CrossSection) = s1.shp≈s2.shp && s1.p≈s2.p  && s1.c≈s2.c
Base.hash(s::CrossSection, h::UInt) = hash(s.shp, hash(s.p, hash(s.c, hash(:CrossSection, h))))

coord3d(x::SVector{2,<:Real}, s::CrossSection) = (y = SVector{3,Float64}(x..., s.c); s.p' * y)

level(x::SVector{2,<:Real}, s::CrossSection) = level(coord3d(x,s), s.shp)

function surfpt_nearby(x::SVector{2,<:Real}, s::CrossSection)
    @error "surfpt_nearby(x,s) is not supported for s::CrossSection."
end

translate(s::CrossSection, ∆::SVector{2,<:Real}) = CrossSection(translate(s.shp, coord3d(∆,s)), s.p, s.c)

# This does not create the tightest bounds.  For the tightest bounds, this function should
# be implemented for CrossSection{S} for each S<:Shape3.
function bounds(s::CrossSection)
    bₙ, bₚ = bounds(s.shp)
    pbₙ = s.p * bₙ
    pbₚ = s.p * bₚ

    return min.(pbₙ,pbₚ), max.(pbₙ,pbₚ)
end
