export Sector, SectoralPrism

#= Sector (for a base shape) =#

mutable struct Sector{D} <: Shape{2,4,D}  # M = 2K
    c::SVector{2,Float64}  # center of circle
    r::Float64  # radius of circle
    ϕ₀::Float64  # center angle bisecting sector: -π ≤ ϕ₀ < π  (π excluded)
    ∆ϕ2::Float64  # "radius" in angle dimension: 0 ≤ ∆ϕ2 ≤ π (sector spans from ϕ₀ - ∆ϕ2 to ϕ₀ + ∆ϕ2)
    data::D  # auxiliary data
    Sector{D}(c,r,ϕ₀,∆ϕ,data) where {D} = new(c,r,ϕ₀,∆ϕ,data)  # suppress default outer constructor
end

function Sector(c::SVector{2,<:Real}, r::Real, ϕ::Real, ∆ϕ::Real, data::D=nothing) where {D}
    r≥0 || throw(ArgumentError("r = $r must be nonnegative."))
    -2π≤∆ϕ≤2π  || throw(ArgumentError("∆ϕ = $∆ϕ must be between -2π and 2π, inclusive."))

    ϕ₀ = rem(ϕ + ∆ϕ/2, 2π, RoundNearest)  # put ϕ₀ in [-π, π)
    ∆ϕ2 = abs(∆ϕ/2)

    return Sector{D}(c, r, ϕ₀, ∆ϕ2, data)
end

Sector(c::AbstractVector{<:Real},  # center of circle
       r::Real,  # radius of circle
       ϕ::Real,  # start angle
       ∆ϕ::Real,  # width in angle dimension; can be negative (-2π ≤ ∆ϕ ≤ 2π)
       data=nothing) =
    Sector(SVector{2}(c), r, ϕ, ∆ϕ, data)

Base.:(==)(s1::Sector, s2::Sector) = s1.c==s2.c && s1.r==s2.r && s1.ϕ₀==s2.ϕ₀ && s1.∆ϕ2==s2.∆ϕ2 && s1.data==s2.data
Base.isapprox(s1::Sector, s2::Sector) = s1.c≈s2.c && s1.r≈s2.r && s1.ϕ₀≈s2.ϕ₀ && s1.∆ϕ2≈s2.∆ϕ2 && s1.data==s2.data
Base.hash(s::Sector, h::UInt) = hash(s.c, hash(s.r, hash(s.ϕ₀, hash(s.∆ϕ2, hash(s.data, hash(:Sector, h))))))

distangle(ϕ::Real, ϕ₀:: Real) = rem(ϕ-ϕ₀, 2π, RoundNearest)  # ϕ measured from ϕ₀; result within [-π, π)

function Base.in(x::SVector{2,<:Real}, s::Sector)
    d = x - s.c
    ld = norm(d)
    ϕ = ld==0 ? s.ϕ₀ : atan(d[2], d[1])  # angle to x with respect to c

    return ld ≤ s.r && abs(distangle(ϕ, s.ϕ₀)) ≤ s.∆ϕ2
end

function surfpt_nearby(x::SVector{2,<:Real}, s::Sector)
    # Basically mimic the same function for Prism, but proceeds in the (ρ,ϕ) domain.
    d = x - s.c
    ld = norm(d)

    # Calculate the closest point in the ρ dimension and outward normal direciton there.
    r2 = s.r / 2
    ρ = ld - r2  # positive if closer to arc; negative if closer to center
    d̂ = ld ≤ Base.rtoldefault(Float64) * r2  ? SVector(cos(s.ϕ₀),sin(s.ϕ₀)) : normalize(d)

    surfρ = ρ<0 ? 0.0 : s.r  # scalar: closest point to x between center and perimeter point
    noutρ = copysign(1.0,ρ) * d̂  # SVector{2}: outward direction normal at surfρ

    absρ = abs(ρ)
    abs∆ρ = abs(r2 - absρ)  # radial distance between x and either center or perimeter, whichever closer to x

    onbndρ = abs∆ρ ≤ Base.rtoldefault(Float64) * r2  # basically r2 ≈ ρ but faster
    isoutρ = (r2 < absρ) || onbndρ

    # Calculate the closest point in the ϕ dimension and outward normal direciton there.
    ϕ = ld==0 ? 0.0 : distangle(atan(d[2], d[1]), s.ϕ₀)  # positive if closer to end side; negative if closer to start side

    ϕtemp = ϕ<0 ? s.ϕ₀ - s.∆ϕ2 : s.ϕ₀ + s.∆ϕ2
    surfϕ = SVector(cos(ϕtemp), sin(ϕtemp))  # SVector{2}: closest point in ϕ dimension
    cosθ = d̂⋅surfϕ  # cosine of angle between d and surfϕ
    ldcosθ = ld*cosθ  # ldconθ .* surfϕ is surface point
    ldsinθ = ld*sqrt(1-cosθ^2)  # always positive

    ϕtemp += copysign(π/2,ϕ)
    noutϕ = SVector(cos(ϕtemp), sin(ϕtemp))  # SVector{2}: outward direction normal at surfϕ

    absϕ = abs(ϕ)
    abs∆ϕ = abs(s.∆ϕ2 - absϕ)  # angular distance between x and closer side of sector

    onbndϕ = abs∆ϕ ≤ Base.rtoldefault(Float64) * s.∆ϕ2  # basically ∆ϕ2 ≈ ϕ but faster
    isoutϕ = (s.∆ϕ2 < absϕ) || onbndϕ

    # Pick the surface point and outward direction normal depending on the location of x.
    if isoutρ && isoutϕ  # x outside in both ρ and ϕ dimensions
        surf = surfρ .* surfϕ  # one of two end points of arc
        nout = (onbndρ && onbndϕ) ? (noutρ + noutϕ) : (d - surf)
        nout = normalize(nout)
    elseif !isoutρ && isoutϕ  # x inside in ρ dimension, but outside in ϕ dimension
        (surf, nout) = (absϕ < s.∆ϕ2+π/2) ? (ldcosθ .* surfϕ, noutϕ) : (SVector(0.0,0.0), d̂)
    elseif isoutρ && !isoutϕ  # x outside in ρ dimension, but inside in ϕ dimension
        (surf, nout) = (surfρ .* d̂, noutρ)
    else  # !isoutρ && !isoutϕ: x strictly inside sector
        (surf, nout) = (abs∆ρ ≤ ldsinθ) ? (surfρ .* d̂, noutρ) : (ldcosθ .* surfϕ, noutϕ)
    end

    return surf+s.c, nout
end

translate(s::Sector{D}, ∆::SVector{2,<:Real}) where {D} = Sector{D}(s.c+∆, s.r, s.ϕ₀, s.∆ϕ2, s.data)

function bounds(s::Sector)
    # Find the minimum and maximum coordinates among the center and two ends of the arc.
    ϕ = SVector(s.ϕ₀ - s.∆ϕ2, s.ϕ₀ + s.∆ϕ2)  # [start angle, end angle]
    v = s.r .* [SVector(0.0,0.0) [cos.(ϕ) sin.(ϕ)]']  # [center, start point of arc, end point of arc]

    # Consider using the code below once https://github.com/JuliaArrays/StaticArrays.jl/issues/498
    # is resolved:
    # l = minimum(v, dims=Val(2))[:,1]
    # u = maximum(v, dims=Val(2))[:,1]

    xs = v[1,:]
    ys = v[2,:]

    xmin = minimum(xs)
    xmax = maximum(xs)

    ymin = minimum(ys)
    ymax = maximum(ys)

    # Consider the extreme points on the arc in the Cartesian directions.
    abs(distangle(0, s.ϕ₀)) ≤ s.∆ϕ2 && (xmax = max(xmax,s.r))  # sector contains +x-direction from center
    abs(distangle(π/2, s.ϕ₀)) ≤ s.∆ϕ2 && (ymax = max(ymax, s.r))  # sector contains +y-direction from center
    abs(distangle(π, s.ϕ₀)) ≤ s.∆ϕ2  && (xmin = min(xmin, -s.r))  # sector contains -x-direction from center
    abs(distangle(3π/2, s.ϕ₀)) ≤ s.∆ϕ2  && (ymin = min(ymin, -s.r))  # sector contains -y-direction from center

    return (SVector(xmin,ymin)+s.c, SVector(xmax,ymax)+s.c)
end


#= Sectoral prism =#
const SectoralPrism = Prism{Sector{Nothing}}

# Below, if we called SectoralPrism(c, ...) in the function body, it would call the inner
# constructor Prism{Sector{Nothing}}(c, ...) because SectoralPrism = Prism{Sector{Nothing}},
# which is not what we want.
# To call the outer constructor of Prism, we should call Prism(c, ...) instead of SectoralPrism(c, ...).
SectoralPrism(c::SVector{3,<:Real},
              r::Real,
              ϕ::Real,
              ∆ϕ::Real,
              h::Real=Inf,
              a::SVector{3,<:Real}=SVector(0.0,0.0,1.0),
              data=nothing) where {K} =
    (â = normalize(a); Prism(c, Sector(SVector(0.0,0.0),r,ϕ,∆ϕ), h, [orthoaxes(â)... â], data))

SectoralPrism(c::AbstractVector{<:Real},  # center of prism
              r::Real,  # radius of sectoral base
              ϕ::Real,  # start angle of sectoral base: 0 ≤ ϕₛ < 2π  (2π excluded)
              ∆ϕ::Real,  # end angle sectoral base: 0 ≤ ϕₑ-ϕₛ ≤ 2π
              h::Real=Inf,  # height of prism
              a::AbstractVector{<:Real}=[0.0,0.0,1.0],  # axis direction of prism
              data=nothing) =
    SectoralPrism(SVector{3}(c), r, ϕ, ∆ϕ, h, SVector{3}(a), data)

function bounds_ctrcut(s::SectoralPrism)
    ax = inv(s.p)  # prism axes: columns are not only unit vectors, but also orthogonal
    if ax ≈ I  # prism axes are aligned with Cartesian directions (this covers most usage)
        l, u = bounds(s.b)  # (SVector{2}, SVector{2})
        a₁₂ = ax[:,SVector(1,2)]  # SMatrix{3,2}

        return a₁₂*l, a₁₂*u
    else  # prism axes are not aligned with Cartesian directions
        b = s.b
        r = b.r

        el = Ellipsoid(SVector(0.0,0.0,0.0), SVector(r,r,0.0), ax)  # center is set at origin to return bounds with respect to prism center
        bp = boundpts(el)  # SMatrix{3,3}: boundary points; all three columns of b are free of NaN because no column of ax is aligned with Cartesian directions
        bp′ = s.p * bp  # SMatrix{3,3}: bp in prism coordinates
        bp2′ = bp′[SVector(1,2),:]  # SMatrix{2,3}: in each column of bp′, third entry is in axis dimension, so must be zero mathematically

        ϕ = atan.(bp2′[2,:], bp2′[1,:])  # SVector{3}: angles of boundary points in base plane
        ϕsym = [ϕ; ϕ.+π]  # SVector{6}: include symmetric points with respect to base center
        ϕall = [ϕsym; SVector(b.ϕ₀-b.∆ϕ2, b.ϕ₀+b.∆ϕ2)]  # SVector{8}: angles of all boundary point candidates, except base center

        bpall′ = r .* [cos.(ϕall) sin.(ϕall) @SVector(zeros(8))]'  # SMatrix{3,8}: each column is boundary point candidate in prism coordinates
        bpall = ax * bpall′  # SMatrix{3,8}: each column is boundary point candidate in external coordinates
        bpallc = [bpall SVector(0.0,0.0,0.0)]  # SMatrix{3,9}: include base center
        ind = abs.(distangle.(ϕsym, b.ϕ₀)) .≤ b.∆ϕ2  # SVector{6}: indices of angles contained in base arc
        indallc = [ind; SVector(true,true,true)]  # SVector{9}: include arc ends and base center

        xs = bpallc[1,:]
        ys = bpallc[2,:]
        zs = bpallc[3,:]

        # It might be possible to implement the following more efficiently using reduce().
        xmin = ymin = zmin = Inf
        xmax = ymax = zmax = -Inf
        for i = 1:9
            if indallc[i]
                xmin = min(xmin, xs[i])
                ymin = min(ymin, ys[i])
                zmin = min(zmin, zs[i])

                xmax = max(xmax, xs[i])
                ymax = max(ymax, ys[i])
                zmax = max(zmax, zs[i])
            end
        end

        return (SVector(xmin,ymin,zmin), SVector(xmax,ymax,zmax))
    end
end
