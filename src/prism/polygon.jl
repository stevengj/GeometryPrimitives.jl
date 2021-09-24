export Polygon, PolygonalPrism
export regpoly, isosceles

#= Polygon (for a base shape) =#

# Assume the followings for the polygon represented by Polygon:
# - The polygon is convex.
# - The vertices are listed in the counter-clockwise order around the origin.
mutable struct Polygon{K,K2,D} <: Shape{2,4,D}  # K2 = 2K
    v::SMatrix{2,K,Float64,K2}  # vertices
    n::SMatrix{2,K,Float64,K2}  # direction normals to edges
    data::D  # auxiliary data
    Polygon{K,K2,D}(v,n,data) where {K,K2,D} = new(v,n,data)  # suppress default outer constructor
end

function Polygon(v::SMatrix{2,K,<:Real}, data::D=nothing) where {K,D}
    # Sort the vertices in the counter-clockwise direction
    w = v .- mean(v, dims=Val(2))  # v in center-of-mass coordinates
    ϕ = mod.(atan.(w[2,:], w[1,:]), 2π)  # SVector{K}: angle of vertices between 0 and 2π; `%` does not work for negative angle
    if !issorted(ϕ)
        # Do this only when ϕ is not sorted, because the following uses allocations.
        ind = MVector{K}(sortperm(ϕ))  # sortperm(::SVector) currently returns Vector, not MVector
        v = v[:,ind]  # SVector{K}: sorted v
    end

    # Calculate the increases in angle between neighboring edges.
    ∆v = hcat(diff(v, dims=Val(2)), SMatrix{2,1}(v[:,1]-v[:,end]))  # SMatrix{2,K}: edge directions
    ∆z = ∆v[1,:] + im * ∆v[2,:]  # SVector{K}: edge directions as complex numbers
    icurr = ntuple(identity, Val(K-1))
    inext = ntuple(x->x+1, Val(K-1))
    ∆ϕ = angle.(∆z[SVector(inext)] ./ ∆z[SVector(icurr)])  # angle returns value between -π and π

    # Check all the angle increases are positive.  If they aren't, the polygon is not convex.
    all(∆ϕ .> 0) || throw("v = $v should represent vertices of convex polygon.")

    n = [∆v[2,:] -∆v[1,:]]'  # SMatrix{2,K}; outward normal directions to edges
    n = n ./ hypot.(n[1,:], n[2,:])'  # normalize

    return Polygon{K,2K,D}(v,n,data)
end

Polygon(v::AbstractMatrix{<:Real}, data=nothing) = (K = size(v,2); Polygon(SMatrix{2,K}(v), data))

Base.:(==)(s1::Polygon, s2::Polygon) = s1.v==s2.v && s1.n==s2.n && s1.data==s2.data  # assume sorted v
Base.isapprox(s1::Polygon, s2::Polygon) = s1.v≈s2.v && s1.n≈s2.n && s1.data==s2.data  # assume sorted v
Base.hash(s::Polygon, h::UInt) = hash(s.v, hash(s.n, hash(s.data, hash(:Polygon, h))))

function level(x::SVector{2,<:Real}, s::Polygon)
    c = mean(s.v, dims=Val(2))  # center of mass
    
    d = sum(s.n .* (x .- c), dims=Val(1))
    r = sum(s.n .* (s.v .- c), dims=Val(1))
    @assert all(r .> 0)

    return maximum(d ./ r) - 1.0
end

function surfpt_nearby(x::SVector{2,<:Real}, s::Polygon{K}) where {K}
    # Calculate the signed distances from x to edge lines.
    ∆xe = sum(s.n .* (x .- s.v), dims=Val(1))[1,:]  # SVector{K}: values of equations of edge lines
    abs∆xe = abs.(∆xe)  # SVector{K}

    # Determine if x is outside of edges, inclusive.
    sz = abs.((-)(bounds(s)...))  # SVector{2}
    onbnd = abs∆xe .≤ Base.rtoldefault(Float64) * max(sz.data...)  # SVector{K}
    isout = (∆xe.>0) .| onbnd  # SVector{K}

    # For x inside the polygon, it is easy to find the closest surface point: we can simply
    # pick the closest edge and find its point closest to x.
    # For x outside the polygon, there are many cases depending on how many edges x lies
    # outside of.  However, x that is sufficiently close to the polygon should be outside of
    # at most two edges.  (Outside two edges near the vertices, and one edge elsewhere.)
    # Therefore, we will consider only cout ≤ 2 below.
    cout = count(isout)
    if cout == 2  # x is outside two edges
        # We could choose to find ind corresponding to the two nonnegative ∆xe, but such an
        # operation leads to allocations because it does not create an SVector.  Instead,
        # find the closest vertex directly, which is the surface point we are looking for
        # for cout = 2.
        ∆xv = x .- s.v
        l∆xv = hypot.(∆xv[1,:], ∆xv[2,:])
        imin = argmin(l∆xv)
        surf = s.v[:,imin]
        imin₋₁ = mod1(imin-1,K)

        if onbnd[imin] && onbnd[imin₋₁]  # x is very close to vertex imin
            nout = s.n[:,imin] + s.n[:,imin₋₁]
        else
            nout = x - s.v[:,imin]
        end
        nout = normalize(nout)
    else  # cout ≤ 1 or cout ≥ 3
        cout ≤ 1 || @warn "x = $x is outside $cout edges: too far from polygon with vertices $(s.v); " *
                            "result could be inaccurate."
        # Choose the closest edge to x.
        # If cout = 0, all ∆xe are negative, so the largest ∆xe is the smallest in magnitude
        # and corresponds to the closest edge.
        # If cout = 1, all but one ∆xe are negative, so again the largest ∆xe is the only
        # nonnegative one and corresponds to the edge outside which x lies.
        # Even for cout = 3, this seems to be a reasonable choice from a simple geometric
        # argument.
        imax = argmax(∆xe)
        vmax, nmax = s.v[:,imax], s.n[:,imax]

        ∆x = (nmax⋅(vmax-x)) .* nmax
        surf = x + ∆x
        nout = nmax
    end

    return surf, nout
end

translate(s::Polygon{K,K2,D}, ∆::SVector{2,<:Real}) where {K,K2,D} = Polygon{K,K2,D}(s.v .+ ∆, s.n, s.data)

function bounds(s::Polygon)
    l = minimum(s.v, dims=Val(2))[:,1]
    u = maximum(s.v, dims=Val(2))[:,1]

    return (l, u)
end

#= Factory methods =#
# Regular polygon
function regpoly(::Val{K},  # number of vertices
                 r::Real,  # distance between center and each vertex
                 θ::Real=π/2,  # angle from +x-direction towards first vertex; π/2 corresponds to +y-direction
                 c::SVector{2,<:Real}=SVector(0.0,0.0),  # center location
                 data=nothing) where {K}
    ∆θ = 2π / K

    θs = θ .+ ∆θ .* SVector(ntuple(k->k-1, Val(K)))  # SVector{K}: angles of vertices
    v = c .+ r .* [cos.(θs) sin.(θs)]'  # SMatrix{2,K}: locations of vertices

    return Polygon(v, data)
end

regpoly(k::Integer,  # number of vertices
        r::Real,  # radius: distance from center to vertices
        θ::Real=π/2,  # angle of first vertex
        c::AbstractVector{<:Real}=[0.0,0.0],  # [x, y]: center of regular polygon
        data=nothing) =
   regpoly(Val(k), r, θ, SVector{2}(c), data)

# Isosceles triangle
function Isosceles(base::NTuple{2,SVector{2,<:Real}},
                   h::Real,
                   data=nothing)
    m = (base[1] + base[2]) / 2  # midpoint of base
    bvec = normalize(base[2] - base[1])  # unit direction of base
    hvec = [-bvec[2], bvec[1]]  # unit direction of height
    p = m + h.*hvec  # apex

    v = [base[1] base[2] p]  # vertices

    return Polygon(v, data)
end

Isosceles(base::NTuple{2,AbstractVector{<:Real}},  # (end point 1, end point 2): two end points of base
          h::Real,  # height drawn normal to base; direction is such that base pt 1, base pt 2, apex are put in counter-clockwise order
          data=nothing) = Isosceles(SVector{2}.(base), h, data)

# To-dos: parallegram, rhombus, isoscles trapezoid, ...


#= Polygonal prism =#
const PolygonalPrism{K,K2} = Prism{Polygon{K,K2,Nothing}}

# Below, if we called PolygonalPrism(c, ...) in the function body, it would call the inner
# constructor Prism{Polygon{K,K2,Nothing}}(c, ...) because PolygonalPrism = Prism{Polygon{K,K2,Nothing}},
# which is not what we want.
# To call the outer constructor of Prism, we should call Prism(c, ...) instead of PolygonalPrism(c, ...).
PolygonalPrism(c::SVector{3,<:Real},
               v::SMatrix{2,K,<:Real},  # 2D coordinates of base vertices in projected prism coordinates
               h::Real=Inf,
               a::SVector{3,<:Real}=SVector(0.0,0.0,1.0),
               data=nothing) where {K} =
    (â = normalize(a); Prism(c, Polygon(v), h, [orthoaxes(â)... â], data))

PolygonalPrism(c::AbstractVector{<:Real},  # center of prism
               v::AbstractMatrix{<:Real},  # vertices of base polygon
               h::Real=Inf,  # height of prism
               a::AbstractVector{<:Real}=[0.0,0.0,1.0],  # axis direction of prism
               data=nothing) =
    (K = size(v,1); PolygonalPrism(SVector{3}(c), SMatrix{2,K}(v), h, SVector{3}(a), data))

# Return the bounds of the center cut with respect to the prism center.
function bounds_ctrcut(s::PolygonalPrism{K}) where {K}
    p = inv(s.p)  # projection matrix to prism coordinates: rows are not only unit vectors, but also orthogonal
    v = [s.b.v; @SMatrix(zeros(1,K))]  # SMatrix{3,K}: 3D vectices in prism axis coordinates
    w = p * v  # SMatrix{3,K}: vertices in external coordinates

    return minimum(w, dims=Val(2))[:,1], maximum(w, dims=Val(2))[:,1]
end
