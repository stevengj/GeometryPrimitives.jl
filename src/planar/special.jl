export Isosceles, Trapezoid  # constructor-like methods (a.k.a factory methods)

#= Special polygons =#
# The constructors for the special polygons are defined as "factory methods" that construct
# vertices and pass them to Polygon().

# To-dos: parallegram, rhombus, ...

# Isosceles triangle
function Isosceles(base::NTuple{2,SVector{2,<:Real}},  # (end point 1, end point 2): two end points of base
                   h::Real)  # height drawn normal to base; direction is π/2 from base[2]-base[1]
    m = (base[1] + base[2]) / 2  # midpoint of base
    bvec = normalize(base[2] - base[1])  # unit direction of base
    hvec = @SVector [-bvec[2], bvec[1]]  # unit direction of height
    p = m + h.*hvec  # apex

    v = [base[1] base[2] p]  # vertices

    return Polygon(v)
end

Isosceles(base::NTuple{2,AbstractVector{<:Real}}, h::Real) = Isosceles(SVector{2}.(base), h)

# Trapezoid
function Trapezoid(base::NTuple{2,SVector{2,<:Real}},  # (end point 1, end point 2): two end points of base
                   h::Real,  # height drawn normal to base; direction is π/2 from base[2]-base[1]
                   θ::NTuple{2,Real})  # (base angle 1, base angle 2)
    bvec = normalize(base[2] - base[1])  # unit direction of base
    hvec = @SVector [-bvec[2], bvec[1]]  # unit direction of height

    t1 = base[1] + (h * cot(θ[1])) .* bvec + h .* hvec
    t2 = base[2] - (h * cot(θ[2])) .* bvec + h .* hvec

    v = [base[1] base[2] t2 t1]

    return Polygon(v)
end

Trapezoid(base::NTuple{2,AbstractVector{<:Real}}, h::Real, θ::NTuple{2,Real}) =
    Trapezoid(SVector{2}.(base), h, θ)

# Isosceles trapezoid
Trapezoid(base::NTuple{2,AbstractVector{<:Real}}, h::Real, θ::Real) = Trapezoid(base, h, (θ,θ))
