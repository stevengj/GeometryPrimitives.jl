export Isosceles, Trapezoid  # constructor-like methods (a.k.a factory methods)

#= Special polygons =#
# The constructors for the special polygons are defined as "factory methods" that construct
# vertices and pass them to Polygon().

# To-dos: parallegram, rhombus, ...

# Isosceles triangle
function Isosceles(base::Tuple2{SReal{2}},  # (end point 1, end point 2): two end points of base
                   h::Real)  # height drawn normal to base; direction is π/2 from base[2]-base[1]
    m = (base[1] + base[2]) / 2  # midpoint of base
    b̂ = normalize(base[2] - base[1])  # unit direction of base
    ĥ = SVec(-b̂[2], b̂[1])  # unit direction of height
    p = m + h.*ĥ  # apex

    v = [base[1] base[2] p]  # vertices

    return Polygon(v)
end

Isosceles(base::Tuple2{AbsVecReal}, h::Real) = Isosceles(SVec{2}.(base), h)

function Isosceles(c::SReal{2},  # midpoint of base
                   b::Real,  # base length
                   h::Real,  # height
                   θ::Real=0.0)  # height direction measured from +y-direction (= base direction measured from +x-direction )
    b2 = 0.5b
    b̂ = SVec(cos(θ), sin(θ))
    base = (c - b2.*b̂, c + b2.*b̂)

    return Isosceles(base, h)
end

Isosceles(c::AbsVecReal, b::Real, h::Real, θ::Real=0.0) = Isosceles(SVec{2}(c), b, h, θ)

# Trapezoid
function Trapezoid(base::Tuple2{SReal{2}},  # (end point 1, end point 2): two end points of base
                   h::Real,  # height drawn normal to base; direction is π/2 from base[2]-base[1]
                   ϕ::Tuple2{Real})  # (base angle 1, base angle 2)
    b̂ = normalize(base[2] - base[1])  # unit direction of base
    ĥ = SVec(-b̂[2], b̂[1])  # unit direction of height

    t₁ = base[1] + (h*cot(ϕ[1])) .* b̂ + h.*ĥ
    t₂ = base[2] - (h*cot(ϕ[2])) .* b̂ + h.*ĥ

    v = [base[1] base[2] t₂ t₁]

    return Polygon(v)
end

Trapezoid(base::Tuple2{AbsVecReal}, h::Real, ϕ::Tuple2{Real}) =
    Trapezoid(SVec{2}.(base), h, ϕ)

function Trapezoid(c::SReal{2},  # midpoint of base
                   b::Real,  # base length
                   h::Real,  # height
                   ϕ::Tuple2{Real},  # (base angle 1, base angle 2)
                   θ::Real=0.0)  # height direction measured from +y-direction (= base direction measured from +x-direction)
    b2 = 0.5b
    b̂ = SVec(cos(θ), sin(θ))
    base = (c - b2.*b̂, c + b2.*b̂)

    return Trapezoid(base, h, ϕ)
end

Trapezoid(c::AbsVecReal, b::Real, h::Real, ϕ::Tuple2{Real}, θ::Real=0.0) =
    Trapezoid(SVec{2}(c), b, h, ϕ, θ)


# Isosceles trapezoid
Trapezoid(base::Tuple2{AbsVecReal}, h::Real, ϕ::Real) = Trapezoid(base, h, (ϕ,ϕ))
Trapezoid(c::AbsVecReal, b::Real, h::Real, ϕ::Real, θ::Real=0.0) = Trapezoid(c, b, h, (ϕ,ϕ), θ)
