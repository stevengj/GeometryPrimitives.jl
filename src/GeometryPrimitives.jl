module GeometryPrimitives
using AbbreviatedTypes
using LinearAlgebra
using Statistics: mean
using Makie

export Shape, Shape1, Shape2, Shape3
export level, surfpt_nearby, normal, bounds, translate

abstract type Shape{N,N²} end # a solid geometric shape in N dimensions (N² = N*N is needed in some shapes, e.g., Cuboid)
const Shape1 = Shape{1,1}
const Shape2 = Shape{2,4}
const Shape3 = Shape{3,9}

Base.ndims(o::Shape{N}) where {N} = N

# The following functions return Any due to the limitations of Julia's dispatch system.
# Therefore, always call them with return type assertions.  See
# https://discourse.julialang.org/t/extending-base-in-type-stably/5341/12
# https://github.com/JuliaLang/julia/issues/23210
level(x::AbsVecReal, o::Shape{N}) where {N} = level(SVec{N}(x), o)
Base.in(x::AbsVecReal, o::Shape{N}) where {N} = level(x,o) ≥ 0
surfpt_nearby(x::AbsVecReal, o::Shape{N}) where {N} = surfpt_nearby(SVec{N}(x), o)
normal(x::AbsVecReal, o::Shape) = surfpt_nearby(x, o)[2]  # outward direction even for x inside o
translate(s::Shape{N}, ∆::AbsVecReal) where {N} = translate(s, SVec{N}(∆))
translate(s::Shape{N}, ∆::SReal{N}) where {N} = (s2 = deepcopy(s); s2.c += ∆; s2)  # default implementation

function orthoaxes(n::SReal{3})
    u_temp = abs(n[3]) < abs(n[1]) ? SVec(0,0,1) : SVec(1,0,0)
    v = normalize(n × u_temp)
    u = v × n

    return u, v
end

function orthoaxes(n::SReal{N}) where {N}
    u_temp = abs(n[3]) < abs(n[1]) ? SVec(0,0,1) : SVec(1,0,0)
    v = normalize(n × u_temp)
    u = v × n

    return u, v
end


include("hyper/hyper.jl")
include("planar/planar.jl")
include("prism/prism.jl")
include("util/util.jl")

end # module
