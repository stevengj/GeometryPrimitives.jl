module GeometryPrimitives
using Compat, StaticArrays
export Shape, surfpt_nearby, normal, bounds

abstract type Shape{N,L,D} end # a solid geometric shape in N dimensions (L = N*N is needed in some shapes, e.g., Box)

Base.ndims(o::Shape{N}) where {N} = N

# The following functions return Any due to the limitations of Julia's dispatch system.
# Therefore, always call them with type assertions.  See
# https://discourse.julialang.org/t/extending-base-in-type-stably/5341/12
# https://github.com/JuliaLang/julia/issues/23210
Base.in(x::AbstractVector, o::Shape{N}) where {N} = SVector{N}(x) in o
surfpt_nearby(x::AbstractVector, o::Shape{N}) where {N} = surfpt_nearby(SVector{N}(x), o)
normal(x::AbstractVector, o::Shape) = surfpt_nearby(x, o)[2]

include("box.jl")
include("cylinder.jl")
include("ellipsoid.jl")
include("sphere.jl")
include("kdtree.jl")
include("vxlcut.jl")

end # module
