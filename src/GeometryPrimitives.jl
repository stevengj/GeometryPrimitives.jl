module GeometryPrimitives

using Compat, StaticArrays

@compat abstract type Shape{N} end # a solid geometric shape in N dimensions
Base.ndims{N}(o::Shape{N}) = N

export Shape, normal, bounds

Base.in{N}(x::AbstractVector, o::Shape{N}) = SVector{N}(x) in o
normal{N}(x::AbstractVector, o::Shape{N}) = normal(SVector{N}(x), o)

include("sphere.jl")
include("box.jl")
include("ellipsoid.jl")
include("cylinder.jl")

include("kdtree.jl")

include("vxlcut.jl")

end # module
