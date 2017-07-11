module GeometryPrimitives

using Compat, StaticArrays

abstract type Shape{N,D} end # a solid geometric shape in N dimensions
Base.ndims(o::Shape{N}) where {N} = N

export Shape, normal, bounds

Base.in(x::AbstractVector, o::Shape{N}) where {N} = SVector{N}(x) in o
normal(x::AbstractVector, o::Shape{N}) where {N} = normal(SVector{N}(x), o)

include("sphere.jl")
include("box.jl")
include("ellipsoid.jl")
include("cylinder.jl")

include("kdtree.jl")

include("vxlcut.jl")

end # module
