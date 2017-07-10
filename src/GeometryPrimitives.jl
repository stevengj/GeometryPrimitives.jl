module GeometryPrimitives

using Compat, StaticArrays

<<<<<<< HEAD
@compat abstract type Shape{N,D} end # a solid geometric shape in N dimensions
Base.ndims{N}(o::Shape{N}) = N
=======
abstract type Shape{N} end # a solid geometric shape in N dimensions
Base.ndims(o::Shape{N}) where {N} = N
>>>>>>> 0f8dd2b... Drop support for Julia 0.5, and use Julia 0.6 syntax (#7)

export Shape, normal, bounds

Base.in(x::AbstractVector, o::Shape{N}) where {N} = SVector{N}(x) in o
normal(x::AbstractVector, o::Shape{N}) where {N} = normal(SVector{N}(x), o)

include("sphere.jl")
include("box.jl")
include("ellipsoid.jl")
include("cylinder.jl")

include("kdtree.jl")

end # module
