export Box

mutable struct Box{N,D,L} <: Shape{N,D}
    c::SVector{N,Float64} # box center
    r::SVector{N,Float64}   # "radius" (semi-axis) in each direction
    p::SMatrix{N,N,Float64,L} # projection matrix to box coordinates
    data::D             # auxiliary data
    Box{N,D,L}(c,r,p,data) where {N,D,L} = new(c,r,p,data)  # suppress default outer constructor
end

Box(c::SVector{N}, d::SVector{N},
    axes::SMatrix{N,N,<:Real,L}=@SMatrix(eye(N)),  # columns are axes unit vectors
    data::D=nothing) where {N,D,L} =
    Box{N,D,L}(c, 0.5d, inv((axes' ./ sqrt.(sum(abs2,axes,Val{1}))[1,:])'), data)
# Use this after StaticArrays issue 242 is fixed:
#    Box{N,D,L}(c, 0.5d, inv(axes ./ sqrt.(sum(abs2,axes,Val{1}))), data)

Box(c::AbstractVector, d::AbstractVector, axes=eye(length(c)), data=nothing) =
    (N = length(c); Box(SVector{N}(c), SVector{N}(d), SMatrix{N,N}(axes), data))

Base.:(==)(b1::Box, b2::Box) = b1.c==b2.c && b1.r==b2.r && b1.p==b2.p && b1.data==b2.data
Base.hash(b::Box, h::UInt) = hash(b.c, hash(b.r, hash(b.p, hash(b.data, hash(:Box, h)))))

function Base.in(x::SVector{N}, b::Box{N}) where {N}
    d = b.p * (x - b.c)
    for i = 1:N
        abs(d[i]) > b.r[i] && return false
    end
    return true
end

function normal(x::SVector{N}, b::Box{N}) where {N}
    d = b.p * (x - b.c)
    (m,i) = findmin(abs.(abs.(d) - b.r))
    return SVector{N}(b.p[i,:]) * sign(d[i])
end

signmatrix(b::Shape{1}) = SMatrix{1,2}(1,-1)
signmatrix(b::Shape{2}) = SMatrix{2,4}(1,1, -1,1, 1,-1, -1,-1)
signmatrix(b::Shape{3}) = SMatrix{3,8}(1,1,1, -1,1,1, 1,-1,1, 1,1,-1, -1,-1,1, -1,1,-1, 1,-1,-1, -1,-1,-1)

function bounds(b::Box)
    # Below, b.p' .* b.r' would have been conceptually better because its "columns"
    # are scaled axes vectors.  However, then the result is not SMatrix because b.r'
    # is not SVector.  Then, we cannot use maximum(..., Val{2}), which is type-stable
    # for SMatrix.  A workaround is to calculate b.p .* b.r and take transpose.
    A = (b.p .* b.r)'  # SMatrix

    m = maximum(A * signmatrix(b), Val{2})[:,1] # extrema of all 2^N corners of the box
    return (b.c-m,b.c+m)
end
