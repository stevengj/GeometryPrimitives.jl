export Box

type Box{N,D} <: Object{N}
    c::SVector{N,Float64} # box center
    p::SMatrix{N,N,Float64} # projection matrix to box coordinates
    r::SVector{N,Float64}   # "radius" (semi-axis) in each direction
    data::D             # auxiliary data
end

function Box(c::AbstractVector, d::AbstractVector,
             axes=eye(length(c),length(c)), # columns are axes unit vectors
             data=nothing)
    length(c) == length(d) == size(axes,1) == size(axes,2) || throw(DimensionMismatch())
    return Box{length(c),typeof(data)}(c, inv(axes ./ sqrt.(sum(abs2,axes,1))), d*0.5, data)
end

function Base.in{N}(x::SVector{N}, b::Box{N})
    d = b.p * (x - b.c)
    for i = 1:N
        abs(d[i]) > b.r[i] && return false
    end
    return true
end

function normal{N}(x::SVector{N}, b::Box{N})
    d = b.p * (x - b.c)
    (m,i) = findmin(abs.(abs.(d) - b.r))
    return SVector{N}(b.p[i,:]) * sign(d[i])
end

signmatrix(b::Object{1}) = SMatrix{1,2}(1,-1)
signmatrix(b::Object{2}) = SMatrix{2,4}(1,1, -1,1, 1,-1, -1,-1)
signmatrix(b::Object{3}) = SMatrix{3,8}(1,1,1, -1,1,1, 1,-1,1, 1,1,-1, -1,-1,1, -1,1,-1, 1,-1,-1, -1,-1,-1)

function bounds(b::Box)
    # Below, b.p' .* b.r' would have been conceptually better because its "columns"
    # are scaled axes vectors.  However, then the result is not SMatrix because b.r'
    # is not SVector.  Then, we cannot use maximum(..., Val{2}), which is type-stable
    # for SMatrix.  A workaround is to calculate b.p .* b.r and take transpose.
    A = (b.p .* b.r)'  # SMatrix

    m = maximum(A * signmatrix(b), Val{2})[:,1] # extrema of all 2^N corners of the box
    return (b.c-m,b.c+m)
end
