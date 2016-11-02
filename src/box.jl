export Box

immutable Box{N,D} <: Object{N}
    c::SVector{N,Float64} # box center
    p::SMatrix{N,N,Float64} # projection matrix to box coordinates
    r::SVector{N,Float64}   # "radius" (semi-axis) in each direction
    data::D             # auxiliary data
end

function Box(c::AbstractVector, d::AbstractVector,
             axes=eye(length(c),length(c)), # columns are axes unit vectors
             data=nothing)
    length(c) == length(d) == size(axes,1) == size(axes,2) || throw(DimensionMismatch())
    return Box{length(c),typeof(data)}(c, inv(axes ./ sqrt(sumabs2(axes,2))), d*0.5, data)
end

function Base.in{N}(x::SVector, b::Box{N})
    d = b.p * (x - b.c)
    for i = 1:N
        abs(d[i]) > b.r[i] && return false
    end
    return true
end

function normal{N}(x::SVector, b::Box{N})
    d = b.p * (x - b.c)
    (m,i) = findmin(abs(abs(d) - b.r))
    return SVector{N}(b.p[i,1:N]) * sign(d[i])
end

function bounds(b::Box)
    a = abs(inv(b.p)*b.r)
    return (b.c-a,b.c+a)
end
