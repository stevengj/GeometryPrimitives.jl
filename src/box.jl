export Box

immutable Box{N,D} <: Object{N}
    c::Point{N,Float64} # box center
    r::Vec{N,Float64}   # "radius" (semi-axis) in each direction
    p::Mat{N,N,Float64} # projection matrix to box coordinates
    data::D             # auxiliary data
end

function Box(c, r, axes=eye(length(c),length(c)), # columns are axes unit vectors
             data=nothing)
    length(c) == length(r) == size(axes,1) == size(axes,2) || throw(DimensionMismatch())
    return Box{length(c),typeof(data)}(c, r, inv(axes ./ sqrt(sumabs2(axes,2))),
                                       data)
end

function Base.in{N}(x::Point, b::Box{N})
    d = b.p * (x - b.c)
    for i = 1:N
        abs(d[i]) > b.r[i] && return false
    end
    return true
end

function normal{N}(x::Point, b::Box{N})
    d = b.p * (x - b.c)
    (m,i) = findmin(abs(abs(d) - b.r))
    return Vec{N}(b.p[i,1:N]) * sign(d[i])
end

function bounds(b::Box)
    a = abs(inv(b.p)*b.r)
    return (b.c-a,b.c+a)
end
