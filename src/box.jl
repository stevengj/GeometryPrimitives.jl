export Box

mutable struct Box{N,L,D} <: Shape{N,L,D}
    c::SVector{N,Float64} # box center
    r::SVector{N,Float64}   # "radius" (semi-axis) in each direction
    p::SMatrix{N,N,Float64,L} # projection matrix to box coordinates
    data::D             # auxiliary data
    Box{N,L,D}(c,r,p,data) where {N,L,D} = new(c,r,p,data)  # suppress default outer constructor
end

Box(c::SVector{N}, d::SVector{N},
    axes::SMatrix{N,N,<:Real,L}=SMatrix{N,N}(1.0I),  # columns are axes unit vectors
    data::D=nothing) where {N,L,D} =
    Box{N,L,D}(c, 0.5d, inv(axes ./ sqrt.(sum(abs2,axes,dims=Val(1)))), data)

Box(c::AbstractVector, d::AbstractVector, axes=Matrix{Float64}(I,length(c),length(c)), data=nothing) =
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

function surfpt_nearby(x::SVector{N}, b::Box{N}) where {N}
    ax = inv(b.p)  # axes
    n = (b.p ./ sqrt.(sum(abs2,b.p,dims=Val(2))[:,1]))  # rows are direction normals; do not take transpose

    # Below, θ[i], the angle between ax[:,i] and n[i,:], is always acute, because the
    # diagonal entries of ax * b.p are positive (= 1).
    cosθ = sum(ax.*n', dims=Val(1))[1,:]
    # cosθ = diag(n*ax)  # faster than SVector(ntuple(i -> ax[:,i]⋅n[i,:], Val(N)))
    # @assert all(cosθ .≥ 0)

    d = b.p * (x - b.c)
    n = n .* copysign.(1.0,d)  # operation returns SMatrix (reason for leaving n untransposed)
    absd = abs.(d)
    onbound = abs.(b.r.-absd) .≤ Base.rtoldefault(Float64) .* b.r  # basically b.r .≈ absd but faster
    isout = (b.r.<absd) .| onbound
    ∆ = (b.r .- absd) .* cosθ
    if count(isout) == 0  # x strictly inside box
        l∆x, i = findmin(∆)
        nout = n[i,:]
        ∆x = l∆x * nout
    else  # x outside box or on boundary in one or multiple directions
        ∆x = n' * (∆ .* isout)  # project only in isout directions
        nout = all(.!isout .| onbound) ? n'*onbound : -∆x  # "if onbound in projected directions"
        nout = normalize(nout)
    end

    return x+∆x, nout
end

translate(b::Box{N,L,D}, ∆::SVector{N}) where {N,L,D} = Box{N,L,D}(b.c+∆, b.r, b.p, b.data)

signmatrix(b::Box{1}) = SMatrix{1,1}(1)
signmatrix(b::Box{2}) = SMatrix{2,2}(1,1, -1,1)
signmatrix(b::Box{3}) = SMatrix{3,4}(1,1,1, -1,1,1, 1,-1,1, 1,1,-1)

function bounds(b::Box)
    A = inv(b.p) .* b.r'
    m = maximum(abs.(A * signmatrix(b)), dims=Val(2))[:,1] # extrema of all 2^N corners of the box
    return (b.c-m,b.c+m)
end
