# Assume that the prism axis is orthogonal to the prism base.  This makes it much easier to
# implement surfpt_nearby by separating the axis and base dimensions (i.e., the closest
# point on the surface can be found independently in axis and base dimensions).
#
# Supporting skewed prisms is not impossible, but this will be done in the future if the
# demand is high.

# It is possible to combine Cylinder and Prism by defining a more general prism type that
# accepts any base shape.  We could use Shape{2} to create a 2D base shape and store it as
# a base shape.  However there is one problem.  We want to use
#
#   surfpt_nearby(x_base::SVec{2}, shape_base::Shape{2})
#
# to find the closest point in the base dimension (such a point is on the side of the prism),
# and compare the distance to it with the distance from x to the base plane along the axis
# dimension.  Then, if x is inside the prism for example, then this comparison will tell us
# which of the base and side is the closest surface from x.  However, this means that we
# have to take the difference between the output of the above function and x_base to figure
# out the distance in the base dimension.  Provided that for many Shape types ∆x is first
# calculated inside surfpt_nearby and x+∆x is returned as the found surface point, having to
# recover ∆x by subtraction is not only inefficient, but also suffers from loss of
# significant digits.
#
# Maybe we should have not implemented surfpt_nearby, but a different function that returns
# ∆x, because the surface point x + ∆x can be easily calculated.  Let's consider fixing this
# in the future.
#
# Hmmm...  Because the surface point is not calculated by adding ∆x to x for Ball and
# Ellipse, maybe implemeting a function returning ∆x is not a very good idea.  Then, it may
# not be a bad idea to calculate ∆x by subtracting x from the surface point.  Let's try to
# implement the general prism.

export Prism

mutable struct Prism{B<:Shape{2}} <: Shape{3,9}
    c::SFloat{3}  # prism center
    b::B  # base shape described in prism coordinates (i.e, when translating prism, do not need to translate b)
    h2::Float  # height * 0.5
    p::S²Float{3,9}  # projection matrix to prism coordinates; must be orthonormal (see surfpt_nearby)
    Prism{B}(c,b,h2,p) where {B} = new(c,b,h2,p)  # suppress default outer constructor
end

Prism(c::SReal{3},
      b::B,
      h::Real=Inf,
      axes::S²Real{3,9}=S²Float{3,9}(I)  # columns are axes vectors: first two columns span prism base, and last column is prism axis
      ) where {B<:Shape{2}} =
    Prism{B}(c, b, 0.5h, inv(axes ./ sqrt.(sum(abs2,axes,dims=Val(1)))))

Prism(c::AbsVecReal, b::Shape{2}, h::Real=Inf, axes::AbsMatReal=MatFloat(I,length(c),length(c))) =
    Prism(SVec{3}(c), b, h, S²Mat{3}(axes))

Base.:(==)(s1::Prism, s2::Prism) = s1.c==s2.c && s1.b==s2.b && s1.h2==s2.h2 && s1.p==s2.p
Base.isapprox(s1::Prism, s2::Prism) = s1.c≈s2.c && s1.b≈s2.b && s1.h2≈s2.h2 && s1.p≈s2.p
Base.hash(s::Prism, h::UInt) = hash(s.c, hash(s.b, hash(s.h2, hash(s.p, hash(:Prism, h)))))

function level(x::SReal{3}, s::Prism)
    y = s.p * (x - s.c)  # coordinates after projection
    ya = y[3]  # scalar: coordinate in axis dimension
    yb = y[SVec(1,2)]  # SVec{2}: coordinate in base dimensions

    return min(1.0 - abs(ya)/s.h2, level(yb,s.b))
end

function surfpt_nearby(x::SReal{3}, s::Prism)
    ax = inv(s.p)  # prism axes: columns are not only unit vectors, but also orthogonal

    y = s.p * (x - s.c)  # x in prism coordinates
    ya = y[3]  # scalar: coordinate in axis dimension
    yb = y[SVec(1,2)]  # SVec{2}: coordinates in base dimensions

    la = abs(ya)
    abs∆a = abs(s.h2 - la)  # scalar: distance between x and base point closest to x
    surfa = SVec(yb.data..., copysign(s.h2, ya))  # SVec{3}: coordinates of base point closest to x
    nouta = SVec(0.0, 0.0, copysign(1.0, ya))  # SVec{3}: outward direction normal at surfa
    onbnda = abs∆a ≤ Base.rtoldefault(Float) * s.h2
    isouta = s.h2<la || onbnda

    surfb2, noutb2 = surfpt_nearby(yb, s.b)  # (SVec{2}, SVec{2}): side point closest to x and outward direction normal to side there
    abs∆b = norm(surfb2 - yb)  # scalar: distance between x and side point closest to x
    surfb = SVec(surfb2.data..., ya)  # SVec{3}: coordinates of side point closest to x
    noutb = SVec(noutb2.data..., 0.0)  # SVec{3}: outward direction normal to side surface at surfb
    basesize = abs.((-)(bounds(s.b)...))  # SVec{2}: size of bounding rectancle of base
    onbndb = abs∆b ≤ Base.rtoldefault(Float) * max(basesize.data...)
    isoutb = yb∉s.b || onbndb

    if isouta && isoutb  # x outside in both axis and base dimensions
        surf = SVec(surfb[1], surfb[2], surfa[3])
        nout = (onbnda && onbndb) ? (noutb + nouta) : (y - surf)
        nout = norm(nout)==Inf ? isinf.(nout) .* sign.(nout) : normalize(nout)  # e.g., return [0,0,-1] for nout = [1,-2,-Inf]
    elseif !isouta && isoutb  # x outside in base dimensions, but inside prism in axis dimension
        (surf, nout) = (surfb, noutb)
    elseif isouta && !isoutb # x outside in axis dimension, but inside prism in base dimensions
        (surf, nout) = (surfa, nouta)
    else  # !isouta && !isoutb: x strictly inside prism
        (surf, nout) = (abs∆a ≤ abs∆b) ? (surfa, nouta) : (surfb, noutb)
    end

    return ax*(surf+s.c), ax*nout
end

function bounds(s::Prism)
    ax = inv(s.p)  # prism axes: columns are not only unit vectors, but also orthogonal
    a = ax[:,3]  # SVec{3}
    h2a = s.h2 * a

    l0, u0 = bounds_ctrcut(s)  # (SVec{3}, SVec{3})
    l1, u1 = l0+h2a, u0+h2a  # (SVec{3}, SVec{3})
    l2, u2 = l0-h2a, u0-h2a  # (SVec{3}, SVec{3})

    return min.(l1,l2)+s.c, max.(u1,u2)+s.c
end


include("cylinder.jl")
include("polygonal.jl")
include("sectoral.jl")
