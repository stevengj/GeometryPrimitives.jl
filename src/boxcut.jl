# Given a plane and an axis-aligned unit cube with a corner at the origin,
# compute the volume of the cube on one side of the plane.
#
# A planar halfspace is defined by p'*x ≤ c for a vector p and a real number c.

immutable Plane{N}
    p::SVector{N,Float64}
    c::Float64
end

Base.in{N}(x::SVector{N}, p::Plane{N}) = dot(x,p.p) ≤ p.c

# corners of the unit cube in different dimensions
const corners3 = [let d = digits(n, 2, 3); SVector(d[3],d[2],d[1]); end for n in 0:7]
const corners2 = [let d = digits(n, 2, 2); SVector(d[2],d[1]); end for n in 0:3]
const corners1 = [SVector(0), SVector(1)]
corners(p::Plane{3}) = corners3
corners(p::Plane{2}) = corners2
corners(p::Plane{1}) = corners1

function boxcut{N}(p::Plane{N})
    # construct a bitmask b in which the (i-1)-th bit is set if
    # the i-th corner of the box is in p.
    c = corners(p)
    b = 0x00
    for i = 1:length(c)
        if c[i] in p
            b += 0x01 << (i-1)
        end
    end
    return _boxcut(b, p)
end

# permute the bits of b, assuming the first bit is 0, to
# the new positions given by b1..b7.  e.g. pbits(b,1,2,3,4,5,6,7) is the identity
pbits(b::UInt8, b1,b2,b3,b4,b5,b6,b7) =
    UInt8(bit(b,0x02)<<b1 + bit(b,0x04)<<b2 + bit(b,0x08)<<b3 +
          bit(b,0x10)<<b4 + bit(b,0x20)<<b5 + bit(b,0x40)<<b6 +
          bit(b,0x80)<<b7)
bit(b, mask) = ifelse(b & mask != zero(b), one(b), zero(b))

# boxcut function but given bitmask b of corners in p
function _boxcut(b::UInt8, p::Plane{3})
    b == 0x00 && return 0.0 # no corners in p
    b == 0xff && return 1.0 # all corners in p
    # if more than 4 corners are in p, then compute volume outside of p
    count_ones(b) > 4 && return 1.0 - _boxcut(~b, Plane(-p.p, -p.c))

    # check if (0,0,0) corner is in p, so we can call __boxcut:
    b&0x01!=0x00 && return __boxcut(b, p)

    # Otherwise, send x -> (1-x) (for correct coords x) to put (0,0,0) inside p.
    # For reference, the table of bit position, vectors, and bitmasks is:
    #     0  [0,0,0]  0x01
    #     1  [0,0,1]  0x02
    #     2  [0,1,0]  0x04
    #     3  [0,1,1]  0x08
    #     4  [1,0,0]  0x10
    #     5  [1,0,1]  0x20
    #     6  [1,1,0]  0x40
    #     7  [1,1,1]  0x80
    b&0x02!=0x00 && return __boxcut(pbits(0,3,2,5,4,7,6), Plane(SVector(p.p[1],p.p[2],-p.p[3]), p.c - p.p[3])
    b&0x04!=0x00 && return __boxcut(pbits(3,0,1,6,7,4,5), Plane(SVector(p.p[1],-p.p[2],p.p[3]), p.c - p.p[2])
    b&0x10!=0x00 && return __boxcut(pbits(5,6,7,0,1,2,3), Plane(SVector(-p.p[1],p.p[2],p.p[3]), p.c - p.p[1])
    b&0x08!=0x00 && return __boxcut(pbits(2,1,0,7,6,5,4), Plane(SVector(p.p[1],-p.p[2],-p.p[3]), p.c - p.p[2] - p.p[3])
    b&0x20!=0x00 && return __boxcut(pbits(4,7,6,1,0,3,2), Plane(SVector(-p.p[1],p.p[2],-p.p[3]), p.c - p.p[1] - p.p[3])
    b&0x40!=0x00 && return __boxcut(pbits(7,4,5,2,3,0,1), Plane(SVector(-p.p[1],-p.p[2],p.p[3]), p.c - p.p[1] - p.p[2])
    #= b&0x80!=0x00 && =# return __boxcut(pbits(6,5,4,3,2,1,0), Plane(SVector(-p.p[1],-p.p[2],-p.p[3]), p.c - p.p[1] - p.p[2] - p.p[3])
end

# boxcut function given bitmask b, canonicalized so that
# first bit is set (0,0,0 in p) and ≤ 4 bits are set (≤ 4 corners in p)
function __boxcut(b::UInt8, p::Plane{3})
    b==0x01 && return boxcut1(s(p,0,1),s(p,0,2),s(p,0,4))
    b==0x03 && return boxcut2(s(p,0,2),s(p,0,4),s(p,1,3),s(p,1,5))
    b==0x05 && return boxcut2(s(p,0,1),s(p,0,4),s(p,2,3),s(p,2,6))
    b==0x11 && return boxcut2(s(p,0,1),s(p,0,2),s(p,4,5),s(p,4,6))
end
