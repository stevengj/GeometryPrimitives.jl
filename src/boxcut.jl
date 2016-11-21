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

# a permutation p of 8 bits.  bit i -> bit p[i+1], where i=0 is the
# least significant bit.  BitPerm(0,1,2,3,4,5,6,7) is the identity.
# (performance doesn't matter since this is only used at compile time)
immutable BitPerm
    p::NTuple{8,UInt8}
    function BitPerm(p::NTuple{8,UInt8})
        # sort(collect(p)) != collect(0:7) && error("$p is not a permutation of 0:7")
        new(p)
    end
    BitPerm(b0,b1,b2,b3,b4,b5,b6,b7) = BitPerm((UInt8(b0),UInt8(b1),UInt8(b2),UInt8(b3),UInt8(b4),UInt8(b5),UInt8(b6),UInt8(b7)))
end
Base.getindex(p::BitPerm, i::Integer) = p.p[i+1]
bit(b, mask) = ifelse(b & mask != zero(b), one(b), zero(b))
Base.:*(p::BitPerm, b::UInt8) = # permute b according to p:
    UInt8(bit(b,0x01)<<p[0] + bit(b,0x02)<<p[1] + bit(b,0x04)<<p[2] + bit(b,0x08)<<p[3] +
          bit(b,0x10)<<p[4] + bit(b,0x20)<<p[5] + bit(b,0x40)<<p[6] + bit(b,0x80)<<p[7])
Base.:*(p1::BitPerm, p2::BitPerm) = BitPerm(p1[p2[0]],p1[p2[1]],p1[p2[2]],p1[p2[3]],p1[p2[4]],p1[p2[5]],p1[p2[6]],p1[p2[7]])

# bit permutations corresponding to x -> (1-x) mirror flips that flip
# the corner of bits 0..7 to (0,0,0) for corners of 3d box
const bitflip3 = [
    BitPerm(0,1,2,3,4,5,6,7), # bit 0  [0,0,0]  0x01
    BitPerm(1,0,3,2,5,4,7,6), # bit 1  [0,0,1]  0x02
    BitPerm(2,3,0,1,6,7,4,5), # bit 2  [0,1,0]  0x04
    BitPerm(3,2,1,0,7,6,5,4), # bit 3  [0,1,1]  0x08
    BitPerm(4,5,6,7,0,1,2,3), # bit 4  [1,0,0]  0x10
    BitPerm(5,4,7,6,1,0,3,2), # bit 5  [1,0,1]  0x20
    BitPerm(6,7,4,5,2,3,0,1), # bit 6  [1,1,0]  0x40
    BitPerm(7,6,5,4,3,2,1,0), # bit 7  [1,1,1]  0x80
]
# bit permutations corresponding to the 6 rotations of the
# cube that preserve (0,0,0) and (1,1,1) (bits 0 and 7)
const rotate3 = [
BitPerm(0,1,2,3,4,5,6,7), # identity
BitPerm(0,2,4,6,1,3,5,7), # rotate 120° counter-clockwise
BitPerm(0,4,1,5,2,6,3,7), # rotate 120° clockwise
BitPerm(0,2,1,3,4,6,5,7), # mirror flip 1
BitPerm(0,1,4,5,2,3,6,7), # mirror flip 2
BitPerm(0,4,2,6,1,5,3,7), # mirror flip 3
]
# Find the permutation, if any, that maps a0 to b.
# composed from the bitflip3 and rotate3 operations.
# Assumes bit 0 of a0 is set.  Returns a nullable.
function findperm3(a0::UInt8, b::UInt8)
    a0 == b && return Nullable(BitPerm(0,1,2,3,4,5,6,7))
    count_ones(a0) != count_ones(b) && return Nullable{BitPerm}()
    # just do a brute-force search, since performance is not critical:
    for i in leading_zeros(b):7
        if (0x01 << i) & b != 0x00 # bit i is set: rotate to it from [0,0,0]
            for k = 1:length(rotate3)
                p = bitflip3[i+1]*rotate3[k]
                p*a0 == b && return Nullable(p)
            end
        end
    end
    return Nullable{BitPerm}()
end

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
    b&0x02!=0x00 && return __boxcut(, Plane(SVector(p.p[1],p.p[2],-p.p[3]), p.c - p.p[3])
    b&0x04!=0x00 && return __boxcut(, Plane(SVector(p.p[1],-p.p[2],p.p[3]), p.c - p.p[2])
    b&0x10!=0x00 && return __boxcut(, Plane(SVector(-p.p[1],p.p[2],p.p[3]), p.c - p.p[1])
    b&0x08!=0x00 && return __boxcut(, Plane(SVector(p.p[1],-p.p[2],-p.p[3]), p.c - p.p[2] - p.p[3])
    b&0x20!=0x00 && return __boxcut(, Plane(SVector(-p.p[1],p.p[2],-p.p[3]), p.c - p.p[1] - p.p[3])
    b&0x40!=0x00 && return __boxcut(, Plane(SVector(-p.p[1],-p.p[2],p.p[3]), p.c - p.p[1] - p.p[2])
    #= b&0x80!=0x00 && =# return __boxcut(, Plane(SVector(-p.p[1],-p.p[2],-p.p[3]), p.c - p.p[1] - p.p[2] - p.p[3])
end

# boxcut function given bitmask b, canonicalized so that
# first bit is set (0,0,0 in p) and ≤ 4 bits are set (≤ 4 corners in p)
function __boxcut(b::UInt8, p::Plane{3})
    b==0x01 && return boxcut1(s(p,0,1),s(p,0,2),s(p,0,4))

    b==0x03 && return boxcut2(s(p,0,2),s(p,0,4),s(p,1,3),s(p,1,5))
    b==0x05 && return boxcut2(s(p,0,1),s(p,0,4),s(p,2,3),s(p,2,6))
    b==0x11 && return boxcut2(s(p,0,1),s(p,0,2),s(p,4,5),s(p,4,6))

    b==0x07 && return boxcut3(s(p,0,4),s(p,1,5),s(p,2,6),s(p,1,3),s(p,2,3))
    b==0x13 && return boxcut3(s(p,0,2),s(p,4,6),s(p,1,3),s(p,4,5),s(p,1,5))
    b==0x15 && return boxcut3(s(p,0,1),s(p,2,3),s(p,4,5),s(p,2,6),s(p,4,6))
    b==0x45 && return boxcut3(s(p,2,3),s(p,6,7),s(p,0,1),s(p,6,4),s(p,0,4))
    b==0x51 && return boxcut3(s(p,4,5),s(p,0,1),s(p,6,7),s(p,0,2),s(p,6,2))
    b==0x0b && return boxcut3(s(p,1,5),s(p,3,7),s(p,0,4),s(p,3,2),s(p,0,2))
    b==0x0d && return boxcut3(s(p,2,6),s(p,0,4),s(p,3,7),s(p,0,1),s(p,3,1))
    b==0x23 && return boxcut3(s(p,1,3),s(p,0,2),s(p,5,7),s(p,0,4),s(p,5,4))
    b==0x31 && return boxcut3(s(p,4,6),s(p,5,7),s(p,0,2),s(p,5,1),s(p,0,1))

    b==0x33 && return boxcut4f(s(p,0,2),s(p,4,6),s(p,1,3),s(p,5,7))
    b==0x55 && return boxcut4f(s(p,0,1),s(p,2,3),s(p,4,5),s(p,6,7))
    b==0xff && return boxcut4f(s(p,0,4),s(p,1,5),s(p,2,6),s(p,3,7))

    b==0x17 && return boxcut4(s(p,1,5),s(p,4,5),s(p,2,3),s(p,1,3),s(p,4,6),s(p,2,6))
    b==0x4d && return boxcut4(s(p,3,1),s(p,0,1),s(p,6,7),s(p,3,7),s(p,0,4),s(p,6,4))
    b==0x71 && return boxcut4(s(p,6,2),s(p,0,2),s(p,5,7),s(p,6,7),s(p,0,1),s(p,5,1))
    b==0x2b && return boxcut4(s(p,5,4),s(p,0,4),s(p,3,7),s(p,5,7),s(p,0,2),s(p,3,2))

    error("unhandled intersection")
end

# return the intersection point, in [0,1] if the edge intersects p,
# of the plane p with the cube edge from corner c1 to c2.
function s{N}(p::Plane{N},c1::SVector{N},c2::SVector{N})
    # Solve the equation p'*(c1 + α*(c2-c1)) = c for α,
    # hence α = (c - p'*c1) / (p'*(c2-c1))
    return (p.c - dot(p.p,c1)) / dot(p.p, c2-c1)
end
function s(p::Plane, c1::Integer, c2::Integer)
    c = corners(p)
    return s(p, c[c1-1], c[c2-1])
end

# box volume when plane contains one vertex, and intersects adjacent
# edges at e1,e2,e3
boxcut1(e1,e2,e3) = 0.16666666666666666 * e1*e2*e3

# box volume when plane contains 2 vertices
boxcut2(e01,e02,e11,e12) = 0.25 * (e01*e02 + e11*e12)

# box volume when plane contains 3 vertices
boxcut3(a,b,c,d,e,f) = 0.16666666666666666 * (a+b+c) + 0.125*(b+c)*(1-(1-d)*(1-e))

# box volume when plane contains 4 vertices, all in the same face
boxcut4f(a,b,c,d) = 0.25*(a+b+c+d)

# box volume when plane contains 4 vertices, not all on same face
boxcut4(a,b,c,d,e,f) = 0.16666666666666666 + ...
