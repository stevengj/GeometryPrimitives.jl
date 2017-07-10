export volfrac

const X, Y, Z = 1, 2, 3
const XYZ, YZX, ZXY = (X,Y,Z), (Y,Z,X), (Z,X,Y)
const UVW = YZX, ZXY, XYZ
const N, P = 1, 2  # negative, positive
const NP = (N, P)

corner(vxl::NTuple{2,SVector{3,<:Real}}, sx::Integer, sy::Integer, sz::Integer) =
    @SVector [vxl[sx][X], vxl[sy][Y], vxl[sz][Z]]

function corner_bits(vxl::NTuple{2,SVector{3}},  # two ends of solid diagonal of voxel
                     nout::SVector{3}, # unit outward normal of plane
                     nr₀)  # equation of plane: nout⋅(r - r₀) = 0, or nout⋅r = nr₀
    # Calculate the bit array that indicates corner contained-ness.

    # The bit array is 8 bits (= 1 byte).  The kth bit is 1 if the kth corner is
    # contained in the half space defined by the plane, an 0 otherwise.  The 8 corners
    # of the voxel are indexed 1 through 8 in the order of (---), (+--), (-+-), (++-),
    # (--+), (+-+), (-++), (+++), where, e.g., (+--) indicates the corner at the
    # +x, -y, -z corner.
    cbits = 0x00
    bit = 0x01
    for sz = NP, sy = NP, sx = NP
        r = corner(vxl, sx, sy, sz)
        if nout⋅r ≤ nr₀  # corner is inside (nout is outward normal)
            cbits |= bit
        end
        bit <<= 1
    end

    return cbits
end

function isquadsect(cbits::UInt8)
    # Determine if the corner contained-ness bit array corresponds to the case where
    # the plane crosses one of the three sets of four parallel edges of the voxel.

    return cbits==0x0F || cbits==~0x0F || cbits==0x33 || cbits==~0x33 || cbits==0x55 || cbits==~0x55

    # Equivalent to
    # n = count_ones(cbits)
    # m = count_ones(cbits & 0x0F)
    # return n==4 && iseven(m)
end

function edgedir_quadsect(cbits::UInt8)
    # For the cases where the plane crosses a set of four parallel edges of the voxel,
    # determine which direction those edges lie.

    if cbits==0x0F || cbits==~0x0F
        dir = Z
    elseif cbits==0x33 || cbits==~0x33
        dir = Y
    else
        assert(cbits==0x55 || cbits==~0x55)
        dir = X
    end

    return dir
end

function rvol_quadsect(vxl::NTuple{2,SVector{3}}, nout::SVector{3}, nr₀, cbits::UInt8)
    # Return the volume fraction for the case where the plane croses a set of four parallel
    # edges.

    const w = edgedir_quadsect(cbits)
    const ∆w = vxl[P][w] - vxl[N][w]

    const u, v, ~ = UVW[w]
    const nu, nv, nw = nout[u], nout[v], nout[w]
    const mean_cepts = 4nr₀
    for sv = NP, su = NP
        mean_cepts -= nu*vxl[su][u] + nv*vxl[sv][v]
    end
    mean_cepts /=  nw * 4∆w

    const sw = nw>0 ? N : P  # nw cannot be 0
    return abs(mean_cepts - vxl[sw][w]/∆w)
end

function rvol_gensect(vxl::NTuple{2,SVector{3}}, nout::SVector{3}, nr₀, cbits::UInt8)
    # Calculate the volume fraction of most general cases, by cutting out corners from a
    # triangular pyramid.
    # Assume count_ones(cbits) ≤ 4.  Othewise, call this function with flipped nout, nr₀,
    # cbits.

    const s = (nout.<0) .+ 1
    const c = corner(vxl, s...)  # corner coordinates
    const ∆ = vxl[P] - vxl[N]  # vxl edges
    const nc = nout .* c
    rmax, rmid, rmin =  abs.((((nr₀-sum(nc)) .+ nc) ./ nout - c) ./ ∆) # (lengths from corner to intercetps) / (voxel edges)

    # const nx, ny, nz = nout
    # const sx, sy, sz = ((nx≥0 ? N : P), (ny≥0 ? N : P), (nz≥0 ? N : P))  # signs of corner
    # const cx, cy, cz = vxl[nX][sx], vxl[nY][sy], vxl[nZ][sz]  # corner coordinates
    # const ∆x, ∆y, ∆z = vxl[nX][nP]-vxl[nX][nN], vxl[nY][nP]-vxl[nY][nN], vxl[nZ][nP]-vxl[nZ][nN]  # voxel edges
    # const nxcx, nycy, nzcz = nx*cx, ny*cy, nz*cz
    # rmax, rmid, rmin =  # (lengths from corner to intercetps) / (voxel edges)
    #     abs((nr₀-nycy-nzcz)/nx-cx)/∆x, abs((nr₀-nzcz-nxcx)/ny-cy)/∆y, abs((nr₀-nxcx-nycy)/nz-cz)/∆z


    # Sort rmax, rmin, rmin properly.
    if rmax < rmid; rmax, rmid = rmid, rmax; end
    if rmid < rmin; rmid, rmin = rmin, rmid; end
    if rmax < rmid; rmax, rmid = rmid, rmax; end

    # Calculate the volume of the triangular pyramid, and cut off appropriate corners.
    const tmax = 1 - 1/rmax
    rvol_core = 1 + tmax + tmax^2
    if rmid > 1
        tmid = 1 - 1/rmid
        rvol_core -= rmax * tmid^3
    end
    if rmin > 1
        tmin = 1 - 1/rmin
        rvol_core -= rmax * tmin^3
    end

    return rvol_core * rmin * rmid / 6
end

"""
    volfrac(voxel, nout, r₀)

Returns the volume fraction `rvol = vol(voxel ⋂ half-space) / vol(voxel)` that indicates how
much portion of the volume of a voxel is included the half-space defined by a plane.  The
result is between 0.0 and 1.0, inclusive.  Inside and outside of the half space is
determined by `nout`, which is the outward normal.  (If a point is on the `nout` side of the
plane, then it is outside.)

The voxel needs to be aligned with the Cartesian axes, and is described by `v =
([xn,yn,zn], [xp,yp,zp])` that specifies the x-, y-, z-boundaries of the voxel.  The
half-space is described by the boundary plane.  The boundary plane is described by the
outward normal vector `nout = [nx, ny, nz]` and a point `r₀ = [rx, ry, rz]` on the plane.
`nout` does not have to be normalized.
"""
function volfrac(vxl::NTuple{2,SVector{3}}, nout::SVector{3}, r₀::SVector{3})
    # Return the volume fraction rvol = vol(voxel ⋂ half-space) / vol(voxel).

    const nr₀ = nout⋅r₀
    const cbits = corner_bits(vxl, nout, nr₀)
    const n_corners = count_ones(cbits)  # number of corners contained

    if n_corners == 8  # voxel is inside half-space
        rvol = 1.
    elseif n_corners == 0  # voxel is outside half-space
        rvol = 0.
    elseif isquadsect(cbits) # plane crosses a set of four parallel edges of voxel
        rvol = rvol_quadsect(vxl, nout, nr₀, cbits)
    elseif n_corners ≤ 4 # general cases with n_corners ≤ 4
        rvol = rvol_gensect(vxl, nout, nr₀, cbits)
    else  # general cases with n_corners ≥ 5
        assert(n_corners ≥ 5)
        rvol = 1. - rvol_gensect(vxl, .-nout, -nr₀, ~cbits)
    end

    return rvol
end

volfrac(vxl::NTuple{2,SVector{2}}, nout::SVector{2}, r₀::SVector{2}) =
    volfrac((@SVector([vxl[N][1],vxl[N][2],0]), @SVector([vxl[P][1],vxl[P][2],1])),
            @SVector([nout[1], nout[2], 0]),
            @SVector([r₀[1], r₀[2], 0]))
