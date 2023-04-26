"""simple_smoothing.jl

This example shows how one can use the underlying `surfpt_nearby()` and
`volfrac()` routines to perform simple subpixel smoothing (in this case, a
linear interpolation between a cylinder and a background void).

Since we only assume a single geometric shape lies in the domain, we already
know what the foreground and background shapes are, and thus don't need to
search the tree at each pixel. If we did have a more complicated geometry
stackup, we would need to carefully determine what the foreground and background
materials are within each voxel (there could be more than 2 materials, in which
case there's no clear convention).
"""

using GeometryPrimitives
using Makie
using StaticArrays

function populate_matrix()
    # set up geometry
    geometry = [Cylinder([2.0, -3.0, 0.0], 1.0, 100.0)]

    # arbitrary material parameters
    low = 1.0
    high = 3.0

    # set up domain
    sx, sy = 10.0, 10.0
    resolution = 5 # voxels/unit

    # build matrix
    Nx = Int(sx * resolution)
    Ny = Int(sy * resolution)
    vxl_size = SVector(1.0 / resolution, 1.0 / resolution, 0.0)
    M = zeros(Nx, Ny)

    # loop over grid points and fill the matrix
    x = range(-sx / 2, stop = sx / 2, length = Nx)
    y = range(-sy / 2, stop = sy / 2, length = Ny)
    for ix = 1:Nx, iy = 1:Ny
        p = SVector(x[ix], y[iy], 0.0) # current voxel center
        vxl = (p - vxl_size / 2.0, p + vxl_size / 2.0) # current voxel on the grid

        # We only have one shape, so no need to "search" for any other shapes.
        # But normally, one would want to search for a "foreground" shape and a
        # "background" shape (within every pixel) and do averaging on those.
        nearby = surfpt_nearby(p, geometry[1])
        # compute the fill fraction
        fill_fraction = volfrac(vxl, nearby[2], nearby[1])
        # linearly interpolate based on the fill fraction
        M[ix, iy] = low + (fill_fraction) * (high - low)
    end

    return M
end

function run_and_plot()
    M = populate_matrix()
    display(heatmap(M))
end

run_and_plot()
