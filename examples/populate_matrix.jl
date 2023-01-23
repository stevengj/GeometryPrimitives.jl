# In this example, we populate a simple matrix based on a geometry definition.
# Specifically, we loop over a series of grid points and decide what kind of
# material should occupy that cell by searching through the geometry
# specification. We show how one can forward the machinery underneath
# GeometryPrimitives.jl so that we can define our own "material" types.

using GeometryPrimitives
using Lazy

# create new material type
struct Material <: Shape{3,9}
    shape::Shape   # material shape
    ε::Float64     # material permittivity
end

# forward relevant methods from our new type
@forward Material.shape Base.in, GeometryPrimitives.bounds

function populate_matrix()
    # set up geometry
    geometry = [
        Material(Cylinder([0.0, 0.0, 0.0], 2.5, 100.0), 3.4)
    ]

    # set up domain
    sx, sy = 10.0, 10.0
    resolution = 20

    # build matrix
    Nx = Int(sx * resolution)
    Ny = Int(sy * resolution)
    M = zeros(Nx,Ny)

    # loop over grid points
    x = range(-sx/2,stop=sx/2,length=Nx)
    y = range(-sy/2,stop=sy/2,length=Ny)
    for ix = 1:Nx, iy = 1:Ny
        p = [x[ix], y[iy], 0.0] # current point on grid
        idx = findfirst(p,geometry) # search through geometry tree
        a = 1.0 # default material
        if !isnothing(idx)
            a = geometry[idx].ε
        end
        M[ix,iy] = a # fill in matrix
    end

    return M
end

M = populate_matrix()