# Put all the graphical tests here.

@testset "periodize" begin
    c = Cylinder([0,0,0], 1/4, 1)  # circular base
    # c = PolygonalPrism([0,0,0], [1/4 1/4; -1/4 1/4; -1/4 -1/4; 1/4 -1/4], 1)  # square base
    # c = Prism([0,0,0], regpoly(5, 1/4), 1)  # regular pentagonal base
    ∆range = Cuboid([0,0,0], [8,2*√3,1])
    A = [1 0 0; 0.5 √3/2 0; 0 0 1]'
    c_array = periodize(c, A, ∆range)

    b = Cuboid([0,0,0], [12,4*√3,1])
    bnd = bounds(b)
    x = range(bnd[1][1], stop=bnd[2][1], length=120*5)
    y = range(bnd[1][2], stop=bnd[2][2], length=70*5)
    Nx, Ny = length(x), length(y)
    X = repeat(x, outer=(1,Ny))
    Y = repeat(y', outer=(Nx,1))
    V = zeros(Bool, Nx, Ny)

    kd = KDTree(c_array)
    for j = 1:Ny, i = 1:Nx
        p = [x[i], y[j], 0.0]
        s = findfirst(p, kd)
        V[i,j] = s≠nothing
    end

    using PyPlot
    clf()
    axis("equal")
    pcolor(X, Y, V, cmap="gray_r")

    using PyCall
    @pyimport matplotlib.patches as patch
    b∆range = bounds(∆range)
    cnr = b∆range[1][1:2]
    w, h = -((-)(b∆range...)[1:2])
    gca()[:add_artist](patch.Rectangle(cnr, w, h, fill=false, edgecolor="r", linestyle="--"))
end  # @testset "periodize"
