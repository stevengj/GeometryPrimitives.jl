@testset "periodize" begin
    # Square lattice
    c = Cylinder([0,0,0], 1, 5, [0,0,1])
    ∆range = Box([0,0,0], [10,10,5])
    A = [1 0 0; 0 1 0; 0 0 5]'
    c_array = periodize(c, A, ∆range)
    # @test length(c_array) == 11^2  # test if correct number of cylinders are generated

    kd = KDTree(c_array)
    result = true
    for py = -5.0:5.0, px = -5.0:5.0
        p = [px, py, 0.0]
        s = findfirst(p, kd)
        result &= s≠nothing  # test if all lattice points within ∆range are included in some cylinder
    end
    @test result
end
