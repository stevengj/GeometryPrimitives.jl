@testset "triangular cylinder 3D" begin
    vxl = (SVec(0,0,0), SVec(1,1,1))
    nout = SVec(1,1,0)

    @test_nowarn (r₀ = SVec(0.5,0,0); @inferred(volfrac(vxl, nout, r₀)))
    @test (r₀ = SVec(0.5,0,0); volfrac(vxl, nout, r₀) ≈ 0.125)
    @test (r₀ = SVec(0.5,0,0); volfrac(vxl, -nout, r₀) ≈ 0.875)
    @test (r₀ = SVec(1,0,0); volfrac(vxl, nout, r₀) ≈ 0.5)
    @test (r₀ = SVec(1,0,0); volfrac(vxl, -nout, r₀) ≈ 0.5)
    @test (r₀ = SVec(1,0,0); nout = SVec(1,2,0); volfrac(vxl, nout, r₀) ≈ 0.25)
    @test (r₀ = SVec(1,0,0); nout = SVec(1,2,0); volfrac(vxl, -nout, r₀) ≈ 0.75)
end  # @testset "triangular cylinder 3D"

@testset "triangular cylinder 2D" begin
    vxl = (SVec(0,0), SVec(1,1))
    nout = SVec(1,1)

    @test_nowarn (r₀ = SVec(0.5,0); @inferred(volfrac(vxl, nout, r₀)))
    @test (r₀ = SVec(0.5,0); volfrac(vxl, nout, r₀) ≈ 0.125)
    @test (r₀ = SVec(0.5,0); volfrac(vxl, -nout, r₀) ≈ 0.875)
    @test (r₀ = SVec(1,0); volfrac(vxl, nout, r₀) ≈ 0.5)
    @test (r₀ = SVec(1,0); volfrac(vxl, -nout, r₀) ≈ 0.5)
    @test (r₀ = SVec(1,0); nout = SVec(1,2); volfrac(vxl, nout, r₀) ≈ 0.25)
    @test (r₀ = SVec(1,0); nout = SVec(1,2); volfrac(vxl, -nout, r₀) ≈ 0.75)
end  # @testset "triangular cylinder 2D"

@testset "quadrangular cylinder 3D" begin
    @test_nowarn @inferred(volfrac((SVec(0,0,0),SVec(1,1,1)), SVec(1,2,0), SVec(0.5,0.5,0.5)))
    @test begin
        result = true
        for i = 1:100
            vxl = (-@SVector(rand(3)), @SVector(rand(3)))
            r₀ = mean(vxl)
            nout = randn(3)
            nout[rand(1:3)] = 0
            result &= volfrac(vxl, SVec{3}(nout), r₀)≈0.5
        end
        result
    end
end  # @testset "quadrangular cylinder 3D"

@testset "quadrangular cylinder 2D" begin
    @test_nowarn @inferred(volfrac((SVec(0,0),SVec(1,1)), SVec(1,2), SVec(0.5,0.5)))
    @test begin
        result = true
        for i = 1:100
            vxl = (-@SVector(rand(2)), @SVector(rand(2)))
            r₀ = mean(vxl)
            nout = @SVector randn(2)
            result &= volfrac(vxl, nout, r₀)≈0.5
        end
        result
    end
end  # @testset "quadrangular cylinder 2D"

@testset "general cases" begin
    # Test random cases.
    @test_nowarn @inferred(volfrac((-@SVector(rand(3)),@SVector(rand(3))), SVec(0,0,0), @SVector(randn(3))))
    @test begin
        result = true
        for i = 1:100
            vxl = (-@SVector(rand(3)), @SVector(rand(3)))
            r₀ = mean(vxl)
            nout = @SVector randn(3)
            result &= volfrac(vxl, nout, r₀)≈0.5
        end
        result
    end

    @test begin
        result = true
        for i = 1:100
            vxl = (-@SVector(rand(3)), @SVector(rand(3)))
            r₀ = @SVector randn(3)
            nout = @SVector randn(3)
            result &= (volfrac(vxl, nout, r₀) + volfrac(vxl, -nout, r₀) ≈ 1)
        end
        result
    end

    # Test boundary cases.
    vxl = (SVec(0,0,0), SVec(1,1,1))
    r₀ = SVec(0,0,0)
    @test (nout = SVec(1,0,0); volfrac(vxl, nout, r₀) ≈ 0)  # completely outside
    @test (nout = SVec(-1,0,0); volfrac(vxl, nout, r₀) ≈ 1)  # completely inside
    @test (nout = SVec(1,1,0); volfrac(vxl, nout, r₀) ≈ 0)  # completely outside
    @test (nout = SVec(-1,-1,0); volfrac(vxl, nout, r₀) ≈ 1)  # completely inside
    @test (nout = SVec(1,1,1); volfrac(vxl, nout, r₀) ≈ 0)  # completely outside
    @test (nout = SVec(-1,-1,-1); volfrac(vxl, nout, r₀) ≈ 1)  # completely inside
    @test (nout = SVec(-1,1,0); volfrac(vxl, nout, r₀) ≈ 0.5)  # rvol_tricyl()
    @test (nout = SVec(-1,-1,1); volfrac(vxl, nout, r₀) ≈ 5/6)  # rvol_gensect()
    @test (nout = SVec(1,-2,1); volfrac(vxl, nout, r₀) ≈ 0.5)  # rvol_quadsect()
    r₀ = SVec(0.5,0.5,0)
    @test (nout = SVec(-2,1,0); volfrac(vxl, nout, r₀) ≈ 0.5)  # rvol_quadcyl()

    # Test tolerance to floating-point arithmetic.
    vxl = (SVec(0,9.5,0), SVec(1,10.5,1))
    r₀ = SVec(0.5,10.0,0.5)
    nout = SVec(0.0,1/√2,1/√2)
    @test volfrac(vxl, nout, r₀) ≈ 0.5

    # Test rvol_quadsect() for nontrivial cases.
    vxl = (SVec(0,0,0), SVec(1,1,2))
    r₀ = SVec(0.5, 0.5, 0.5)
    @test begin
        result = true
        for i = 1:100
            nout = SVec(randn()/20, randn()/20, 1)
            result &= volfrac(vxl, nout, r₀)≈0.5/2
        end
        result
    end
end  # @testset "general cases"
