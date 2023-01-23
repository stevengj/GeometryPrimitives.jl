@testset "triangular cylinder 3D" begin
    vxl = (SVector(0,0,0), SVector(1,1,1))
    nout = SVector(1,1,0)

    @test_nowarn (r₀ = SVector(0.5,0,0); @inferred(volfrac(vxl, nout, r₀)))
    @test (r₀ = SVector(0.5,0,0); volfrac(vxl, nout, r₀) ≈ 0.125)
    @test (r₀ = SVector(0.5,0,0); volfrac(vxl, -nout, r₀) ≈ 0.875)
    @test (r₀ = SVector(1,0,0); volfrac(vxl, nout, r₀) ≈ 0.5)
    @test (r₀ = SVector(1,0,0); volfrac(vxl, -nout, r₀) ≈ 0.5)
    @test (r₀ = SVector(1,0,0); nout = SVector(1,2,0); volfrac(vxl, nout, r₀) ≈ 0.25)
    @test (r₀ = SVector(1,0,0); nout = SVector(1,2,0); volfrac(vxl, -nout, r₀) ≈ 0.75)
end  # @testset "triangular cylinder 3D"

@testset "triangular cylinder 2D" begin
    vxl = (SVector(0,0), SVector(1,1))
    nout = SVector(1,1)

    @test_nowarn (r₀ = SVector(0.5,0); @inferred(volfrac(vxl, nout, r₀)))
    @test (r₀ = SVector(0.5,0); volfrac(vxl, nout, r₀) ≈ 0.125)
    @test (r₀ = SVector(0.5,0); volfrac(vxl, -nout, r₀) ≈ 0.875)
    @test (r₀ = SVector(1,0); volfrac(vxl, nout, r₀) ≈ 0.5)
    @test (r₀ = SVector(1,0); volfrac(vxl, -nout, r₀) ≈ 0.5)
    @test (r₀ = SVector(1,0); nout = SVector(1,2); volfrac(vxl, nout, r₀) ≈ 0.25)
    @test (r₀ = SVector(1,0); nout = SVector(1,2); volfrac(vxl, -nout, r₀) ≈ 0.75)
end  # @testset "triangular cylinder 2D"

@testset "quadrangular cylinder 3D" begin
    @test_nowarn @inferred(volfrac((SVector(0,0,0),SVector(1,1,1)), SVector(1,2,0), SVector(0.5,0.5,0.5)))
    @test begin
        result = true
        for i = 1:100
            vxl = (-@SVector(rand(3)), @SVector(rand(3)))
            r₀ = mean(vxl)
            nout = randn(3)
            nout[rand(1:3)] = 0
            result &= volfrac(vxl, SVector{3}(nout), r₀)≈0.5
        end
        result
    end
end  # @testset "quadrangular cylinder 3D"

@testset "quadrangular cylinder 2D" begin
    @test_nowarn @inferred(volfrac((SVector(0,0),SVector(1,1)), SVector(1,2), SVector(0.5,0.5)))
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
    @test_nowarn @inferred(volfrac((-@SVector(rand(3)),@SVector(rand(3))), SVector(0,0,0), @SVector(randn(3))))
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
    vxl = (SVector(0,0,0), SVector(1,1,1))
    r₀ = SVector(0,0,0)
    @test (nout = SVector(1,0,0); volfrac(vxl, nout, r₀) ≈ 0)  # completely outside
    @test (nout = SVector(-1,0,0); volfrac(vxl, nout, r₀) ≈ 1)  # completely inside
    @test (nout = SVector(1,1,0); volfrac(vxl, nout, r₀) ≈ 0)  # completely outside
    @test (nout = SVector(-1,-1,0); volfrac(vxl, nout, r₀) ≈ 1)  # completely inside
    @test (nout = SVector(1,1,1); volfrac(vxl, nout, r₀) ≈ 0)  # completely outside
    @test (nout = SVector(-1,-1,-1); volfrac(vxl, nout, r₀) ≈ 1)  # completely inside
    @test (nout = SVector(-1,1,0); volfrac(vxl, nout, r₀) ≈ 0.5)  # rvol_tricyl()
    @test (nout = SVector(-1,-1,1); volfrac(vxl, nout, r₀) ≈ 5/6)  # rvol_gensect()
    @test (nout = SVector(1,-2,1); volfrac(vxl, nout, r₀) ≈ 0.5)  # rvol_quadsect()
    r₀ = SVector(0.5,0.5,0)
    @test (nout = SVector(-2,1,0); volfrac(vxl, nout, r₀) ≈ 0.5)  # rvol_quadcyl()

    # Test tolerance to floating-point arithmetic.
    vxl = (SVector(0,9.5,0), SVector(1,10.5,1))
    r₀ = SVector(0.5,10.0,0.5)
    nout = SVector(0.0,1/√2,1/√2)
    @test volfrac(vxl, nout, r₀) ≈ 0.5

    # Test rvol_quadsect() for nontrivial cases.
    vxl = (SVector(0,0,0), SVector(1,1,2))
    r₀ = SVector(0.5, 0.5, 0.5)
    @test begin
        result = true
        for i = 1:100
            nout = SVector(randn()/20, randn()/20, 1)
            result &= volfrac(vxl, nout, r₀)≈0.5/2
        end
        result
    end
end  # @testset "general cases"
