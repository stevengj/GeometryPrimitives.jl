export SectoralPrism

const SectoralPrism = Prism{Sector}

# Below, if we called SectoralPrism(c, ...) in the function body, it would call the inner
# constructor Prism{Sector{Nothing}}(c, ...) because SectoralPrism = Prism{Sector{Nothing}},
# which is not what we want.
# To call the outer constructor of Prism, we should call Prism(c, ...) instead of SectoralPrism(c, ...).
SectoralPrism(c::SReal{3},
              r::Real,
              ϕ::Real,
              ∆ϕ::Real,
              h::Real=Inf,
              a::SReal{3}=SVec(0.0,0.0,1.0)
              ) where {K} =
    (â = normalize(a); Prism(c, Sector(SVec(0.0,0.0),r,ϕ,∆ϕ), h, [orthoaxes(â)... â]))

SectoralPrism(c::AbsVecReal,  # center of prism
              r::Real,  # radius of sectoral base
              ϕ::Real,  # start angle of sectoral base: 0 ≤ ϕₛ < 2π  (2π excluded)
              ∆ϕ::Real,  # end angle sectoral base: 0 ≤ ϕₑ-ϕₛ ≤ 2π
              h::Real=Inf,  # height of prism
              a::AbsVecReal=[0.0,0.0,1.0]) =  # axis direction of prism
    SectoralPrism(SVec{3}(c), r, ϕ, ∆ϕ, h, SVec{3}(a))

function bounds_ctrcut(s::SectoralPrism)
    ax = s.p'  # prism axes: columns are not only unit vectors, but also orthogonal
    if ax ≈ I  # prism axes are aligned with Cartesian directions (this covers most usage)
        l, u = bounds(s.b)  # (SVec{2}, SVec{2})
        a₁₂ = ax[:,SVec(1,2)]  # SMat{3,2}

        return a₁₂*l, a₁₂*u
    else  # prism axes are not aligned with Cartesian directions
        b = s.b
        r = b.r

        el = Ellipsoid(SVec(0.0,0.0,0.0), SVec(r,r,0.0), ax)  # center is set at origin to return bounds with respect to prism center
        bp = boundpts(el)  # S²Mat{3}: boundary points; all three columns of b are free of NaN because no column of ax is aligned with Cartesian directions
        bp′ = s.p * bp  # S²Mat{3}: bp in prism coordinates
        bp2′ = bp′[SVec(1,2),:]  # SMat{2,3}: in each column of bp′, third entry is in axis dimension, so must be zero mathematically

        ϕ = atan.(bp2′[2,:], bp2′[1,:])  # SVec{3}: angles of boundary points in base plane
        ϕsym = [ϕ; ϕ.+π]  # SVec{6}: include symmetric points with respect to base center
        ϕall = [ϕsym; SVec(b.ϕ₀-b.∆ϕ2, b.ϕ₀+b.∆ϕ2)]  # SVec{8}: angles of all boundary point candidates, except base center

        bpall′ = r .* [cos.(ϕall) sin.(ϕall) @SVector(zeros(8))]'  # SMat{3,8}: each column is boundary point candidate in prism coordinates
        bpall = ax * bpall′  # SMat{3,8}: each column is boundary point candidate in external coordinates
        bpallc = [bpall SVec(0.0,0.0,0.0)]  # SMat{3,9}: include base center
        ind = abs.(distangle.(ϕsym, b.ϕ₀)) .≤ b.∆ϕ2  # SVec{6}: indices of angles contained in base arc
        indallc = [ind; SVec(true,true,true)]  # SVec{9}: include arc ends and base center

        xs = bpallc[1,:]
        ys = bpallc[2,:]
        zs = bpallc[3,:]

        # It might be possible to implement the following more efficiently using reduce().
        xmin = ymin = zmin = Inf
        xmax = ymax = zmax = -Inf
        for i = 1:9
            if indallc[i]
                xmin = min(xmin, xs[i])
                ymin = min(ymin, ys[i])
                zmin = min(zmin, zs[i])

                xmax = max(xmax, xs[i])
                ymax = max(ymax, ys[i])
                zmax = max(zmax, zs[i])
            end
        end

        return (SVec(xmin,ymin,zmin), SVec(xmax,ymax,zmax))
    end
end
