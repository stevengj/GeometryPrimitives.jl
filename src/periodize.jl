export periodize

periodize(o::Shape{N}, A::AbstractMatrix, ∆range::Shape{N}) where {N} = periodize(o, SMatrix{N,N}(A), ∆range)

function periodize(o::S,  # shape to periodize
                   A::SMatrix{N,N,<:Real,L},  # columns are primitive vectors of lattice
                   ∆range::Shape{N,L}  # range of shift
                  ) where {N,L,S<:Shape{N,L}}
    ∆bound = bounds(∆range)  # bounds of shift range
    nmax_fl = SVector(ntuple(k->-Inf,Val{N}))
    nmin_fl = SVector(ntuple(k->Inf,Val{N}))
    cr = CartesianRange(ntuple(k->2, Val{N}))  # (2,2,2) for N = 3
    for ci = cr
        cnr = map((∆b1,∆b2,i)-> i==1 ? ∆b1 : ∆b2, ∆bound[1], ∆bound[2], SVector(ci.I))  # SVector: corner of ∆bound (ci.I[k] = 1 or 2)
        n_fl = A \ cnr  # SVector: A * n_fl = corner
        nmax_fl = max.(nmax_fl, n_fl)  # SVector
        nmin_fl = min.(nmin_fl, n_fl)  # SVector
    end

    r = 1.0 - Base.rtoldefault(Float64)  # scale factor to multiply to exclude points on boundary
    # Below, use `r.*` for fusion after StaticArrays Issue 390 is resolved.
    nmax = floor.(Int, r * nmax_fl)  # SVector
    nmin = ceil.(Int, r * nmin_fl)  # SVector

    nshape = prod((nmax.-nmin).+1)  # tentative number of shapes
    shp_array = Vector{S}(nshape)  # preallocation (will be trimmed later)

    ishape = 0
    nrange = map((n1,n2)->n1:n2, nmin.data, nmax.data)  # e.g., (1:10, 1:10, 1:10) for N = 3
    for n = CartesianRange(nrange)
		∆ = A * SVector(n.I)  # lattice vector
		if ∆ ∈ ∆range
			shape = translate(o, ∆)
			ishape += 1
			shp_array[ishape] = shape
		end
    end

    resize!(shp_array, ishape)  # return resized shp_array
end
