export periodize

periodize(shp::Shape{N}, A::AbstractMatrix{<:Real}, ∆range::Shape{N}) where {N} = periodize(shp, SMatrix{N,N}(A), ∆range)

function periodize(shp::S,  # shape to periodize
                   A::SMatrix{N,N,<:Real},  # columns: primitive vectors of Bravais lattice
                   ∆range::Shape{N}  # range of translation vectors; boundaries included
				   ) where {N,S<:Shape{N}}
    # R = n₁a₁ + n₂a₂ + n₃a₃ is a translation vector for A = [a₁ a₂ a₃].  Find the ranges of
    # n₁, n₂, n₃ to test the inclusion of R in the Shape ∆range.
    #
    # The minimum and maximum nᵢ's to examine can be found by findfirstg nᵢ's at the corners of
    # the bounding box of ∆range.  This can be easily seen by unitary transformation of aᵢ's
    # into the Cartesian directions.  The transformation transforms the lattice planes into
    # the x-, y-, z-normal planes, and the bounding box into a parallelpiped.  The nᵢ's at
    # the corners of the original bounding box corresponds to the Cartesian indices of the
    # corners of the parallelpiped after transformation, and the minimum and maximum nᵢ's
    # are the Cartesian indices of the corners of the bounding box of this post-transform
    # parallelpiped.  Because the post-transform parallelpiped is enclosed by its bounding
    # box, all the lattice points within the post-transform parallelpiped are also within
    # the bounding box of the parallelpiped, whose corner indices are the minimum and
    # maximum nᵢ's.
    ∆bound = bounds(∆range)  # (SVector{3}, SVector{3}): bounds of translation range
    nmax_fl = SVector(ntuple(k->-Inf, Val(N)))  # [-Inf, -Inf, -Inf] for N = 3
    nmin_fl = SVector(ntuple(k->Inf, Val(N)))  # [Inf, Inf, Inf] for N = 3
    rcart = CartesianIndices(ntuple(k->2, Val(N)))  # (2,2,2) for N = 3
    for icart = rcart  # icart: CartesianIndex (see, e.g., https://julialang.org/blog/2016/02/iteration)
        cnr = map((∆b1,∆b2,i)-> i==1 ? ∆b1 : ∆b2, ∆bound[1], ∆bound[2], SVector(icart.I))  # SVector: corner of ∆bound (note icart.I[k] = 1 or 2)
        n_fl = A \ cnr  # SVector: A * n_fl = corner
        nmax_fl = max.(nmax_fl, n_fl)  # SVector
        nmin_fl = min.(nmin_fl, n_fl)  # SVector
    end

    # Find the integral nᵢ's.
    nmax = ceil.(Int, nmax_fl)  # SVector
    nmin = floor.(Int, nmin_fl)  # SVector

    # For nᵢ's within the range calculated above, check if R = n₁a₁ + n₂a₂ + n₃a₃ is within
    # the Shape ∆range.  If it is, then create a shape transformed by R.
    nshape = prod((nmax.-nmin).+1)  # tentative number of shapes
    shp_array = Vector{S}(undef, nshape)  # preallocation (will be truncated later)
    ishape = 0
    nrange = map((n1,n2)->n1:n2, nmin.data, nmax.data)  # e.g., (1:10, 1:10, 1:10) for N = 3
    for n = CartesianIndices(nrange)  # n: CartesianIndex
		∆ = A * SVector(n.I)  # lattice vector R
		if ∆ ∈ ∆range
			shape = translate(shp, ∆)
			ishape += 1
			shp_array[ishape] = shape
		end
    end

    resize!(shp_array, ishape)  # return resized shp_array
end
