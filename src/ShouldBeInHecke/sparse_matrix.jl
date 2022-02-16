


import Base: findfirst

# Something in Hecke's Sparse matrix logic fails because this function is missing.
function findfirst(A::AbstractArray{T}, x::T) where T
    return findfirst(isequal(x), A)
end

@doc Markdown.doc"""
    setindex!(A::SMat{T}, a::T, i::Int64, j::Int64) -> nothing
Given a sparse matrix $A$, set $A[i,j] = a$.
"""
function Base.setindex!(A::SMat{T}, a::T, i::Int64, j::Int64) where {T<:Hecke.RingElem}

    srow = A[i]    
    p = findfirst(isequal(j), srow.pos)
    if p === nothing
        irange = searchsorted(srow.pos,j)
        splice!(srow.pos, irange, [j])
        splice!(srow.values, irange, [a])
        A.nnz += 1
    else
        srow.values[p] = a
    end
    A[i] = srow
    return
end

@doc Markdown.doc"""
    sparse_identity(R::T, n::Int64) where T<:Hecke.Ring -> Id

Given a ring `R`, and an integer `n`, return the `n x n`-identity matrix over the ring R.

#-------------------

INPUTS:
R         -- type ::Hecke.Ring
col_pivot -- a type, either Val(true) or Val(false), indicating whether column permutations
             should be used to move p-adically large entries to the pivot position.

"""
function sparse_identity(R::T, n::Int64) where T<:Hecke.Ring
    II = Hecke.sparse_matrix(R)
    for i=1:n
        srow = sparse_row( R, Array{Int64,1}([i]), Array{padic,1}([R(1)]) )
        push!(II, srow)
    end
    return II
end


