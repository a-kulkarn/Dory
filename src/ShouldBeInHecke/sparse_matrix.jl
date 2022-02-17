
import Base: findfirst

# Something in Hecke's Sparse matrix logic fails because this function is missing.
function findfirst(A::AbstractArray{T}, x::T) where T
    return findfirst(isequal(x), A)
end

import Base.size
function size(A::SMat{T} where T, i::Int64)
    return i <= 2 ? size(A)[i] : 1
end

function rows(A::SMat{T} where T)
    return A.rows
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
    sparse_identity(R::T, n::Int64) where T<:Hecke.Ring -> SMat{T}

Given a ring `R`, and an integer `n`, return the `n x n`-identity matrix over the ring R.
"""
function sparse_identity(R::T, n::Int64) where T<:Hecke.Ring
    I = Hecke.sparse_matrix(R)
    for i=1:n
        srow = sparse_row(R, [i], [one(R)])
        push!(I, srow)
    end
    return I
end


