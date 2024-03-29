##############################################################################################
#                                                                                            #
#                               p-adic linear algebra                                        #
#                                                                                            #
##############################################################################################

##############################################################################################
#
#    Types
#
##############################################################################################

struct QRNonArchimedeanPivoted{T <: DiscreteValuedFieldElem}
    Q::Hecke.Generic.MatElem{T}
    R::Hecke.Generic.MatElem{T}
    p::Array{Int64,1}
    q::Array{Int64,1}
end

struct QRNonArchimedeanSparsePivoted{T <: DiscreteValuedFieldElem}
    Q::Hecke.SMat{T}
    R::Hecke.SMat{T}
    p::Array{Int64,1}
    q::Array{Int64,1}
end


##############################################################################################
#
#    Norms and normalization
#
##############################################################################################

function norm(A::Hecke.Generic.MatElem{T}) where T <: DiscreteValuedFieldElem
    return maximum(abs, A)
end

function norm_valuation(A::Hecke.Generic.MatElem{T}) where T <: DiscreteValuedFieldElem
    return minimum(valuation, A)
end

function normalize_matrix(A)

    Qp = base_ring(A)
    vals_of_A = valuation.(A)
    min_val = minimum(vals_of_A)

    iexp = -Integer(min_val)
    scale_factor = uniformizer(Qp)^iexp
    return scale_factor * A, scale_factor
end


##############################################################################################
#
#    Iwasawa decomposition (QR)
#
##############################################################################################

# Compute a factorization of a padic matrix A = QR, where R is
# upper triangular (Borel) and Q lies in the maximally compact
# subgroup of SL(Qp) = SL(Zp).
#
# It turns out that the algorithm to do this is just the LU factorization
# with pivoting.

@doc Markdown.doc"""
    padic_qr(A :: Hecke.Generic.MatElem{<:DiscreteValuedFieldElem} ; col_pivot :: Union{Val{true}, Val{false}}) -> F :: QRDiscreteValuedFieldElemPivoted

Compute the p-adic QR factorization of `A`. More precisely, compute matrices `Q`,`R`, and an arrays `p`, `q` such that 

    A[F.p,F.q] = F.Q * F.R

If `col_pivot=Val(false)`, then `F.q = [1,2,...,size(A,2)]`.

#-------------------

INPUTS:

`A`         -- a matrix over Qp

`col_pivot` -- a type, either `Val(true)` or `Val(false)`, indicating whether column permutations
             should be used to move p-adically large entries to the pivot position.

"""
function padic_qr(A::Hecke.Generic.MatElem{T} where T<:DiscreteValuedFieldElem;
                  col_pivot=Val(false) :: Union{Val{true},Val{false}},
                  hessenberg=Val(false) :: Union{Val{true},Val{false}})

    # Check parameters
    if hessenberg == col_pivot == Val(true)
        throw(IncompatibleOptionsError("Cannot use both `col_pivot` and `hessenberg` options in padic_qr."))
    end
    
    # Set constants
    R = base_ring(A)
    T = elem_type(R)
    n = size(A,1)::Int64
    m = size(A,2)::Int64
    basezero = zero(R)
    
    L = unaliased_identity_matrix(R, n)
    Lent = L.entries::Array{T, 2}  
    Umat = deepcopy(A)
    U = Umat.entries
    
    P = Array(1:n)
    Pcol = Array(1:m)

    # We cache the maximum value of the matrix at each step, so we save an iteration pass
    # through the matrix.
    val_list = float64_valuation.(U)
    min_val, min_val_index = findmin(val_list);
    
    # Allocate specific working memory for multiplications.
    container_for_swap    = R()
    container_for_product = R()
    container_for_div     = R()
    
    # Allocate a 2-element array to hold the index of the maximum valuation.
    min_val_index_mut = [x for x in min_val_index.I]
    
    for k=1:(min(n,m)::Int64)

        if col_pivot == Val(true)
            col_index = min_val_index_mut[2]
            if col_index != k
                # interchange columns m and k in U
                for r=1:n
                    U[r,k], U[r,col_index] = U[r,col_index], U[r,k]
                end
                
                # interchange entries m and k in Pcol
                Pcol[k], Pcol[col_index] = Pcol[col_index], Pcol[k]
            end
        end

        # Determing the row to stop row operations based on the matrix shape.
        if hessenberg == Val(true)
            row_stop = min(k+1, n)
        else
            row_stop = n
        end
        
        val_list = float64_valuation.(U[k:row_stop,k])
        minn, row_pivot_index = findmin(val_list);
        if minn == Inf continue end

        row_pivot_index=row_pivot_index+k-1;
        if row_pivot_index!=k

            # interchange rows `row_pivot_index` and `k` in U
            for r=1:m
                U[k,r], U[row_pivot_index,r] = U[row_pivot_index,r], U[k,r]
            end               
            
            # interchange entries `row_pivot_index` and k in P
            P[k],P[row_pivot_index] = P[row_pivot_index],P[k]

            # swap columns corresponding to the row operations already done.
            swap_prefix_of_row!(Lent, k, row_pivot_index)
        end

        # Reset min_valuation for selecting new pivot.
        min_val = Inf

        # Note to self: with Julia, the optimal thing to do is split up the row operations
        # and write a for loop.
        # The entries left of the k-th column are zero, so skip these.
        # Cache the values of L[j,k] first.
        #
        if iszero(U[k,k])
            # If col_pivot == true, then we don't need to perform further column swaps
            # in this case, since the element of largest valuation in the lower-right
            # block is zero. In fact, no further operations need to be performed.
            continue
        end 

        # The use of the inversion command preserves relative precision. By row-pivoting,
        # the extra powers of p cancel to give the correct leading term.
        # the "lost" digits of precision for L[j,k] can simply be set to 0.
        container_for_inv = inv(U[k,k]) 
        
        for j=k+1:row_stop
            Hecke.mul!(L[j,k],U[j,k], container_for_inv)
            L[j,k] = setprecision!(L[j,k], precision(parent(L[j,k])))  # L[j,k] is really an integer.
        end

        # zero the subdiagonal the of column.
        for j = k+1:row_stop
            zero!(U[j,k])
        end

        # Perform the row operations.
        for r=k+1:m
            for j=k+1:row_stop
                # Compute U[j,r] = U[j,r] - L[j,k]*U[k,r]                
                Hecke.mul!(container_for_product, L[j,k], U[k,r])
                U[j,r] = _unsafe_minus!(U[j,r], container_for_product)
                
                # Update the smallest valuation element
                if float64_valuation(U[j,r]) < min_val
                    min_val = float64_valuation(U[j,r])
                    min_val_index_mut[1] = j
                    min_val_index_mut[2] = r
                end
            end
        end
    end

    return QRNonArchimedeanPivoted(L,Umat,P,Pcol)
end

function swap_prefix_of_row!(Lent, k::Int64, i::Int64)
    # The index of the diagonal point is (k,k)
    for r=1:(k-1)
        container_for_swap = Lent[k,r]
        Lent[k,r] = Lent[i,r] 
        Lent[i,r] = container_for_swap
    end
    return
end


#######################################################################################
#
# Sparse QR-algorithm
#
#######################################################################################


# The index of the diagonal point is (k,k)
function swap_prefix_of_column!(L, diagonal_index::Int64, i::Int64)
    k = diagonal_index
    for r = 1:k-1
        L[r,k], L[r,i] = L[r,i], L[r,k]
    end
    return
end

function padic_qr(A::Hecke.SMat{T} where T<:DiscreteValuedFieldElem;
                  col_pivot=Val(false) :: Union{Val{true},Val{false}})
    
    # Set constants
    Qp = base_ring(A)
    T = elem_type(Qp)
    n = size(A,1)::Int64
    m = size(A,2)::Int64
    
    # We store the ***transpose*** of L as a sparse matrix, and flip everything at the end.
    # Allocate the rows of Ltrans ahead of time.
    Ltrans = sparse_identity(Qp, n)
    
    U = deepcopy(A)    
    P = Array(1:n)
    Pcol = Array(1:m)

    # Function to pivot and return the list of rows with a non-zero entry at index k.
    # in the subdiagonal.
    function pivot_and_select_row_indices(U, Ltrans, k, piv)

        # Scan through the matrix to check if a column swap is needed.
        if col_pivot==Val(true)
            minn = Inf
            mindex = piv
            for j=k:n
                srow = U[j]
                if isempty(srow) break end
                
                rowmin, rowmindex = findmin(valuation.(srow.values))
                if rowmin < minn
                    minn = rowmin
                    mindex = srow.pos[rowmindex]
                end
            end

            if mindex != piv
                Hecke.swap_cols!(U,piv,mindex)
                Pcol[piv], Pcol[mindex] = Pcol[mindex], Pcol[piv]
            end
        end

        # Scan through the matrix to check if a rowswap is needed.
        valuation_index_pairs = [(valuation(U[j, piv]), j) for j=k:n if !iszero(U[j, piv])]
            
        if !isempty(valuation_index_pairs)

            row_pivot_index = last(minimum(valuation_index_pairs))
            rows_with_entry_at_piv = last.(valuation_index_pairs)
            
            if row_pivot_index!=k
                Hecke.swap_rows!(U,k,row_pivot_index)
                P[k],P[row_pivot_index] = P[row_pivot_index],P[k]

                # swap columns corresponding to the row operations already done.
                # Do not swap the diagonal elements.
                swap_prefix_of_column!(Ltrans, k, row_pivot_index)
            end

            # Remove k or leading zeros from the list of rows to iterate over.
            # (Leading zeros can be introduced by a swap.)
            filter!(j->(j!=k && !iszero(U[j,piv])), rows_with_entry_at_piv)
        else        
            return last.(valuation_index_pairs)
        end
    end ### END PIVOT FUNCTION

    # Initialize the shift index to determine the critial pivot location.
    shift=0
    k=1

    while k <= (min(n,m)::Int64) && k+shift <= m::Int64

        ### Pivot logic ###
        
        # set the pivot column and determine the rows to apply elimination.
        piv = k + shift
        rows_with_entry_at_piv = pivot_and_select_row_indices(U, Ltrans, k, piv)

        # If everything is zero, shift the algorithm to operate on the
        # right rectangular submatrix window.
        if isempty(rows_with_entry_at_piv) shift+=1; continue; end   

        
        ### Row operations loop ###
        container_for_inv = inv(U[k,piv])
        
        for j in rows_with_entry_at_piv
                        
            # The "lost" digits of precision for L[j,k] can simply be set to 0.
            # as L[j,k] is really an integer.
            Ltrans[k,j] = U[j,piv]*container_for_inv

            if precision(Ltrans[k,j]) < precision(Qp)
                Ltrans[k,j] = setprecision!(Ltrans[k,j], precision(Qp))
            end
            
            if Ltrans[k,j] != 0
                Hecke.add_scaled_row!(U, k, j, -Ltrans[k,j])
            elseif valuation(Ltrans[k,j]) < precision(Qp)
                error("Problem related to Hecke's `add_scaled_row` function encountered.")
            end
        end

        #Update loop counter.
        k += 1        
    end
    
    # The test is actually quite expensive, but we keep it for now.
    @vtime :local_QR @assert iszero(matrix(A)[P,Pcol] - transpose(matrix(Ltrans))*matrix(U))

    return QRNonArchimedeanSparsePivoted(transpose(Ltrans), U, P, Pcol)
end


########################################################################################

# IMPORTANT!
# We deviate slightly from LinearAlgebra's SVD structure by putting a diagonal matrix for S.
struct SVDNonArchimedean{T<:DiscreteValuedFieldElem}
    U::Hecke.Generic.MatElem{T}
    S::Hecke.Generic.MatElem{T}
    Vt::Hecke.Generic.MatElem{T}
    p::Array{Int64,1}
    q::Array{Int64,1}
end

# A padic analogue for svd
@doc Markdown.doc"""
    svd(A :: Hecke.Generic.Matelem{<:DiscreteValuedFieldElem}) -> SVDNonArchimedean

  Compute the singular value decomposition (SVD) of A and return an SVDNonArchimedean object.

  A pAdic singular value decomposition is a factorization of the form

  A = U * S * Vt

  where U,Vt are matrices in GL_n(Zp). For efficiency reasons, we give a factorization
  as well as two arrays `F.p` and `F.q` such that

  A[p, q] = U*S*Vt

  where `U` is lower triangular with ones on the diagonal and `Vt` is upper triangular
  with ones on the diagonal. The singular values in S are sorted in descending order of padic 
  absoute value.

  The return types are
     U, S, Vt :: Hecke.Generic.MatElem{<:DiscreteValuedFieldElem}
     p,q      :: Array{Int64,1}

"""
function svd(A::Hecke.Generic.MatElem{T}) where T <: DiscreteValuedFieldElem

    F = padic_qr(A, col_pivot=Val(true))
    G = padic_qr(transpose(F.R))

    @assert G.p == [i for i=1:length(G.p)]

    U = deepcopy(F.Q)
    S = transpose(G.R)
    Vt= transpose(G.Q)
    
    @assert iszero(A[F.p,F.q] - U*S*Vt)

    return SVDNonArchimedean(U, S, Vt, F.p, F.q)
end

# stable version of rank for padic matrices.
@doc Markdown.doc"""
    rank(A::Hecke.Generic.MatElem{<:DiscreteValuedFieldElem})

Compute the rank of a padic matrix by counting how many singular values satisfy `iszero(a)`.
"""
function rank(A::Hecke.MatElem{T}) where T <: DiscreteValuedFieldElem
    n = nrows(A)
    m = ncols(A)
    F = padic_qr(A)

    rank=0
    for i=1:min(n,m)
        if !iszero(F.R[i,i])
            rank += 1
        end
    end
    return rank
end

# Returns the p-adic singular values of a matrix
@doc Markdown.doc"""
    singular_values(A::Hecke.MatElem{<:DiscreteValuedFieldElem}) -> Array{<:DiscreteValuedFieldElem, 1}

Returns the list of diagonal elements in the singular value decomposition of the matrix `A`.
"""
function singular_values(A::Hecke.MatElem{<:DiscreteValuedFieldElem})
    F = padic_qr(A, col_pivot=Val(true))
    return [F.R[i,i] for i=1:minimum(size(A))]
end


@doc Markdown.doc"""
    nullspace(A::Hecke.MatElem{DiscreteValuedFieldElem}) -> (nu, N)

Computes the nullspace of a padic matrix `A`. The dimension of the nullspace is dependent
on the number of singular values of `A` for which `iszero(a)` is true.

nu -- An Int64, which is the dimension of the nullspace.
N  -- A matrix whose columns generate the nullspace of A. Type ::Hecke.MatElem{DiscreteValuedFieldElem}. 
"""
function nullspace(A::Hecke.MatElem{<:DiscreteValuedFieldElem})

    m = nrows(A)
    n = ncols(A)
    F = padic_qr(transpose(A), col_pivot=Val(true))

    col_list = Array{Int64,1}()
    for i=1:min(n,m)
        if iszero(F.R[i,:])
            push!(col_list, i)
        end
    end

    Pinv = invperm(F.p)   
    
    Q = F.Q
    inv_unit_lower_triangular!(Q)
    Qinvt = transpose(Q)[Pinv,:]
    
    return length(col_list) + max(0,n-m), hcat(Qinvt[:, col_list], Qinvt[:,(m+1):n])
end

function nullspace(A::Hecke.SMat{<:DiscreteValuedFieldElem})

    m = nrows(A)
    n = ncols(A)
    F = padic_qr(transpose(A), col_pivot=Val(true))

    # Sparse elimination sorts zero singular values to the bottom.
    # However, without column pivoting there is potential precision loss.

    
    # This function is really bad from a precision analysis standpoint. It isn't clear
    # how to deal with the insane precision gains we get from sparse elimination.
    function bad_iszero(srow)
        if iszero(srow)
            return true
        end
        return minimum(valuation.(srow.values)) >= precision(base_ring(A))
    end
    
    col_list = Array{Int64,1}()
    for i=1:min(n,m)
        if bad_iszero(F.R[i])
            push!(col_list, i)
        end
    end

    Pinv = invperm(F.p)   
    
    Q = F.Q
    
    badQ = matrix(Q) # The "matrix" call makes things dense. This is not ideal.
    inv_unit_lower_triangular!(badQ)   # It is probably a good idea to have a specialized QR method
                                       # that computes the inverse of Q directly, instead of Q.
    Qinvt = transpose(badQ)[Pinv,:]
    
    return length(col_list) + max(0,n-m), hcat(Qinvt[:, col_list], Qinvt[:,(m+1):n])
end


# stable version of inverse for p-adic matrices
function inv(A::Hecke.MatElem{<:DiscreteValuedFieldElem})
    check_square(A)
    id = unaliased_identity_matrix(base_ring(A), size(A,2))
    return rectangular_solve(A, id)
end

function inv_unit_lower_triangular!(L::Hecke.Generic.MatElem{T} where T)

    m = size(L,1)::Int64
    n = size(L,2)::Int64    

    Qp = parent(L[1,1])
    container_for_mul = Qp()
    container_for_result = Qp()
    
    for k = 1:n
        for i = k+1:n
            container_for_result=zero(Qp)
            for r=k:i-1
                Hecke.mul!(container_for_mul, L[i,r], L[r,k])
                addeq!(container_for_result,  container_for_mul)
            end
            L[i,k] = -container_for_result
        end
    end

    return L
end

@doc Markdown.doc"""
    inv_unit_lower_triangular(A::Hecke.MatElem{<:DiscreteValuedFieldElem}) -> Hecke.MatElem{<:DiscreteValuedFieldElem}

Matrix inverse, specialized to invert a lower triangular matrix with ones on the diagonal.
"""
function inv_unit_lower_triangular(L)
    L2 = deepcopy(L)
    return inv_unit_lower_triangular!(L2)
end

@doc Markdown.doc"""
    rectangular_solve(A::Hecke.MatElem{T}, b_input::Hecke.MatElem{T}; stable::Bool=false) where T<:DiscreteValuedFieldElem
                                                                        -> (nu :: Int64,N::Hecke.MatElem{DiscreteValuedFieldElem})

Solves the linear system A*N = b. The output `nu` is the dimension of the nullspace. Parameter `stable` determines whether `padic_qr` or `svd` method is used. Default is qr (for speed).

WARNINGS:
If `A,b_input` have different precisions, maximal precision output is not guarenteed.
Underdetermined solve not implemented.
"""
function rectangular_solve(A::Hecke.MatElem{<:DiscreteValuedFieldElem}, b::Hecke.MatElem{<:DiscreteValuedFieldElem}; stable::Bool=false)
    return rectangular_solve!(A, deepcopy(b), stable=stable)
end

@doc Markdown.doc"""
    rectangular_solve!(A::Hecke.MatElem{<:DiscreteValuedFieldElem}, b_input::Hecke.MatElem{<:DiscreteValuedFieldElem}; stable::Bool=false)
                                                                        --> (nu :: Int64,N::Hecke.MatElem{padic})

Solves the linear system A*N = b inplace. The output `nu` is the dimension of the nullspace. Parameter `stable` determines whether `padic_qr` or `svd` method is used. Default is qr (for speed).

WARNINGS:
If `A,b_input` have different precisions, maximal precision output is not guarenteed.
Underdetermined solve not implemented.
"""
function rectangular_solve!(A::Hecke.MatElem{T}, b::Hecke.MatElem{T}; stable::Bool=false) where T <: DiscreteValuedFieldElem
    if stable
        return _svd_rectangular_solve(A::Hecke.MatElem{T}, b::Hecke.MatElem{T})
    else
        return _lu_rectangular_solve(A::Hecke.MatElem{T}, b::Hecke.MatElem{T})
    end
end

# Specialization to lu-solve
function _lu_rectangular_solve(A::Hecke.MatElem{T}, b_input::Hecke.MatElem{T}) where T <: DiscreteValuedFieldElem

    m = nrows(A)
    n = ncols(A)
    if nrows(b_input) != m
        throw(DimensionMismatch("`A` and `b` must have the same number of rows."))
    end
    b = deepcopy(b_input)

    if m < n
        error("System is underdetermined. Use `underdetermined_solve` instead.")
    end

    F = padic_qr(A)
    b = b[F.p,:]

    # forward substitution, all diag entries are scaled to 1
    for i in 1:m
        for j in 1:(i-1)
            b[i,:] = b[i,:] - b[j,:]* F.Q[i,j]
        end
    end

    # consistency check for overdetermined systems
    if m > n
        for i in (n+1):m
            for j in 1:ncols(b)
                if !iszero(b[i, j])
                    throw(InconsistentSystemError(A, b_input))
                end
            end
        end
    end
    b = b[1:n, :]   # truncate zero rows if consistent

    # backward substitution
    for i in n:-1:1
        for j in (i+1):n
            b[i,:] = b[i,:] - b[j,:]*F.R[i,j]
        end
        #scale = A[i, i]
        
        if !iszero(b[i,:]) && iszero(F.R[i,i])
            throw(InconsistentSystemError(A, b_input))
            
        elseif !iszero(F.R[i,i])
            b[i,:] *= inv(F.R[i,i])
        end
    end

    return b
end

# Specialization to svd-solve
function _svd_rectangular_solve(A::Hecke.MatElem{T}, b_input::Hecke.MatElem{T}) where T <: DiscreteValuedFieldElem

    m = nrows(A)
    n = ncols(A)
    if nrows(b_input) != m
        throw(DimensionMismatch("`A` and `b` must have the same number of rows."))
    end
    b = deepcopy(b_input)

    if m < n
        error("System is underdetermined. Use `underdetermined_solve` instead.")
    end

    F = svd(A)
    b = b[F.p,:]

    # forward substitution, all diag entries are scaled to 1
    for i in 1:m
        for j in 1:(i-1)
            b[i,:] = b[i,:] - b[j,:]* F.U[i,j]
        end
    end

    # consistency check for overdetermined systems
    if m > n
        for i in (n+1):m
            for j in 1:ncols(b)
                if !iszero(b[i, j])
                    throw(InconsistentSystemError(A, b))
                end
            end
        end
    end
    b = b[1:n, :]   # truncate zero rows if consistent

    # Scaling step
    for i in 1:n
        if !iszero(b[i,:]) && iszero(F.S[i,i])
            throw(InconsistentSystemError(A,b))

        elseif !iszero(F.S[i,i])
            b[i,:] *= inv(F.S[i,i])
        end
    end

    # backward substitution
    for i in n:-1:1
        for j in (i+1):n
            b[i,:] = b[i,:] - b[j,:]*F.Vt[i,j]
        end
    end

    return b[F.q,:]
end
    
function underdetermined_solve()
    throw(NotImplemented)
    return
end


##############################################################################################
#
# Eigenvector iteration methods.
#
##############################################################################################


#################################################
#  Basic inverse iteration
#################################################

# REMARK: Invese iteration is somehow theoretically doomed to result in precision loss.
#         Inverting (A - λ * I) is usually ill-conditioned.
#
#         The loss of precision is more related to the issue of the maximum size of
#         Av - λv, given that Bv = μv + O(p^N). 

@doc Markdown.doc"""
    function inverse_iteration!(A, shift, V)
    
Solve for an eigenvector using inverse iteration.
Note that the algorithm will not converge to a particular vector in general, but the norm of
A*w - λ*w converges to zero. Here, λ is the unique eigenvalue closest to `shift`, (if it is unique).

"""
function inverse_iteration!(A, shift, V)

    # Note: If A is not known to precision at least one, really bad things happen.
    Qp = base_ring(A)
    N  = precision(Qp)
    In = unaliased_identity_matrix(Qp, size(A,1))
    B  = A - shift * In
    
    if rank(B) < ncols(B)
        @vprint :local_inverse_iteration "Value `shift` is exact eigenvalue: shift = $shift"
        return [nullspace(B)[2]], [shift]
    end

    # Compute a pseudo-inverse for B
    F = svd(B)
    maxdval = maximum(valuation, diagonal(F.S))

    Dpseudoinv = let
        scale = uniformizer(Qp)^maxdval
        diag = diagonal(F.S)
        
        newd = fill(zero(Qp), size(F.S, 1))
        for i in 1:size(F.S,1)
            a = diag[i]
            b = inv(a) * scale
            b = setprecision!(b, N)
            newd[i] = b
        end
        
        diagonal_matrix(newd)
    end

    # Set the precision on V based on the condition number for inversion.
    for j=1:size(V,2)
        for i=1:size(V,1)
            V[i,j] = setprecision!(V[i,j], N - maxdval)
        end
    end
    
    pow = let
        Y = inv(F.Vt) * Dpseudoinv * inv(F.U)
        Y[invperm(F.q), F.p]
    end

    ### Iteration Loop ###
    for i=1:N * size(A, 1)
        Vprev = V
        V, _ = normalize_matrix(pow * V)
    end

    # Test for convergence and calculate the restricted transform on the eigenspace.
    X = try
        # !iszero(V[3,1]) && @info " " V/V[3,1] - pow * V / (pow * V)[3,1]
        rectangular_solve(V, A*V, stable=true)        
    catch e
        throw(ConvergenceFailureError("Error in inverse iteration. Likely a stability issue."))
    end

    nu = trace(X) // size(X,2)

    # TODO: It is possible that we can diagonalize X and refine the eigenspaces.
    #       However, the theoretical details haven't been worked out here.
    #
    # First step: check:
    #     Y  = X - nu * unaliased_identity_matrix(Qp,size(X,2))
    #     iszero(Y)
    
    return [V],[nu]
end


@doc Markdown.doc"""
    inverse_iteration(A::Hecke.MatElem{<:DiscreteValuedFieldElem}, shift, v ::Hecke.MatElem{<:DiscreteValuedFieldElem}) -> Hecke.MatElem{<:DiscreteValuedFieldElem}

Iterate `v = (A-shift*I)^(-1) * v`. The inverse is cached at the beginning of the computation. The columns of the entry `v` define a subspace.

If subspace iteration does not satisfy `A*v ⊆ v`, an error is raised.
"""
function inverse_iteration(A, shift, v)
    w = deepcopy(v)
    wlist, nulist = inverse_iteration!(A, shift, w)
    return wlist, nulist
end

@doc Markdown.doc"""
    inverse_iteration_decomposition(A::Hecke.MatElem{<:DiscreteValuedFieldElem}, Amp::Hecke.MatElem{nmod_mat}) -> values, spaces

Return types.
    values :: Array{<:DiscreteValuedFieldElem,1}
    spaces :: Array{ Hecke.MatElem{<:DiscreteValuedFieldElem}, 1}

Internal function. Compute an invariant subspace decomposition of `A` using its reduction mod-p `Amp`.
Makes one call to `inverse_iteration` for each eigenvalue of `Amp`. 
"""
function inverse_iteration_decomposition(A, Amp)

    Qp = A.base_ring
    E = eigspaces(Amp)

    values_lift = fill(zero(Qp), 0)
    spaces_lift = fill(zero(parent(A)), 0)

    for i in 1:length(E.values)

        # Approximate input data
        appx_eval = Qp(lift(E.values[i]))
        appx_espace =  change_base_ring(Qp, lift(E.spaces[i]))

        # Apply inverse iteration step.
        wlist, nulist = inverse_iteration(A, appx_eval, appx_espace)

        # Append refined data to the main list.
        values_lift = vcat(values_lift, nulist)
        spaces_lift = vcat(spaces_lift,  wlist)
    end

    return values_lift, spaces_lift
end


###############################################################################
#
#   Power iteration decomposition
#
###############################################################################


function power_iteration_decomposition(A, Amp)

    Qp = base_ring(A)
    N = precision(Qp)
    E = eigspaces(Amp)

    restricted_maps = Array{typeof(fill(zero(Qp), 0)), 1}()
    spaces_lift = Array{typeof(fill(zero(parent(A)), 0)), 1}()

    roots_and_mults = roots_with_multiplicities(Hecke.charpoly(Amp))

    if length(E.values) > 0
        M = maximum([a[2] for a in roots_and_mults])
    end
    
    for i in 1:length(E.values)
       
        # Approximate input data
        appx_eval = Qp(lift(E.values[i]))
        appx_espace =  change_base_ring(Qp, lift(E.spaces[i]))

        # Apply power iteration step.

        B = A - appx_eval * unaliased_identity_matrix(Qp,size(A,1))

        for j=1:ceil(log2(M*N))
            B = B^2
        end

        nu,V = nullspace(B)

        if nu == 0
            display(valuation.(singular_values(A)))
            display(valuation.(singular_values(B)))
            error("Matrix 'B' is invertible. Iteration did not converge.")
        end
        
        X = rectangular_solve(V, A*V, stable=true)
        
        # Append refined data to the main list.
        restricted_maps = vcat(restricted_maps, [X])
        spaces_lift = vcat(spaces_lift,  [V])
    end

    return restricted_maps, spaces_lift

end

###############################################################################
#
#   "Classical Algorithm"
#
###############################################################################

function _eigenspaces_by_classical(A::Generic.MatElem{T}) where T

    # TODO: There will be precision errors since `charpoly` is only computed
    #       at a flat precision.
    f = charpoly(A)
    rts = roots(f)
    n = size(A,2)
    I = unaliased_identity_matrix(base_ring(A), n)
    
    return EigenSpaceDec(base_ring(A), T[rt for rt in rts], Generic.MatElem{T}[nullspace(A - rt*I)[2] for rt in rts])
end

###############################################################################
#
#   Hessenberg form
#
###############################################################################

# Also return a basis by default.
@doc Markdown.doc"""
    hessenberg!(A::Hecke.Generic.Mat{<:DiscreteValuedFieldElem}; basis=Val(true)) --> nothing or B::Hecke.Generic.Mat{T}

Computes the Hessenberg form of `A` inplace. If `basis=Val(true)`, also return the matrix B such that
    AV = VB
"""
function hessenberg!(A::Hecke.Generic.Mat{T} where T <: DiscreteValuedFieldElem; basis=Val(true))
    !issquare(A) && DimensionMismatch("Dimensions don't match in hessenberg")
    R = base_ring(A)
    N = precision(R)
    n = nrows(A)
    u = R()
    t = R()

    if basis == Val(true)
        # Initialize the identity matrix.
        B = parent(A)()
        for i = 1:size(B, 1)
            B[i,i] = one(R)
        end
    end

    for m = 1:n - 2

        val_list = float64_valuation.(A[m+1:n , m])
        minn, row_pivot_index = findmin(val_list);

        if minn==Inf continue end
        i = row_pivot_index[1] + m;            
        
        # Perform a row/column swap to move the pivot to the subdiagonal
        if i > m+1
            for j = m:n
                A[i, j], A[m+1, j] = A[m+1, j], A[i, j]
            end
            for j = 1:n
                A[j, i], A[j, m+1] = A[j, m+1], A[j, i]
            end

            if basis==Val(true)
                for j = 1:n
                    B[i, j], B[m+1, j] = B[m+1, j], B[i, j]
                end
            end
        end

        # cache the inverted pivot.
        h = -inv(A[m+1, m])
        
        # Perform the elimination.
        for i = m + 2:n
            iszero(A[i, m]) && continue
            
            u = Hecke.mul!(u, A[i, m], h)
            u = setprecision!(u, N) # The pivot is always the largest entry.

            # Row operatons
            for j = 1:n
                if j > m
                    t = Hecke.mul!(t, u, A[m+1, j])
                    A[i, j] = addeq!(A[i, j], t)
                end
                    
                if basis==Val(true)
                    t = Hecke.mul!(t, u, B[m+1,j])
                    B[i,j] = addeq!(B[i,j], t)
                end
            end
            u = -u

            # Column eliminations
            for j = 1:n
                t = Hecke.mul!(t, u, A[j, i])
                A[j, m+1] = addeq!(A[j, m+1], t)
            end
            A[i, m] = R()            
        end
    end

    if basis==Val(true)
        return A, B
    else
        return A
    end
end

@doc Markdown.doc"""
    hessenberg(A::Generic.MatrixElem{T}; basis::Union{Val{true}, Val{false}}) where T<:DiscreteValuedFieldElem -> Generic.MatElem{T} [, Generic.MatElem{T}]

Returns the Hessenberg form of A, i.e. an upper Hessenberg matrix
which is similar to A. The upper Hessenberg form has nonzero entries
above and on the diagonal and in the diagonal line immediately below the
diagonal.

The transformation matrix can be returned as the second output. The transformation matrix will always
lie in $GL_{n}(Zp)$.
"""
function hessenberg(A::Hecke.Generic.Mat{T} where T <: DiscreteValuedFieldElem, basis=Val(true))
    check_square(A)
    M = deepcopy(A)
    A, B = hessenberg!(M, basis=Val(true))
    return M, B
end


###############################################################################
#
#   QR-iteration
#
###############################################################################

#################################################
#  QR-iteration helper functions
#################################################

# Default function for selecting a rayleigh shift of a block.
function default_rayleigh_shift(B)
    Qp = base_ring(B)
    m = size(B,2)
    
    value = trace(B)/Qp(m)
    setprecision!(value, precision(Qp))
    return value
end

function default_iter_bound(m::Integer, N)
    if m == 1
        return Int(ceil(log(2,N))) + 5
    else
        return m * N + 2 * m + 3
    end
end

#################################################
#  QR-iteration 
#################################################

@doc Markdown.doc"""
    block_schur_form(A::Hecke.Generic.Mat{<:DiscreteValuedFieldElem}) --> B,V :: Hecke.Generic.Mat{T}

Computes the block schur form `B` of a padic matrix `A`, where the
blocks correspond to the different eigenvalues of `A modulo p`. The outputs satisfy
`VA = BV`

NOTE: 
Presently, `block_shur_form` does not attempt to further refine the blocks recursively. Theoretical
details need to be worked out to make the best practical improvements of the algorithm. 
"""
function block_schur_form(A::Hecke.Generic.Mat{T} where T <: DiscreteValuedFieldElem;
                          shift = default_rayleigh_shift,
                          iter_bound = default_iter_bound)

    Qp = base_ring(A)
    N = precision(Qp)

    iszero(A) && return A, unaliased_identity_matrix(base_ring(A), ncols(A))
    
    # Extract data from the reduction modulo p
    Aint, scale_factor  = normalize_matrix(A)
    Amp   = modp.(Aint)
    chiAp = charpoly(Amp)
    inv_scale_factor = inv(scale_factor)

    
    ####################################################
    # Part 1: computation of the sorted form. (Skipped)
    ####################################################
    
    # The first shift needs to happen before the hessenberg form.
    # TODO: New changes:
    # 1. Figure out the kernel + orthogonal complement of Amp
    # r,p,L,U = lu(Amp)
    
    # 2. Make the basis change. We need to lift the matrices back to Qp
    # invL = inv(L)
    # B = lift(invL) * A * change_base_ring(Qp, L)

    # TODO: also update the permutation.
    # V = invL

    ####################################################
    # Part 2: Main QR iteration.
    ####################################################

    B, V = hessenberg(A)    
    id = unaliased_identity_matrix(Qp, size(B,1))

    # @info "block_schur_form, Hessenberg" precision.(B)
    
    bottom_block_end = size(A,2)
    rts_and_muls = roots_with_multiplicities(chiAp)
    sort!(rts_and_muls, lt=(x,y)->(x[2]<y[2]))
    
    # rts = map(x->x[1], rts_and_muls)
    # m = isempty(rts) ? 0 : maximum(map(x->x[2], rts_and_muls))

    for (rt, m) in rts_and_muls
        
        # Regarding convergence. Some extra time is needed as QR is not always
        # rank revealing.

        num_iters = iter_bound(m, N)
        
        for i in 1:num_iters

            bottom_block_range = bottom_block_end-m+1:bottom_block_end
            Bview = view(B, bottom_block_range, bottom_block_range)

            shift_value = shift(Bview)
            
            # Ensure that the rayleigh shift actually helps convergence.
            if modp(shift_value * scale_factor) != rt
                lambdaI = lift(rt) * id
            else
                lambdaI = shift_value * id
            end

            # QR-step.
            # Note about Julia's syntax. A[:,F.p] = A*inv(P), for a permutation P.
            #
            F = padic_qr(B - lambdaI, hessenberg=Val(true))
            B = F.R[:,F.p] * F.Q + lambdaI
            V = inv_unit_lower_triangular(F.Q)*V[F.p,:]

            # i == 1 && @info "block_schur_form, update" precision.(B)

        end

        bottom_block_end = bottom_block_end - m
    end
    
    return B,V
end


###############################################################################
#
#   Eigenvalues and Eigenvectors
#
###############################################################################

@doc Markdown.doc"""
    eigspaces(A::Hecke.Generic.Mat{<:DiscreteValuedFieldElem}; method=:power)

Compute the eigenvectors of a padic matrix iteratively.

The `method` parameter selects the method to be used to compute the eigenvectors.
The options are:

- :classical
- :power
- :qr
- :schur
- :inverse (Caveat: precision loss guarenteed)

The default is `inverse`, since at the moment this is the one that is implemented.

"""
function eigspaces(A::Hecke.Generic.Mat{<:DiscreteValuedFieldElem}; method=:power)

    check_square(A)
    Qp = base_ring(A)

    if method isa String
        method = Symbol(method)
    end
    
    if iszero(A)        
        return EigenSpaceDec(Qp, [zero(Qp)] , [unaliased_identity_matrix(Qp, size(A,1))])
    elseif size(A,1) == 1
        return EigenSpaceDec(Qp, [A[1,1]] , [unaliased_identity_matrix(Qp, size(A,1))])
    end
    
    if method == :classical
        return _eigenspaces_by_classical(A)
        
    elseif method == :inverse
        return _eigenspaces_by_inverse_iteration(A)
        
    elseif method == :schur || method == :qr
        error("Not Implemented. However, block_schur_form is available to compute the schur form.")
        
    elseif method == :power
        return  _eigenspaces_by_power_iteration(A)
        
    else
        @info " " method
        error("Not Implemented.")
    end
end

function _modp_charpoly_data(A::Hecke.Generic.Mat{T} where T <: DiscreteValuedFieldElem)
    Aint, scale_factor  = normalize_matrix(A)
    Amp   = modp.(Aint)
    chiAp = charpoly(Amp)

    return Aint, Amp, chiAp
end

function _eigenspaces_by_inverse_iteration(A::Hecke.Generic.Mat{T} where T <: DiscreteValuedFieldElem)
    
    # Extract data from the reduction modulo p
    Qp = base_ring(A)
    T = elem_type(Qp)
    Aint, Amp, chiAp = _modp_charpoly_data(A)    
    factors_chiAp = Hecke.factor(chiAp)

        
    if isirreducible(chiAp)
        empty_array = Array{T,1}()
        empty_spaces_array = Array{Hecke.Generic.Mat{T}, 1}()
        
        return EigenSpaceDec(Qp, empty_array , empty_spaces_array)
    end
    
    # FAILSAFE DURING DEVELOPMENT...
    # Fail automatically if there are large invariant subspaces mod p
    if any(e >= 2 for (f,e) in factors_chiAp if degree(f)==1)
        error("Development Failsafe: Not implemented when roots are not squarefree") 
    end

    # Iteration call
    values_lift, spaces_lift = inverse_iteration_decomposition(Aint, Amp)

    return EigenSpaceDec(Qp, values_lift, spaces_lift)
    
end

function _eigenspaces_by_power_iteration(A::Hecke.Generic.Mat{T} where T <: DiscreteValuedFieldElem)

    Qp = base_ring(A)
    T = elem_type(Qp)
    Aint, Amp, chiAp = _modp_charpoly_data(A)
    factors_chiAp = Hecke.factor(chiAp)
    
    empty_array = Array{T,1}()
    empty_spaces_array = Array{Hecke.Generic.Mat{T}, 1}()    
        
    if isirreducible(chiAp)        
        return EigenSpaceDec(Qp, empty_array, empty_spaces_array)
    end

    factor_multiplicities = collect(Base.values(factors_chiAp.fac))

    #println(Amp)    
    #println(factor_multiplicities)
    
    # Check to ensure chiAp is not an n-th power
    if length(factors_chiAp) == 1 && factor_multiplicities[1] == size(A,1)
        return try
            _eigenspaces_by_classical(A)
        catch e
            println()
            println("WARNING: Classical Algorithm not implemented. Skipping this eigenvalue...")
            println()
            EigenSpaceDec(Qp, empty_array , empty_spaces_array)
        end
    end

    # Iteration call
    restricted_maps, invariant_blocks = power_iteration_decomposition(Aint, Amp)

    # Postprocessing
    values = fill(zero(Qp), 0)
    spaces = fill(zero(parent(Aint)), 0)    
    
    for i = 1:length(restricted_maps)

        X = restricted_maps[i]

        if size(X,1) == 1
            push!(values, X[1,1])
            push!(spaces, invariant_blocks[i])
        else
            # Recursive call
            E = _eigenspaces_by_power_iteration(X)

            # Merge
            for j = 1:length(E.spaces)
                push!(values, E.values[j])
                push!(spaces, invariant_blocks[i]*E.spaces[j])
            end
        end
        
    end
    
    return EigenSpaceDec(Qp, values, spaces)
end


############################################################################################

"""
    block_data(A)
(Non-critical testing function). Print the valuations of the main/sub diagonal.
"""
function block_data(A)
    n = size(A,2)
    data = Array{Any,2}(nothing, 2, n)
    
    for i=1:n-1
        data[:,i] = valuation.([A[i,i], A[i+1,i]])
    end

    data[:,n] = [valuation(A[n,n]), nothing]
    
    return data
end

@doc Markdown.doc"""
    diagonal_block_ranges(A)

Given a square matrix `A`, return the ranges `1:k` so that `A[1:k, 1:k]` is a diagonal block of `A`
"""
function diagonal_block_ranges(A; pad=0)

    check_square(A)
    n = size(A,1)
    n == 0 && return Vector{UnitRange{Int}}()

    # Search from the bottom-right corner, and update the row of highest index
    # known to contain an entry from the block.
    mark = (1,1)
    i = n
    j = 1
    while j <= mark[1]
        while i > mark[1]
            if !iszero(A[i, j])
                mark = (i, j)
                i = n
                break
            end
            i -= 1
        end
        j += 1
    end

    block_size = mark[1]

    # Find the blocks 
    blocks = diagonal_block_ranges(view(A, block_size+1:n, block_size+1:n), pad=pad+block_size)
    insert!(blocks, 1, pad+1 : pad+block_size)    
    return blocks
end


function isweak_block_schur_hessenberg(S)
    # Returns true if and only if S is a (weak) block Schur form. It is
    # assumed that S is also a Hessenberg matrix.

    n = size(S,2)
    Qp = base_ring(S)
    f = charpoly(modp.(S))
    # rts_and_muls = roots_with_multiplicities(f)

    # 1. Determine the blocks
    blocks = []
    block_start = 1
    
    for j=1:n
        if j==n || iszero(S[j+1,j])
            push!(blocks, S[block_start:j, block_start:j])
            block_start = j+1
        end
    end

    block_charpolys = [charpoly(parent(f), modp.(B)) for B in blocks]

    # 2. Ensure that if g is a characteristic polynomial of a block, then either g has no roots
    #    or g = (t-lambda)^m.

    for g in block_charpolys
        g_rt_muls = roots_with_multiplicities(g)

        if length(g_rt_muls) == 0
            continue
        elseif length(g_rt_muls) > 1
            return false
        elseif g_rt_muls[1][2] != degree(g)
            return false
        end
    end

    # 3. Check that the product of the block charpolys is the original charpoly.
    #
    # Note with the computation above this implies that all roots of f mod p have a set
    # of associated blocks.

    return prod(block_charpolys) == f ? true : false
end
