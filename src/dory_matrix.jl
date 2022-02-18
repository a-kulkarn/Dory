
## Make some generic extensions to the matrix utilities in Nemo/Hecke


# I don't actually want LinearAlgebra as a dependency, but I do want to mimick the syntax
# to provide familiarity for the user.
#import LinearAlgebra: eigen, Eigen, eigvals, eigvecs

import Base: getindex, setindex, similar, /
import Hecke: matrix, identity_matrix, diagonal_matrix
import Hecke.AbstractAlgebra: check_square, issquare

##############################################################################################
#                                                                                            #
#                             Basic interface                                                #
#                                                                                            #
##############################################################################################

######################################
# Broadcasting
######################################


struct MyStyle <: Base.BroadcastStyle end

bcstyle = Broadcast.Broadcasted{MyStyle,
                                 Tuple{Base.OneTo{Int64},Base.OneTo{Int64}},
                                 <:Function,
                                 <:Tuple{<:Hecke.Generic.MatElem}}

    

Base.BroadcastStyle(::Type{<:Hecke.Generic.MatElem{T}} where T) = MyStyle()
Base.broadcastable(A::Hecke.Generic.MatElem{T} where T) = deepcopy(A)


function Base.similar(bc::bcstyle, output_eltype::Type{<:NCRingElem}, tup::Any)

    # Scan the inputs for a nemo matrix:
    A = find_nemo_mat(bc)

    # Determine finer information about the output of the broadcasted function.
    R = base_ring(A)
    y = bc.f(zero(R))
    return matrix(parent(y), fill(y, axes(A)))
end

# In the case the output type is not an NCRingElem, return a Julia array.
function Base.similar(bc::bcstyle, output_eltype::Type{T}, tup::Any) where T

    # Scan the inputs for a nemo matrix:
    A = find_nemo_mat(bc)
    return Matrix{output_eltype}(undef, size(A,1), size(A,2))
end

function Base.copyto!(X::Hecke.Generic.MatElem, bc::bcstyle)
    Y = bc.args[1]
    for i = 1:size(X,1)
        for j = 1:size(X,2)
            X[i,j] = bc.f(Y[i,j])
        end
    end    
    return X
end

# function Base.copyto!(X::Array{T,2} where T, bc::bcstyle)
#     Y = bc.args[1]
#     X = bc.f.(Y.entries)
#     return X
# end

# # We need to do something else for nmod_mats again.
# function Base.copyto!(X::Hecke.nmod_mat, bc::bcstyle)
#     Y = bc.args[1]
#     X = matrix(X.base_ring, bc.f.(Y.entries))
#     return X
# end

# # ... and gfp_mat
# function Base.copyto!(X::Hecke.gfp_mat, bc::bcstyle)
#     Y = bc.args[1]
#     X = matrix(X.base_ring, bc.f.(Y.entries))
#     return X
# end



find_nemo_mat(bc::Base.Broadcast.Broadcasted) = find_nemo_mat(bc.args)
find_nemo_mat(args::Tuple) = find_nemo_mat(find_nemo_mat(args[1]), Base.tail(args))
find_nemo_mat(x) = x
find_nemo_mat(a::Hecke.Generic.MatElem, rest) = a
find_nemo_mat(::Any, rest) = find_nemo_mat(rest)


###################################################################
#
#  Indexing
#
###################################################################

#=
function Base.getindex(A::Hecke.Generic.MatSpaceElem{T} where T, koln::Colon, I::Array{Int64,1})
    return matrix(A.base_ring, A.entries[koln,I])
end

function Base.getindex(A::Hecke.Generic.MatSpaceElem{T} where T, I::Array{Int64,1}, koln::Colon)
    return matrix(A.base_ring, A.entries[I,koln])
end

function Base.getindex(A::Hecke.Generic.MatSpaceElem{T} where T, I::Array{Int64,1}, J::Array{Int64,1})
    return matrix(A.base_ring, A.entries[I,J])
end

function Base.getindex(A::Hecke.Generic.MatSpaceElem{T} where T, I::CartesianIndex{2})
    return A[I[1],I[2]]
end

function Base.setindex!(A::Hecke.Generic.MatSpaceElem{T} where T, x, I::CartesianIndex{2})
    return setindex!(A,x,I[1],I[2])
end


function Base.collect(A::Hecke.Generic.MatSpaceElem{T}, state=1) where T
    return A.entries
end

function Hecke.matrix(A::Array{T,2} where T <: Hecke.NCRingElem)
    @assert reduce(==, [parent(x) for x in A]) 
    return matrix(parent(A[1,1], A))
end

function Hecke.matrix(A::Array{Array{T,1},1} where T <: Hecke.NCRingElem)
    return matrix(hcat(A...))
end

=#

function /(A :: Hecke.Generic.Mat{T}, x::T)  where T
    return deepcopy(A) * inv(x)
end


# Typesafe version of hcat-splat. Apparently there is a way to make this more efficient.
# function colcat(L::Array{T,1} where T <: Hecke.Generic.Mat{S} where S)
#     if isempty(L)
#         return T
#     end
# end

function check_square(A::AbstractMatrix)
    issquare(A) || throw(DomainError(A, "matrix must be square"))
    return A
end

function issquare(A::AbstractMatrix)
    return size(A,1) == size(A,2)
end

# import Hecke.*
# function *(a::fmpz, A::Hecke.MatElem)
#     return base_ring(A)(a) * A
# end


##############################################################################################
#                                                                                            
#    diagonal matrix fix
#                                                                                            
##############################################################################################

# An actually sensible version of identity_matrix
function unaliased_identity_matrix(R, n::Int)
    D = Hecke.MatrixSpace(R, n, n)()
    for i = 1:n
        D[i,i] = deepcopy(one(R))
    end
    return D
end

function _has_any_shared_refs(B)
    n = size(B,1)
    m = size(B,2)
    boo = any(B[i,j] === B[k,l] for i=1:n, j=1:m, k=1:n, l=1:m if (i != k || j != l))
    return boo
end

##############################################################################################
#                                                                                            #
#                          Correction to nullspace                                           #
#                                                                                            #
##############################################################################################


@doc Markdown.doc"""
    my_nullspace(A :: T) where T <: Union{nmod_mat, fmpz_mat, gfp_mat}

Computes the nullspace of a matrix `A` with the indicated types. Fixes the bug in Flint with the
correct return for the nullspace of the zero matrix.

(Should just do a pull-request to Nemo to fix this.)
"""
function my_nullspace(A::Hecke.Generic.MatElem)
    if iszero(A)
        return size(A,2), identity_matrix(base_ring(A), size(A,2))
    end
    return nullspace(A)
end


##############################################################################################
#                                                                                            #
#                          Generic Eigenvalue/Eigenvector functions                          #
#                                                                                            #
##############################################################################################


struct MyEigen{T}
    base_ring::Hecke.Generic.Ring
    values::Array{T,1}
    vectors::Hecke.Generic.MatElem{T}
end

struct EigenSpaceDec{T}
    base_ring::Hecke.Generic.Ring
    values::Array{T,1}
    spaces::Array{S, 1} where S <: Hecke.Generic.MatElem{T}
end

@doc Markdown.doc"""
    eigspaces(A::Hecke.Generic.MatElem{T}) where T -> EigenSpaceDec{T}

Computes the eigenspaces of a generic matrix, and returns a list of
matrices whose columns are generators for the eigen spaces.
"""
function eigspaces(A::Hecke.Generic.MatElem{T}) where T
    R,_ = PolynomialRing(A.base_ring)
    g = charpoly(R, A)
    rts = roots(g)
    if isempty(rts)
        rts = Array{T,1}()
    end
    
    Imat = identity_matrix(A.base_ring, size(A,1))

    return EigenSpaceDec( A.base_ring, rts, [ my_nullspace(A-r*Imat)[2] for r in rts])
end
    
# Returns an eigen factorization structure like the default LinearAlgebra.eigen function.
#
"""
    eigen(A::nmod_mat)

Computes the Eigenvalue decomposition of `A`. Requires factorization of polynomials implemented
over the base ring.

(Depreciated. `eigspaces` is better to use.)
"""
function eigen(A::Hecke.Generic.MatElem{T}) where T
    E = eigspaces(A)
    eig_vals = Array{T,1}(vcat([fill( E.values[i] , size(E.spaces[i],2) ) for i=1:size(E.values,1)]...))
    eig_vecs = _spacecat(E)
    return MyEigen(E.base_ring, eig_vals, eig_vecs)
end

# Typesafe version of hcat-splat
function _spacecat(E::EigenSpaceDec)
    if isempty(E.spaces)
        return matrix(E.base_ring, fill(zero(FlintZZ),0,0))
    else
        return hcat(E.spaces...)
    end
end


@doc Markdown.doc"""
    eigvecs( A :: Hecke.Generic.MatElem{T}) where T -> A :: Hecke.Generic.MatElem{T}

Return a matrix `M` whose columns are the eigenvectors of `A`. (The kth eigenvector can be obtained from the
slice `M[:, k]`.)
"""
function eigvecs(A::Hecke.Generic.MatElem{T}) where T
    return _spacecat(eigspaces(A))
end


@doc Markdown.doc"""
    eigvals(A::Hecke.Generic.MatElem{T}) where T -> values :: Array{T,1}
Return the eigenvalues of `A`.
"""
function eigvals(A::Hecke.Generic.MatElem{T}) where T
    return eigen(A).values
end


# Needs to be more robust. Also applied to the situation A is square but not of rank 1.
@doc Markdown.doc"""
    rectangular_solve(A::Hecke.MatElem{T}, b::Hecke.MatElem{T}; suppress_error=false) where T 
                                                                                    --> x ::Hecke.MatElem{T}

Solve the possibly overdetermined linear equation Ax = b. If no solution exists returns an error.
Generally not intended for use.

WARNINGS:
---------
This function is very unsafe. It does not do basic sanity checks and will fail if the top
nxn block is singular.
"""
function rectangular_solve(A::Hecke.MatElem{T}, b::Hecke.MatElem{T}; suppress_error=false) where T

    if !suppress_error
        error("Not implemented correctly. If you really want to use this function, set 'suppress_error=true'.")
    end
    
    rows(A) < cols(A) && error("Not implemented when rows(A) < cols(A)")

    # Extract top nxn block and solve using standard method.
    B = A[1:cols(M),:]
    x = solve(B,b[1:cols(M),:])

    !iszero(A*x - b) && error("Linear system does not have a solution")
    return x
end

