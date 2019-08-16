
## Make some generic extensions to the matrix utilities in Nemo/Hecke


# I don't actually want LinearAlgebra as a dependency, but I do want to mimick the syntax
# to provide familiarity for the user.
#import LinearAlgebra: eigen, Eigen, eigvals, eigvecs

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
                                 S,
                                 Tuple{Hecke.Generic.MatSpaceElem{T}}} where S<:Function where T

    

Base.BroadcastStyle(::Type{<:Hecke.Generic.MatSpaceElem{T}} where T) = MyStyle()
Base.broadcastable(A::Hecke.Generic.MatSpaceElem{T} where T) = deepcopy(A)


function Base.similar(bc::bcstyle,t::Type{T} where T<:NCRingElem)
    # Scan the inputs for a nemo matrix:
    A = find_nemo_mat(bc)
    
    # Use the data fields to create the output
    if isempty(A.entries)
        return similar(A.entries)
    end
    
    val = bc.f(A[1,1])
    init = fill(val, size(A)[1], size(A)[2])

    if typeof(val) <: NCRingElem
        return matrix(val.parent, init)
    else
        return init
    end
end

# In the case the output type has no parent, or doesn't make sense as a matrix, or is a
# Julia default type (like Int64), return an Array
function Base.similar(bc::bcstyle,t::Type{T} where T)
    # Scan the inputs for a nemo matrix:
    A = find_nemo_mat(bc)
    
    # Use the data fields to create the output
    if isempty(A.entries)
        return similar(A.entries)
    end
    
    val = bc.f(A[1,1])
    init = fill(val, size(A)[1], size(A)[2])
    return init
end


function Base.copyto!(X::Hecke.Generic.MatSpaceElem{T} where T, bc::bcstyle)
    Y = bc.args[1]
    X.entries = bc.f.(Y.entries)
    X.base_ring = Y.base_ring    
    return X
end

# We need to do something else for nmod_mats again.
function Base.copyto!(X::Hecke.nmod_mat, bc::bcstyle)
    Y = bc.args[1]
    X = matrix(X.base_ring, bc.f.(Y.entries))
    return X
end


function Base.copyto!(X::Array{T,2} where T, bc::bcstyle)
    Y = bc.args[1]
    X = bc.f.(Y.entries)
    return X
end


find_nemo_mat(bc::Base.Broadcast.Broadcasted) = find_nemo_mat(bc.args)
find_nemo_mat(args::Tuple) = find_nemo_mat(find_nemo_mat(args[1]), Base.tail(args))
find_nemo_mat(x) = x
find_nemo_mat(a::Hecke.Generic.MatSpaceElem, rest) = a
find_nemo_mat(::Any, rest) = find_nemo_mat(rest)


#### end broadcast interface.


###################################################################
#  Conveinence Interface
###################################################################

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
    return matrix( hcat(A...) )
end

function Hecke.matrix(R::Hecke.Nemo.AbstractAlgebra.NCRing, A::Array{Array{T,1},1} where T)
    return matrix( R, hcat(A...) )
end


import Base./
function /(A :: Hecke.Generic.Mat{T}, x::T)  where T
    return deepcopy(A) * inv(x)
end


# Typesafe version of hcat-splat. Apparently there is a way to make this more efficient.
function colcat(L::Array{T,1} where T <: Hecke.Generic.Mat{S} where S)
    if isempty(L)
        return T
    end
end

##############################################################################################
#                                                                                            #
#                          Cosmetic override to nullspace                                    #
#                                                                                            #
##############################################################################################

# BUG:
# nullspace(x::fmpz_mat) in Nemo at /Users/avinash/.julia/packages/Nemo/XNd31/src/flint/fmpz_mat.jl:889
# Fails to compute nullspace for zero matrix.

## Make things a little more consistent with the other Julia types
@doc Markdown.doc"""
    my_nullspace(A :: T) where T <: Union{nmod_mat, fmpz_mat, gfp_mat}

Computes the nullspace of a matrix `A` with the indicated types. Fixes the bug in Flint with the
correct return for the nullspace of the zero matrix.

(Should just do a pull-request to Nemo to fix this.)
"""
function my_nullspace(A :: T) where T <: Union{nmod_mat, fmpz_mat, gfp_mat}
    if iszero(A)
        return size(A,2), identity_matrix(A.base_ring, size(A,2))
    end
    nu,N = nullspace(A)
    return nu, nu==0 ? matrix(A.base_ring, fill(0,size(A,2),0)) : N[:,1:nu]
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
    eigspaces(A::Hecke.Generic.MatElem{T}) where T --> EigenSpaceDec{T}

Computes the eigen spaces of a generic matrix, and returns a list of
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
#
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


