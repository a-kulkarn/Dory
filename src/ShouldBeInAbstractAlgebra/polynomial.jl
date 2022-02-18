
###############################################################################
#
#   Monomials of degree
#
###############################################################################

export monomials_of_degree, dense_coefficients

@doc Markdown.doc"""
    monomials_of_degree(X::Vector{<:Hecke.MPolyElem}, itr)

Given a list of variables X from a polynomial ring, return all monomials of degrees 
specified by `itr`, which can be either a nonnegative integer or an iterable collection
of nonnegative integers.
"""
function monomials_of_degree(X::Vector{T}, itr)::Vector{T} where T<:Hecke.MPolyElem
    return [mon for d in itr for mon in monomials_of_degree(X,d)]
end

function monomials_of_degree(X::Vector{T}, d::Int64)::Vector{T} where T<:Hecke.MPolyElem

    # Edge case returns
    isempty(X) && error("Variables required in input")
    length(X) == 1 && return [X[1]^d]
    d < 0 && return Array{S,1}()
    d ==0 && return [parent(X[1])(1)]
    
    all_monomials = typeof(X)()
    n = length(X)
    Y = X[1:n-1]
    for j=0:d
        for m in monomials_of_degree(Y, d-j)
            push!(all_monomials, X[n]^j * m)
        end
    end
    return all_monomials
end

@doc Markdown.doc"""
    monomials_of_degree(R::Hecke.MPolyRing{T} where T, itr)

Given a polynomial ring, return all monomials of degrees specified by `itr`, 
which can be either a nonnegative integer or an iterable collection of nonnegative integers.
"""
function monomials_of_degree(R::Hecke.MPolyRing{T} where T, ds)
    X = gens(R)
    return monomials_of_degree(X, ds)
end


###############################################################################
#
#   dense_coefficients (Recommended for user use only)
#
###############################################################################


@doc Markdown.doc"""
    dense_coefficients(f::Hecke.MPolyElem{T}) where T <: RingElement

Return the list of monomial coefficients of the polynomial `f` for every monomial of degree up to
the total degree of `f`.
in the list L.
"""
function dense_coefficients(f::Hecke.MPolyElem{T}) where T <: RingElement
    R = parent(f)
    mons = monomials_of_degree(gens(R), 0:total_degree(f))
    return [coeff(f, m) for m in mons]
end

@doc Markdown.doc"""
    coeff(f::Hecke.MPolyElem{T}, L) where T <: RingElement

Return the list of monomial coefficients of the polynomial `f` specified by the monomials
in the list L.
"""
function coeff(f::Hecke.MPolyElem{T}, L) where T <: RingElement
    return [Hecke.coeff(f,m) for m in L]
end


###############################################################################
#
#   Roots and factoring
#
###############################################################################

@doc Markdown.doc"""
    roots_with_multiplicities(f)

Returns a list of roots of the polynomial `f`, with the multiplicities occuring 
in the factorization.
"""
function roots_with_multiplicities(f)
    F = Hecke.factor(f)
    return [(-g(0), m) for (g,m) in F if Hecke.degree(g) == 1]
end
