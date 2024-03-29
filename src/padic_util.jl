##############################################################################################
#                                                                                            
#                          Input types for linear algebra
#                                                                                            
##############################################################################################

DiscreteValuedField = Union{NonArchLocalField, Hecke.Generic.LaurentSeriesField}
DiscreteValuedFieldElem = Union{NonArchLocalFieldElem, Hecke.Generic.LaurentSeriesFieldElem}

##############################################################################################
#                                                                                            
#                          Basic extension to padic numbers                                  
#                                                                                            
##############################################################################################

import Base: +, /, abs
import Hecke: inv, lift, nullspace, valuation


/(x::T, y::T) where T <: DiscreteValuedFieldElem = x // y

function _unsafe_minus!(x::padic, y::padic)
    x.N = min(x.N, y.N)
    ccall((:padic_sub, Hecke.:libflint), Nothing,
          (Ref{padic}, Ref{padic}, Ref{padic}, Ref{FlintPadicField}),
          x, x, y, parent(x))
    return x
end


function isapprox(x::DiscreteValuedFieldElem, y::DiscreteValuedFieldElem;
                  valuation_atol::Real = min(precision(x), precision(y)),
                  atol::Real=0,
                  norm::Function=abs)

    # TODO: Implement the relative tolerance functionality, similar to Julia.
    x == y && return true
    z = x - y
    
    atol == 0 && return valuation(z) >= min(valuation_atol, precision(z))
    return abs(z) <= atol
end

function isapprox_zero(x::DiscreteValuedFieldElem;
                       valuation_atol::Real = precision(x),
                       atol::Real=0,
                       norm::Function=abs)

    # TODO: Implement the relative tolerance functionality, similar to Julia.    
    iszero(x) && return true
    
    atol == 0 && return valuation(x) >= min(valuation_atol, precision(x))
    return abs(x) <= atol
end


# typesafe version
function float64_valuation(x::DiscreteValuedFieldElem)
    iszero(x) && return Inf
    return Float64(valuation(x))
end


function abs(x::NonArchLocalFieldElem)
    p = Hecke.prime(parent(x))
    return Float64(p)^(-valuation(x))
end

function abs(x::Hecke.Generic.LaurentSeriesFieldElem)
    R = base_ring(parent(x))

    e = characteristic(R) == 0 ? exp(1) : characteristic(R)
    return Float64(e)^(-valuation(x))
end

function modp(x::NonArchLocalFieldElem)
    Qp = parent(x)
    Fp = ResidueRing(FlintZZ, Qp.p)
    return Fp(lift(x))
end

function modp(x::Hecke.Generic.LaurentSeriesFieldElem)

    if valuation(x) < 0
        throw(DomainError(x, "Cannot reduce element of negative valuation to residue field."))
    end
    return Hecke.coeff(x, 0)
end

# Function to extract information from a local field element.
elt_info(x) = (iszero(x), valuation(x), precision(x))


##############################################################################################
#                                                                                            #
#                          Randomization                                                     #
#                                                                                            #
##############################################################################################


function randint(K::Hecke.Generic.LaurentSeriesField)
    pi = uniformizer(K)
    N = precision(K)
    F = base_ring(K)
    
    return sum(pi^i * F(rand(Int)) for i=0:N-1)
end

function randint(Qp::Hecke.FlintPadicField)
    p = prime(Qp)
    N = precision(Qp)
    return Qp(rand(1:BigInt(p)^N))
end

function randunit(K::DiscreteValuedField)
    for i=1:100
        x = randint(K)
        valuation(x) == 0 && return x
    end
    error("Cannot generate random unit in the ring of integers.") 
end

function random_test_matrix(Qp, n=4, m=n)
    A = matrix(Qp, fill(zero(Qp),n,m))
    for i=1:n
        for j=1:m
            A[i,j] = randint(Qp)
        end
    end
    return A
end

function random_rotation_matrix(K, n=5)
    A = matrix(K, fill(zero(K), n, n))

    for i=1:100
        A = random_test_matrix(K, n)
        valuation(det(A)) == 0 && return A
    end
    error("Cannot generate a random matrix in GL_n(OK).")
end

function random_matrix_with_eigenblock(Qp, n=4, block_size=2; jordan=false)

    # Specify a diagonal matrix with a large block.
    a = randint(Qp)

    rand_eigvals = Array{elem_type(Qp)}([randint(Qp) for i=1:(n-block_size)])
    diag = vcat(rand_eigvals, fill(a, block_size))

    D = diagonal_matrix(diag)
    if jordan
        for i = n-block_size+1 : n-1
            D[i,i+1] = 1
        end
    end
    
    # Choose a random GL_n(Zp) matrix.
    B = random_rotation_matrix(Qp, n)

    return inv(B) * D * B, diag, B
end

##############################################################################################
#                                                                                            #
#                          Polynomials over p-adic fields                                    #
#                                                                                            #
##############################################################################################

# Getting coefficients of a flint polynomial is not intuitive.
# Flint system crashes if coefficients are not integers.
# Flint system crashes if leading coefficients are divisible by p.

# Lift termwise to a polynomial over the Flintegers.
function lift(f :: Hecke.Generic.Poly{<:DiscreteValuedFieldElem})
    R,_ = PolynomialRing(FlintZZ)
    return R([lift(c) for c in f.coeffs])
end

# This function is...kind of a hack.
# It is also very buggy since FLINT can only handle a specific case
# (integer polynomial, non-vanishing leading coefficient mod p)
function factor(f :: Hecke.Generic.Poly{<:DiscreteValuedFieldElem})
    QpX = f.parent
    Qp = QpX.base_ring
    N = prec(Qp)
    
    f_int = lift(f)
    H = factor_mod_pk_init(f_int,Qp.p)
    D = factor_mod_pk(H,N)

    return Dict(change_base_ring(Qp,k)=>D[k] for k in keys(D))   
end

