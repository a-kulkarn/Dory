

##############################################################################################
#                                                                                            
#    Laurent Series Functionality
#                                                                                            
##############################################################################################

import Hecke: precision, setprecision!, uniformizer, ResidueField

function precision(R::Hecke.Generic.LaurentSeriesField)
    return R.prec_max
end

function uniformizer(R::Hecke.Generic.LaurentSeriesField)
    return gen(R)
end

function ResidueField(Q::Hecke.Generic.LaurentSeriesField)
    k = base_ring(Q)
    T = elem_type(Q)
    S = elem_type(k)
    pro = function(x)
        v = valuation(x)
        v < 0 && error("elt non integral")
        z = coeff(x, 0)
        return z
    end
    lif = function(x)
        z = Q(x)
        return z
    end
    return k, MapFromFunc(pro, lif, Q, k)
end

function (R::Hecke.Generic.LaurentSeriesField)(a::Int64)
    return R(base_ring(R)(a))
end

##############################################################################################
#                                                                                            
#    Element functionality
#                                                                                            
##############################################################################################

function setprecision!(a::Hecke.Generic.LaurentSeriesFieldElem, N)
    a.prec = N
    return a
end

function _unsafe_minus!(a::Hecke.Generic.LaurentSeriesFieldElem, b::Hecke.Generic.LaurentSeriesFieldElem)
    z = add!(a, a, -b)
    return z
end

import Hecke.Generic: promote_rule
promote_rule(::Type{fmpz}, ::Type{T}) where {T<:NCRingElem} = T
