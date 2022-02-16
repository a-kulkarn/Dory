

##############################################################################################
#                                                                                            
#    Laurent Series Functionality
#                                                                                            
##############################################################################################

import Hecke: precision, setprecision!, uniformizer

function precision(R::Hecke.Generic.LaurentSeriesField)
    return R.prec_max
end

function uniformizer(R::Hecke.Generic.LaurentSeriesField)
    return gen(R)
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
