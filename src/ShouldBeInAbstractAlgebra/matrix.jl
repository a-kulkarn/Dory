##############################################################################################
#                                                                                            
#    Iterators
#                                                                                            
##############################################################################################

import Base: eachrow, eachcol
eachrow(A::Hecke.Generic.MatElem) = (view(A, i, :) for i in axes(A, 1))
eachcol(A::Hecke.Generic.MatElem) = (view(A, :, i) for i in axes(A, 2))

