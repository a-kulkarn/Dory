module Dory


## NOTE: It is possible to resolve the impending namespace conflict for the names of linear algebra functions
# via conditional importing (see https://github.com/mikeInnes/Requires.jl). For now, we put the burden of use on the user if they want to use both packages at the same time.

## In general, the issue is a Julia developer desicion. See https://discourse.julialang.org/t/function-name-conflict-adl-function-merging/10335/30

# Functions with the same name as LinearAlgebra functions.
#
#import LinearAlgebra.eigvecs
#import LinearAlgebra.svd
#
#

export broadcast, iterate, collect, matrix, rectangular_solve, my_nullspace, eigen, eigvecs, eigspaces, MyEigen

export /, valuation, abs, modp, test_rings, rand, rand_padic_int, random_test_matrix, padic_qr, inverse_iteration, singular_values

using Hecke, Distributed, Markdown
import Hecke: Generic.Mat, Generic.MatElem, nmod_mat
import Hecke: AbstractAlgebra

## Export the namespace of Hecke for use after Dory
#  In expectation of being rolled into Hecke, the export system is designed so that Dory mimics
# the behaviours of a submodule of Hecke.

#include("ShouldBeInHecke/HeckeTypes.jl")
#include("ShouldBeInHecke/LocalLaurent.jl")
include("ShouldBeInHecke/sparse_matrix.jl")
include("Errors.jl")
include("dory_matrix.jl")
include("sparse_dory_matrix.jl")
include("padic_util.jl")
include("polynomial_ext.jl")


function __init__()

    # Setup Verbose flags
    Hecke.add_verbose_scope(:local_QR)
    
    if myid() == 1
        println("")
        print("You've found \n")
        printstyled("DORY", color = :red)
        println()
        print("Version")
        printstyled(" $VERSION_NUMBER ", color = :green)
        print("... \n ... which comes with absolutely no warranty whatsoever")
        println()
        println("(c) 2019 by Avinash Kulkarni")
        println()
    else
        println("DORY $VERSION_NUMBER ...")
    end

end

##################
# Version Number
global const VERSION_NUMBER = "0.1"



end
