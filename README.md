
# Dory.jl

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://a-kulkarn.github.io/Dory/dev)

Package to extend the functionality of Nemo/Hecke. Notable additions include:

## Broadcasting for AbstractAlgebra matrices
- Broadcast `f.(A)`of a function over an AbstractAlgebra matrix `A` will now return an AbstractAlgebra matrix if the image of `f` is an AbstractAlgebra ring.

## padic linear algebra:
- padic qr-factorization.
- padic singular value decomposition.
- padically stable solving of linear systems.
- padically stable hessenburg form.
- eigenvector solver (power and inverse iterations) over local fields. 
- block schur form. [Only implemented for matrices defined over Qp]

## Requirements
- Hecke
- Nemo
- AbstractAlgebra
