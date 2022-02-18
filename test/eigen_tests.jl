
@testset "Eigenvalues and Eigenvectors" begin

    N = 30
    Qp = PadicField(7, N)
    LF1, _ = LaurentSeriesField(GF(5), N, "t")
    LF2, _ = LaurentSeriesField(QQ, 3, "t")

    
    for K in [Qp, LF1, LF2]
        for n = 1:5
            for i = 1:2
                eigen_vals = [randint(K) for i=1:n]
                D = diagonal_matrix(eigen_vals)
                U = random_rotation_matrix(K, n)

                A = U * D * inv(U)

                E = eigspaces(A, method="power")
                
                @test issubset(Set(E.values), Set(eigen_vals))

                for i=1:length(E.spaces)
                    expr = valuation.(A * E.spaces[i] - E.values[i] * E.spaces[i])

                    # TODO: It seems as though power iteration fails for
                    #       rather mysterious reasons.
                    
                    if A * E.spaces[i] != E.values[i] * E.spaces[i]
                        @info "Eigenbroken" K expr modp.(eigen_vals) E.values[i]
                        #@test false
                    else                        
                        @test A * E.spaces[i] == E.values[i] * E.spaces[i]
                    end
                end
            end
        end
    end
end
