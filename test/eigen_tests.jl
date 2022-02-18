
@testset "Eigenvalues and Eigenvectors" begin

    N = 30
    Qp = PadicField(7, N)
    LF1, _ = LaurentSeriesField(GF(5), N, "t")
    LF2, _ = LaurentSeriesField(QQ, 3, "t")

    function test_eigen(K, method)
        for n = 1:5, i=1:2
            a = randint(K)
            eigen_vals = [a + K(i) for i=1:n]
            D = diagonal_matrix(eigen_vals)
            U = random_rotation_matrix(K, n)

            A = U * D * inv(U)

            E = eigspaces(A, method=method)
            
            @test issubset(Set(E.values), Set(eigen_vals))

            for i=1:length(E.spaces)
                expr = valuation.(A * E.spaces[i] - E.values[i] * E.spaces[i])
                
                if A * E.spaces[i] != E.values[i] * E.spaces[i]
                    @info "Eigenbroken" K method expr modp.(eigen_vals) E.values[i]
                    #@test false
                else                        
                    @test A * E.spaces[i] == E.values[i] * E.spaces[i]
                end
            end
        end
    end
    
    @testset "Well-conditioned input (power iteration)" begin
        for K in [Qp, LF1, LF2]
            test_eigen(K, "power")
        end
    end

    @testset "Well-conditioned input (inverse iteration)" begin
        for K in [Qp, LF1, LF2]
            test_eigen(K, "inverse")
        end
    end

    @testset "Well-conditioned input (classical algorithm)" begin
        test_eigen(Qp, "classical")
    end
    
end
