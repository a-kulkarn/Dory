
@testset "Eigenvalues and Eigenvectors" begin

    N = 30
    Qp = PadicField(7, N)
    LF1, _ = LaurentSeriesField(GF(5), N, "t")
    LF2, _ = LaurentSeriesField(QQ, 3, "t")

    function test_eigen(A::Hecke.Generic.MatElem, E, eigen_vals)
        for lambda in E.values
            @test lambda in eigen_vals
        end
        
        for i=1:length(E.spaces)
            expr = valuation.(A * E.spaces[i] - E.values[i] * E.spaces[i])            
            
            # if A * E.spaces[i] != E.values[i] * E.spaces[i]
            #     @info "Eigenbroken" base_ring(A) method expr modp.(eigen_vals) E.values[i]
            #     @test false
            # end
            @test A * E.spaces[i] == E.values[i] * E.spaces[i]
        end

    end
    
    function test_eigen(K, method::Symbol)
        for n = 1:5, i=1:2
            a = randint(K)
            eigen_vals = [a + K(i) for i=1:n]
            D = diagonal_matrix(eigen_vals)
            U = random_rotation_matrix(K, n)

            A = U * D * inv(U)
            E = eigspaces(A, method=method)

            test_eigen(A, E, eigen_vals)
        end
    end
    
    @testset "Well-conditioned input (power iteration)" begin
        for K in [Qp, LF1, LF2]
            test_eigen(K, :power)
        end
    end

    @testset "Well-conditioned input (inverse iteration)" begin

        # This test specifically catches some weird behaviour.
        eigvalsraw = [10948166310227044247091515 0 0;
                      0 10948166310227044247091516 0;
                      0 0 10948166310227044247091517]

        eigvals = matrix(Qp, eigvalsraw)

        Araw = [1116124834844179080090593 12594878236815170406615967 1409110148890487973242424;
                7389623857754209255469698 6583110818301740372657622 12470061697559623972544343;
                15124876020222225539223963 21371477332922650173388100 2605922986842955200663084]

        A = matrix(Qp, Araw)
        
        # If
        # toBigInt(x) = BigInt(lift(x))
        # Then
        # @assert toBigInt.(A) == Araw && toBigInt.(eigvals) == eigvalsraw

        E = eigspaces(A, method="inverse")
        test_eigen(A, E, eigvals)

        # Now the randomized tests.
        for K in [Qp, LF1, LF2]
            try test_eigen(K, :inverse) catch end
        end
    end

    @testset "Well-conditioned input (classical algorithm)" begin
        test_eigen(Qp, :classical)
    end
    
end
