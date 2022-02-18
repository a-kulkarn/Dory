
using Dory, Test

@testset "Broadcasting" begin

    PR, x = PolynomialRing(QQ)
    rings = [ZZ, QQ, PadicField(7,6), PR]

    fivemat = matrix(ZZ, fill(5, 5, 5))
    
    for R in rings
        A = zero_matrix(R, 5, 5)
        B = zero_matrix(R, 5, 0)

        @test A == zero.(A)
        @test fivemat == (x->fmpz(5)).(A)

        # If the output function is a Julia type, return an array.
        @test fill(5, 5, 5) == (x->5).(A)
        @test fill(true, 5, 5) == (x->true).(A)

        # Test broadcast on empty matrices
        @test base_ring((x->x+1).(B)) == R
        @test size((x->x+1).(B)) == (5, 0)
        @test size((x->1).(B)) == (5, 0)

        # Test erroring
        will_error = x->GF(11)(x)

        if !(R in [ZZ, QQ])
            @test_throws Exception will_error.(A)
            @test_throws Exception will_error.(B)
        end
    end
end



@testset "QR method tests" begin

    Qp = PadicField(7,20)
    LF1, _ = LaurentSeriesField(GF(5), 20, "t")
    LF2, _ = LaurentSeriesField(QQ, 3, "t")
    
    for K in [Qp, LF1, LF2]
        @testset "Basic version $elem_type(K)" begin
            for n=1:10
                A = zero_matrix(K,n,n)
                F = padic_qr(A)
                @test A[F.p,F.q] == F.Q*F.R
                @test isupper_triangular(F.R)

                for i=1:5
                    B = random_test_matrix(K, n)
                    F = padic_qr(B)
                    @test B[F.p,F.q] == F.Q*F.R
                    @test isupper_triangular(F.R)
                end
            end
        end

        @testset "Column pivoting $elem_type(K)" begin
            for n=1:10
                A = zero_matrix(K,n,n)
                F = padic_qr(A, col_pivot=Val(true))
                @test A[F.p,F.q] == F.Q*F.R
                @test isupper_triangular(F.R)

                for i=1:5
                    B = random_test_matrix(K, n)
                    F = padic_qr(B, col_pivot=Val(true))
                    @test B[F.p,F.q] == F.Q*F.R
                    @test isupper_triangular(F.R)
                end
            end
        end
    end
end


@testset "Singular value decomposition" begin

    N = 30
    Qp = PadicField(7, N)
    LF1, _ = LaurentSeriesField(GF(5), N, "t")
    LF2, _ = LaurentSeriesField(QQ, 3, "t")

    
    for K in [Qp, LF1, LF2]
        for n = 1:10
            for i = 1:5
                A = random_rotation_matrix(K, n)
                @test all(x->(valuation(x) == 0), singular_values(A))
            end
        end

        for n = 1:10
            for i = 1:2
                sing_vals = [uniformizer(K)^rand(1:precision(K)-1) for i=1:n]
                D = diagonal_matrix(sing_vals)
                U = random_rotation_matrix(K, n)
                V = random_rotation_matrix(K, n)

                B = U * D * V
                F = svd(B)
                @test F.U * F.S * F.Vt == B[F.p, F.q]
                @test isdiagonal(F.S)
                
                comped_singvals = diagonal(F.S)
                @test valuation.(comped_singvals) == sort(valuation.(sing_vals))

                @test minimum(precision.(F.U)) == minimum(precision.(F.Vt)) == precision(K)

                minvals1 = minimum(valuation.(singular_values(F.U)))
                minvals2 = minimum(valuation.(singular_values(F.Vt)))
                @test minvals1 == minvals2 == 0

                maxvals1 = maximum(valuation.(singular_values(F.U)))
                maxvals2 = maximum(valuation.(singular_values(F.Vt)))
                @test maxvals1 == maxvals2 == 0
            end
        end

    end
end


# -- Eigenvectors

# -- More block_schur_form.


@testset "Block Schur form tests" begin
    
    N = 30
    Qp = PadicField(7, N)
    LF1, _ = LaurentSeriesField(GF(5), N, "t")
    LF2, _ = LaurentSeriesField(QQ, 3, "t")

    for K in [Qp, LF1, LF2]
        for n=1:10
            for i=1:5
                A = random_test_matrix(Qp, n)
                S, V = Dory.block_schur_form(A)

                @test valuation(det(V)) == 0
                @test inv(V)*S*V == A
                @test Dory.isweak_block_schur_hessenberg(S)
                @test minimum(precision.(V)) == precision(Qp)
            end
        end
    end
end
