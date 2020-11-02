using MUC16Analysis
using BioSequences
using Test

@testset "MUC16Analysis.jl" begin
    @test 10 == 10 #just to make sure test is working
    @test letter_list(5) == ["A", "B", "C", "D", "E"]
    SUBMAT_ALI = sub_mat_ali()
    @test SUBMAT_ALI[AA_A,AA_R] === 90
end
