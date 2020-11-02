using MUC16Analysis
using Test

@testset "MUC16Analysis.jl" begin
    @test 10 == 10 #just to make sure test is working
    @test letter_list(5) == ["A", "B", "C", "D", "E"]
end
