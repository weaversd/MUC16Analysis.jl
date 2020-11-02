using MUC16Analysis
using BioSequences
using DataFrames
using Test
using CSV

@testset "MUC16Analysis.jl" begin
    @test 10 == 10 #just to make sure test is working
    
    @test letter_list(5) == ["A", "B", "C", "D", "E"] 
    
    SUBMAT_ALI = sub_mat_ali() #create the ali et al substitution matrix
    @test SUBMAT_ALI[AA_A,AA_R] === 90 #check random values in the substitution matrix
    @test SUBMAT_ALI[AA_G,AA_D] === 65
    @test SUBMAT_ALI[AA_P,AA_P] === 0
    
    example_dmdf = CSV.read(joinpath(dirname(@__FILE__), "data_files", "example.csv"))
    println(example_dmdf)
end
