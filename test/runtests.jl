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
    
    #run upgma on example dataframe
    example_dmdf = CSV.read(joinpath(dirname(@__FILE__), "example.csv"))
    add_missing(example_dmdf)
    example_newick = upgma(example_dmdf, output = "example.tree", verbose = false)
    println(example_newick)
    @test example_newick == "((((A:4.0,D:4.0):4.25,((B:0.5,F:0.5):5.75,G:6.25):2.0):6.25,C:14.5):2.5,E:17.0);"
    
    #run upgma on example distance matrix by different method
    example_dmdf2 = create_dm_df(joinpath(dirname(@__FILE__), "example.csv"))
    example_newick2 = upgma(example_dmdf2, output = "example2.tree", verbose = true)
    @test example_newick2 == "((((A:4.0,D:4.0):4.25,((B:0.5,F:0.5):5.75,G:6.25):2.0):6.25,C:14.5):2.5,E:17.0);"
    
    #create distance matrix
    example_dm_fromfasta = dist_mat_gen(joinpath(dirname(@__FILE__), "example.fasta"), "example_dm_fromfasta.csv")
    @test nrow(example_dm_fromfasta) == 10
    add_missing(example_dm_fromfasta)
    example_newick3 = upgma(example_dm_fromfasta, output = "example3.tree")
    @test example_newick3 == "((((1:2068.0,6:2068.0):138.70833333333348,(((2:1948.5,4:1948.5):155.0,10:2103.5):66.77777777777783,(3:2102.5,(5:2001.0,8:2001.0):101.5):67.77777777777783):36.43055555555566):12.229166666666515,7:2218.9375):63.17361111111131,9:2282.1111111111113);"
end
