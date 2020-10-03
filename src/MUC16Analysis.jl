module MUC16Analysis

include("test_function_file.jl")
include("distance_matrix.jl")
include("letter_list.jl")
include("subset_fasta.jl")
include("upgma.jl")

export sub_mat_ali
export dist_mat_gen
export letter_list
export subset_fasta
export add_missing
export minimum_cell
export combine_strings
export combine_strings_dist
export combine_df
export upgma
export create_dm_df









#testing, not used:
function test()
  println("This works")
end

export test
export test2

end #module
