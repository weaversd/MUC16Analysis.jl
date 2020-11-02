module MUC16Analysis

include("distance_matrix.jl")
include("letter_list.jl")
include("subset_fasta.jl")
include("upgma.jl")
include("data_files/ali_submat.csv")

export sub_mat_ali
export dist_mat_gen
export letter_list
export subset_sequence
export add_missing
export minimum_cell
export combine_strings
export combine_strings_dist
export combine_df
export upgma
export create_dm_df
export fasta_to_text

end #module
