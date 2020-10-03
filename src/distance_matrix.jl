#functions to generate a distance matrix from fasta
#using packages
using CSV
using BioSequences
using BioAlignments
using DataFrames
using FASTX

#create the substitution matrix, derived from Ali et al.
function sub_mat_ali(;submat_csv::String="data_files/ali_submat.csv")
    ali_df = CSV.read(submat_csv)
    n = nrow(ali_df)

    #copy dataframe scale up and round so the values can be converted to integers
    ali_df_100 = copy(ali_df)
    for i in 1:n
        for j in 1:n
            ali_df_100[i,j] =  round(ali_df[i,j]*100, digits = 0)
        end
    end
 


    #copy PAM70 substitution matrix to new variable to use as a template for Ali matrix
    ALI_SUBMAT = copy(PAM70)

    #fill new submat based on df of Ali matrix, needs to use integers
    for i in 1:n
        for j in 1:n
            ALI_SUBMAT[i-1, j-1] = ali_df_100[i,j]
        end
    end
    return ALI_SUBMAT
end

#generate a distance matrix from an input fasta file
#uses substitution matrix from Ali et al, based on chemical differences and similiarities in amino acids
#fasta_file = file of fasta sequences to be compared
#output_file = location to store distance matrix as csv
#insert = penalty for insertion (default is 150, for reference substitutions range from 0-100)
#delete = penalty for deletion (default is 150, for reference, substitutions range from 0-100)
#describer = whether to have "identifer" from fasta file or "describer" from fasta file as the name of each sequence in the dataframe (default is identifer)
function dist_mat_gen(fasta_file::String, output_file::String="distance_matrix.csv";
    insert::Int64=150, delete::Int64=150, describer::String="identifier")
    #create empty matrices for the sequences, indentifies and descriptions
    repeat_matrix = LongSequence{AminoAcidAlphabet}[]
    identifier_matrix = String[]
    description_matrix = String[]

    #write sequence and identifier to respective arrays
    open(FASTA.Reader, fasta_file) do reader
        for record in reader
            push!(repeat_matrix, sequence(record))
            push!(description_matrix, description(record))
            push!(identifier_matrix, identifier(record))
        end
    end

    #import distance matrix:
    ALI_SUBMAT = sub_mat_ali()

    #create cost model
    costmodel = CostModel(ALI_SUBMAT, insertion = insert, deletion = delete)

    #count number of repeats
    n = length(repeat_matrix)


    #create DataFrame to hold distances
    distance_array = fill(0.0, (n, n))
    distance_matrix = convert(DataFrame, distance_array)

    #name the columns in the data frame as either the identifier or the descritpion from fasta
    if describer == "identifier"
        rename!(distance_matrix, identifier_matrix)
    else   
        rename!(distance_matrix, description_matrix)
    end

    for i in 1:n
        for j in 1:n
            alignment = distance(pairalign(EditDistance(), repeat_matrix[i], repeat_matrix[j], costmodel))
            distance_matrix[i,j] = alignment
        end
    end

    #write distance matrix to csv
    CSV.write(output_file, distance_matrix)
    println(string("output file is ", output_file))

    #return distance matrix
    return distance_matrix
end
