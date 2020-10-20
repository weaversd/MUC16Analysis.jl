#!/usr bin/env julia

using FASTX
using BioSequences
using CSV
using DataFrames

#function to subset a text AA sequence file (txt) based on repeat domain index txt file
#file: input text file with AA sequence
#index_file: txt file with start and end AA numbers for each domain, one per line
#output: output file to save fasta file of repeats
#prot_name: text of protein name to put into the fasta file as description
#alpha is default, identifies repeats with unique letter, false identifies with number
function subset_sequence(file::String, index_file::String, output::String, prot_name::String, alpha::Bool = true)
    #import sequence from text file
    prot_txt_string = read(file, String)

    #convert to BioSequence
    prot_aa_sequence = LongSequence{AminoAcidAlphabet}(prot_txt_string)

    #import the indexes of each repeat domain
    index_df = CSV.read(index_file, header = 0)

    n_repeat = nrow(index_df)
    println("Number of Repeats:")
    println(n_repeat)

    #create matrix of empty aa BioSequences, length of the number of repeats
    repeat_matrix = fill(aa"", n_repeat)

    #put each repeat into a new position in the array
    for i in 1:n_repeat
        global repeat_matrix[i,1] = prot_aa_sequence[Int(index_df[i,1]):Int(index_df[i,2])]
    end

    #count the length of each repeat
    length_array = fill(1, n_repeat)

    for i in 1:n_repeat
        global length_array[i] = length(repeat_matrix[i,1])
    end

    #create empty string matrix
    title_matrix = fill("", n_repeat)

    #fill title matrix with repeat numbers 
    for i in 1:n_repeat
        global title_matrix[i] = string(prot_name, "_Repeat_Number", i)
    end

    #check whether identifier is number or letter, create alphabet matrix:
    if alpha == true
        alpha_matrix = letter_list(n_repeat)
    end

    #Write each aa repeat sequence to a record in a fasta file, either with alphabet or number
    if alpha == true
        open(FASTA.Writer, output) do writer
            for i in 1:n_repeat
                write(writer, FASTA.Record(string(alpha_matrix[i], " ", title_matrix[i]), string(length_array[i], " Residues"), repeat_matrix[i,1]))
            end
        end
    else
        open(FASTA.Writer, output) do writer
            for i in 1:n_repeat
                write(writer, FASTA.Record(string(i, " ", title_matrix[i]), string(length_array[i], " Residues"), repeat_matrix[i,1]))
            end
        end
    end

    println("Individual repeat sequences saved to fasta file:")
    println(output)
end

#function to convert fasta file to text file
function fasta_to_text(fasta_file::String, output_file::String)

    #read the fasta file
    r = open(FASTA.Reader, fasta_file)

    #store sequence as Biosequence variable
    for record in r
	    global seq = sequence(record)
    end


    #convert to string
    seq_str = convert(String, seq)


    open(output_file, "w") do f
	    write(f, seq_str)
    end
    println(string("saved sequence as ", output_file))
end
