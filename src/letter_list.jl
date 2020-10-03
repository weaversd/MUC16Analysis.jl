#!/usr bin/env julia

#function to generate alphabet list for naming purposes

function letter_list(letter_count::Int64)

    #set alphabet
    alpha_array = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P",
    "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"]

    alpha_array_n = length(alpha_array)
    
    #initialize array to return_array
    return_array = String[]

    #set iteration counter
    x = 1
    j = 1
    for i in 1:letter_count
        if j > 26
            j = 1
            x = x + 1
        end
        value = repeat(alpha_array[j], x)
        #println(string("i, j, x =", i, j, x))
        #println(string("value is: ", value))
        push!(return_array, value)
        j = j + 1
    end
    
    return return_array
end
