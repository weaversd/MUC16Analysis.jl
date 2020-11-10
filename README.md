# MUC16Analysis - UPGMA of amino acid sequence for phlyogenetic analysis


To import MUC16Analysis julia package

```
using Pkg
Pkg.add(PackageSpec(url="https://github.com/weaversd/MUC16Analysis.jl.git"))
using MUC16Analysis
```

To test the package

```
Pkg.test("MUC16Analysis")
```

**<ins>Main Workflow:</ins>** 
fasta file of aa sequences to compare -> distance matrix -> format distance matrix -> upgma -> newick file

It is possible to run UPGMA clustering on any distance matrix (doesn't have to be from protein sequences). To do this, the distance matrix must be in a csv file where the first row is column names (these are the leaf names in the tree). The rows don't have labels, but the order must be the same as the columns (which will result in 0's in the diaganol). See below for an example:


If this is the case you can run `create_dm_df()` with the csv file, followed by `upgma()` which will perform the clustering. Alternatively, you can create a distance matrix in the same format described above as a DataFrame object in julia, and run `add_missing()` on it, followed by `upgma()`.

## Functions for Main Workflow  
*Note that arguments prior to ';' do not have keywords and are positional. Arguments after ';' need the keyword specified. This is most important in `upgma()`. 
(e.g. `upgma(dist_mat_df, ["A", "C", "J"], output = "tree_file.tree", header = false, verbose = false, distances = true`)*

**dist_mat_gen()**
*to generate a distance matrix based on fasta. This returns the distance matrix as a dataframe, and saves it as a csv file*  
*uses a substitution matrix from Ali et al., 2016 that is based on chemical differences between residues*
```
dist_mat_gen(fasta_file, output_file; insert, delete, describer)
```
* fasta_file (string) is file path of fasta of AA sequences  
* output_file (string) is the location where the csv of the distance matrix will be saved  
* insert (integer) is the penalty applied for an insertion. Default is 150 (for reference: substitutions are 0-100)  
* delete (integer) is the penalty applied for a deletion. Default is 150 (for reference: substitutions are 0-100)  
* desciber (string) is what to name the sequences (ends up as column names in the dist. mat.) default is "identifier". Other option is "description"  


**add_missing()**
*to format the distance matrix so it can be used in the upgma analysis. This makes it a "lower left" matrix where the diaganol and everything above and to the right of the diaganol is 'missing'. Essentially eliminates duplicate values.*
```
add_missing(dataframe)
```
* dataframe (data frame object) is the data frame that will be  modified.


**create_dm_df()**
*to take a csv file of a distance matrix (first row is column headers) and turn it into a dataframe that can be used in upgma analysis. automatically runs `add_missing()` so the returned dataframe can be run directely with `upgma`.*
```
create_dm_df(csv_file)
```
* csv_file is the file path of csv of distance matrix to turn into a dataframe


**upgma()**
*runs unweighted pair group method with arithmetic mean (upgma) on the distance matrix generated with `create_dm_df()` or `dist_mat_gen()` + `add_missing()`. Returns newick file as a string, and saves newick as text file*
```
upgma(distance_matrix, labels(optional); output, header, verbose, distances)
```
* distance_matrix (data frame) is a distance matrix generated by one of the aformentioned methods.
* labels (array of strings) is an optional variable that is an array of labels to be used as identifiers for the samples in the newick file. It can be omitted, in which case *header* must be *true*, and identifiers will come from the column names of the dataframe.
* output (string) is the file path of the output newick file (a text file in newick format). Default is "file.newick"
* header (boolean) is whether the names of the samples will be taken from the column names of the dataframe. Default is true, and should work if the dataframe is generated with the aformentioned methods. If true, will override any array of strings provided in *labels*.
* verbose (boolean) is whether to print all the dataframes during the production of the newick file. This can be a lot of dataframes, so only set to true for troubleshooting and if the labels are short (one or two characters) and the data frame is small (< 10 rows/colums). Default is false
* distances (boolean) is whether to include the distance calculations in the newick file. Default is true. If false, the file will only have sample labels.


## Other Included Functions

**fasta_to_text()**
*converts a one record fasta file to a text file to use with `subset_sequence()`*
```
fasta_to_text(fasta_file, output_file)
```
* fasta_file (string) is the filepath to the fasta file to be converted
* output_file (string) is the filepath where the output text file is saved

**subset_sequence()**
*takes a text amino acid sequence along with a list of amino acid positions (integers), and creates a fasta file that has each subset of the amino acid sequence in its own record based on the indeces provided*
```
subset_sequence(file, index_file, output, prot_name, alpha)
```
* file (string) the filepath to the text file with the amion acid sequence to be subsetted.
* index_file (string) is the filepath to the text file of sequence indeces. In this file, each line should have two integers separated by a comma (e.g. `34,79`). These coorespond to the first and the last amino acid positions in one of the sequences to be subsetted. The lines should be in sequential order.
* output (string) the filepath where the output fasta sequence will be saved. This will be a fasta file where each record will be one of the subsetted sequences as defined by the first two arguments
* prot_name (string) the name of the protein being subsetted (this goes in the description of each fasta record)
* alpha (Boolean) whether the identifier in the fasta file is letters, or numbers. Default is `true` which results in letters, using the `letter_list()` function. In this case, identifiers will be (A, B, C, D...X, Y, Z, AA, BB, CC..., AAA, BBB...). If false, identifiers will be ascending integers.

**letter_list()**
*produces an array of strings of letters (to use as labels in subset_sequence...or elsewhere) with specified length. In the order (A, B, C..., Y, Z, AA, BB..., ZZ, AAA, BBB...)*
```
letter_list(n)
```
* n (integer) is the length of the array

