# MUC16Analysis


To import MUC16 julia package

```
using Pkg
Pkg.add(PackageSpec(url="https://github.com/weaversd/MUC16Analysis.jl.git"))
using MUC16Analysis
```


**main workflow:**  
fasta file of aa sequences to compare -> distance matrix -> format distance matrix -> upgma -> newick file


## Functions

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
*to format the distance matrix so it can be used in the upgma analysis. This makes it a "lower left" matrix where the diaganol and everything above and to the right of the diaganol is 'missing'. Essentially eliminates duplicate values.
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



