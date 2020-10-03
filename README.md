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
*fasta_file (string) is file path of fasta of AA sequences  
*output_file (string) is the location where the csv of the distance matrix will be saved  
*insert (integer) is the penalty applied for an insertion. Default is 150 (for reference: substitutions are 0-100)  
*delete (integer) is the penalty applied for a deletion. Default is 150 (for reference: substitutions are 0-100)  
*desciber (string) is what to name the sequences (ends up as column names in the dist. mat.) default is "identifier". Other option is "description"  



