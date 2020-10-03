# MUC16Analysis


To import MUC16 julia package


using Pkg

Pkg.add(PackageSpec(url="https://github.com/weaversd/MUC16Analysis.jl.git"))

using MUC16Analysis



workflow:

fasta sequence file of aa sequences to compare -> distance matrix -> format distance matrix -> upgma -> newick file

