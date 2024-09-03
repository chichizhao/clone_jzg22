# Build the phylogenetic tree and map the clonal development routes
## 1. Build the phylogenetic tree
####  here we use the pythoon scipt to achieve the goal
#### prepare the input file
    solid somatic mutation sites: typha_somatic_mutation.vcf
    
    python evolution_tree_structure.py -i typha_somatic_mutation.vcf -o typha_somatic_mutation_tree_structure.txt

#### plot the phylogenetic tree using the tree structure file
    python plot_tree.py
## 2. Map the clonal development routes
    matlab map_clonal_development_routes.m