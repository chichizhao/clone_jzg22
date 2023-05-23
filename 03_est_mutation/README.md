# estimate the mutation rate and check the non-synonymous mutation sites
## 1.build the phylogenetic tree
#### here we base on the solid somatic call mutation build the phylogenetic tree with the principle samples come from the same mother ramet shared more same somatic mutation sites
### 1.1 build the phylogentic tree structure
#### here we use the python script to achieve the function
### prepare the input file
    solid somatic mutation sites: typha_somatic_mutation.vcf
    
    python evolution_tree_structure.py -i typha_somatic_mutation.vcf -o typha_somatic_mutation_tree_structure.txt
### 1.2 estimate the mutation rate and convert the tree structure time scale
#### here we use the python script to achieve the function
    each 216955151 > mutational_length.txt
    python evolution_tree_time.py -tree typha_somatic_mutation_tree_structure.txt -mut 5.52e-9 -ml mutational_length.txt -ref typha.fa -out typha_somatic
#### after this step, we get the tree structure with time scale, in the meantime, we get the mutation rate file marked as typha_somatic_mutation_rate.txt in a tree structure
## 2. check the non-synonymous mutation sites
#### above file provide the mutation rate for each branch, and classify the mutation sites which marked a time range of each mutation site
#### here we use the python script to achieve the function.
#### here are some preciple for the non-synonymous mutation sites check
    1. the mutation sites should be in the coding region
    2. the mutation sites should be non-synonymous mutation
    3. the mutation sites should be check in time order and one sample by one sample
#### here we use the python script to achieve the function
#### first we need to annotate the mutation sites with the gene information time order and one sample by one sample
    python annotate_mutation_type.py -vcf solid_somatic_mutation.vcf -gff annotation.gff3 -o mutation_type.txt
    python non_synonymous_check.py -vcf solid_somatic_mutation.vcf -tree typha_somatic_mutation_tree_time.txt -gff annotation.gff3 -mut_type mutation_type.txt -ref typha.fa -out typha_somatic
#### after this step,we got the non_synonymous_mutation_site_matrix.txt which record the non-synonymous mutation information for each sample.
    non_synonymous_mutation_site_matrix.txt
#### the information we plot as mutation_num_heatmap.png.
    mutation_num_heatmap.png
#### the non-synonymous mutation sites information record in non_synonymous_mutation_sample.txt.
    non_synonymous_mutation_sample.txt
## 3 scripts we used for plots
### plot the same and different mutation site that we use when build evolution tree sctruture.
    python plot_same_site.py -i1 same_site.csv -i2 different_site.csv -o same_diff_site.png
### plot the phylogenetic tree with mutation rate and non-synonymous mutation sites
    python plot_tree.py
### plot the heatmap of expanding rate 
    python plot_expanding_rate.py
### plot the pattern of mutation sites
    matlab plot_mutation_site.m
    matlab plot_non_synonymous_mutation_site.m


