# Count the characteristics of somatic mutations character in each sample
## 1.count the position of somatic mutations in each sample
#### here we use python to count the position of somatic mutations in each sample
    python check_pos_somatic_mutations.py
    python check_pos_somatic_mutations_all.py

## 2.count the direction of somatic mutations in each sample
### here we use python to count the direction of somatic mutations in each sample and all samples
    python check_direction_somatic_mutations.py
    
## 3. check the non-synonymous mutation sites

#### here we use the python script to achieve the function.
#### here are some preciple for the non-synonymous mutation sites check
    1. the mutation sites should be in the coding region
    2. the mutation sites should be non-synonymous mutation
    3. the mutation sites should be check in time order and one sample by one sample
#### here we use the python script to achieve the function
#### first we need to annotate the mutation sites with the gene information time order and one sample by one sample
    python annotate_mutation_type.py -vcf solid_somatic_mutation.vcf -gff annotation.gff3 -o mutation_type.txt
    python non_synonymous_check.py -vcf solid_somatic_mutation.vcf -tree typha_somatic_mutation_tree.txt -gff annotation.gff3 -mut_type mutation_type.txt -ref typha.fa -out typha_somatic
    
## 4. check the VAF of somatic mutations
#### here we use the python script to check the VAF of somatic mutations
    python check_vaf_germline_calls.py
    python check_vaf_somatic_mutations.py
