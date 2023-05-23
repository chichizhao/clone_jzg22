#!/py3
# -- coding: utf-8 --
# author: zhaochichi github: https://github.com/chichizhao/clone_jzg22
# function: get the mutationsite gene name if the mutationsite does t belong to the gene mark it as noncoding
# usage: python annotate_mutation_type.py -vcf solid_somatic_mutation.vcf -gff annotation.gff3 -o mutation_type.txt

import sys
import argparse
import re
import gzip

def get_args():
    parser = argparse.ArgumentParser(description="get the mutationsite gene name if the mutationsite does t belong to the gene mark it as noncoding")
    parser.add_argument("-vcf", "--vcf", help="vcf file", required=True)
    parser.add_argument("-gff", "--gff", help="gff file", required=True)
    parser.add_argument("-o", "--output", help="output file", required=True)
    args = parser.parse_args()
    return args

def get_mutaion_site(vcf_file):
    mutation_site = {}
    for line in vcf_file:
        if line.startswith("#"):
            continue
        line = line.strip().split("\t")
        chr = line[0]
        pos = line[1]
        mutation_site.setdefault(chr, []).append(pos)
    return mutation_site

def get_gene(gff_file):
    gene = {}
    for line in gff_file:
        line = line.strip().split("\t")
        type = line[2]
        if type == "gene":
            chr = line[0]
            start = line[3]
            end = line[4]
            id = line[8].split(";")[0].split("=")[1]
            gene.setdefault(chr, []).append([start, end, id])
    return gene

def get_gene_name(mutation_site, gene):
    mutation_site_add_gene_name = {}
    for chr in mutation_site:
        for pos in mutation_site[chr]:
            key = 0
            for i in gene[chr]:
                if int(pos) >= int(i[0]) and int(pos) <= int(i[1]):
                    mutation_site_add_gene_name.setdefault(chr, []).append([pos, i[2]])
                    key = 1
                    break
            if key == 0:
                mutation_site_add_gene_name.setdefault(chr, []).append([pos, "noncoding"])
    return mutation_site_add_gene_name

def write_file(mutation_site, output_file):
    for chr in mutation_site:
        for i in mutation_site[chr]:
            output_file.write(chr + "," + i[0] + "," + i[1] + "\n")

def main():
    vcf = args.vcf
    gff = args.gff
    output = args.output
    if vcf.endswith(".gz"):
        mutation_site = get_mutaion_site(gzip.open(vcf, "rt"))
    else:
        mutation_site = get_mutaion_site(open(vcf, "r"))
    if gff.endswith(".gz"):
        gene = get_gene(gzip.open(gff, "rt"))
    else:
        gene = get_gene(open(gff, "r"))
    mutation_site = get_gene_name(mutation_site, gene)
    write_file(mutation_site, open(output, "w"))

if __name__ == "__main__":
    args = get_args()
    main()

