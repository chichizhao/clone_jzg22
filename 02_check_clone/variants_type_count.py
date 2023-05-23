#!/py3
# -- coding: utf-8 --
# author: zhaochichi github: https://github.com/chichizhao/clone_jzg22
# function: mutation site distribution from vcf file
# goal: to get the mutation site distribution of each gene 
# here is title of the output head csv chr, pos, ref, alt, 0/1, 1/1, 0/0, ./., number_0/1, mutation_type_num
# usage: python mutation_site_distribution.py -i input -o output

import sys
import os
import argparse
import gzip

def get_args():
    parser = argparse.ArgumentParser(description='mutation site distribution')
    parser.add_argument('-i', '--input', help='input file', required=True)
    parser.add_argument('-o', '--output', help='output file', required=True)
    args = parser.parse_args()
    return args

def cal_mutation_site_distribution(input):
    mutation_site_distribution = {}
    for line in input:
        if line.startswith('#'):
            continue
        line = line.strip()
        line = line.split('\t')
        chr = line[0]
        pos = line[1]
        ref = line[3]
        alt = line[4]
        # there are many samples in the vcf file, so we need to each mutation site of all samples
        # here we use the first sample
        mutation_site_distribution.setdefault(chr, {}).setdefault(pos, {}).setdefault(ref, {}).setdefault(alt, {}).setdefault('0/1', 0)
        mutation_site_distribution.setdefault(chr, {}).setdefault(pos, {}).setdefault(ref, {}).setdefault(alt, {}).setdefault('1/1', 0)
        mutation_site_distribution.setdefault(chr, {}).setdefault(pos, {}).setdefault(ref, {}).setdefault(alt, {}).setdefault('0/0', 0)
        mutation_site_distribution.setdefault(chr, {}).setdefault(pos, {}).setdefault(ref, {}).setdefault(alt, {}).setdefault('./.', 0)
        for sample in line[9:]:
            sample = sample.split(':')[0]
            if sample == '0/1':
                mutation_site_distribution[chr][pos][ref][alt]['0/1'] += 1
            elif sample == '1/1':
                mutation_site_distribution[chr][pos][ref][alt]['1/1'] += 1
            elif sample == '0/0':
                mutation_site_distribution[chr][pos][ref][alt]['0/0'] += 1
            else:
                mutation_site_distribution[chr][pos][ref][alt]['./.'] += 1
        # we need calculate the number of 0/1 and number of mutation type of (0/1, 1/1, 0/0, ./.), then add to the line
        mutation_site_distribution.setdefault(chr, {}).setdefault(pos, {}).setdefault(ref, {}).setdefault(alt, {}).setdefault('number_0/1', 0)
        mutation_site_distribution.setdefault(chr, {}).setdefault(pos, {}).setdefault(ref, {}).setdefault(alt, {}).setdefault('mutation_type_num', 0)
        mutation_site_distribution[chr][pos][ref][alt]['number_0/1'] = mutation_site_distribution[chr][pos][ref][alt]['0/1']
        if mutation_site_distribution[chr][pos][ref][alt]['1/1'] != 0:
            mutation_site_distribution[chr][pos][ref][alt]['mutation_type_num'] += 1
        if mutation_site_distribution[chr][pos][ref][alt]['0/0'] != 0:
            mutation_site_distribution[chr][pos][ref][alt]['mutation_type_num'] += 1
        if mutation_site_distribution[chr][pos][ref][alt]['./.'] != 0:
            mutation_site_distribution[chr][pos][ref][alt]['mutation_type_num'] += 1
        if mutation_site_distribution[chr][pos][ref][alt]['0/1'] != 0:
            mutation_site_distribution[chr][pos][ref][alt]['mutation_type_num'] += 1
    return mutation_site_distribution

def write_mutation_site_distribution(mutation_site_distribution, output):
    with open(output+'_snps.csv', 'w') as f:
        with open(output+'_indels.csv', 'w') as f1:
            f.write('chr,pos,ref,alt,0/1,1/1,0/0,./.,number_0/1,mutation_type_num' + '\n')
            f1.write('chr,pos,ref,alt,0/1,1/1,0/0,./.,number_0/1,mutation_type_num' + '\n')
            for chr in mutation_site_distribution:
                for pos in mutation_site_distribution[chr]:
                    for ref in mutation_site_distribution[chr][pos]:
                        for alt in mutation_site_distribution[chr][pos][ref]:
                            if len(alt) == 1 and len(ref) == 1:
                                f.write(chr + ',' + pos + ',' + ref + ',' + alt + ',' + str(mutation_site_distribution[chr][pos][ref][alt]['0/1']) + ',' + str(mutation_site_distribution[chr][pos][ref][alt]['1/1']) + ',' + str(mutation_site_distribution[chr][pos][ref][alt]['0/0']) + ',' + str(mutation_site_distribution[chr][pos][ref][alt]['./.']) + ',' + str(mutation_site_distribution[chr][pos][ref][alt]['number_0/1']) + ',' + str(mutation_site_distribution[chr][pos][ref][alt]['mutation_type_num']) + '\n')
                            else:
                                f1.write(chr + ',' + pos + ',' + ref + ',' + alt + ',' + str(mutation_site_distribution[chr][pos][ref][alt]['0/1']) + ',' + str(mutation_site_distribution[chr][pos][ref][alt]['1/1']) + ',' + str(mutation_site_distribution[chr][pos][ref][alt]['0/0']) + ',' + str(mutation_site_distribution[chr][pos][ref][alt]['./.']) + ',' + str(mutation_site_distribution[chr][pos][ref][alt]['number_0/1']) + ',' + str(mutation_site_distribution[chr][pos][ref][alt]['mutation_type_num']) + '\n')

def main():
    input = get_args().input
    output = get_args().output
    if input.endswith('.gz'):
        mutation_site_distribution = cal_mutation_site_distribution(gzip.open(input, 'rt'))
    else:
        mutation_site_distribution = cal_mutation_site_distribution(open(input, 'r'))
    write_mutation_site_distribution(mutation_site_distribution, output)

if __name__ == '__main__':
    main()

        


