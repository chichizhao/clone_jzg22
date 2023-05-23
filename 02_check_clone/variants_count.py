#!/py3
#--coding:utf-8--
# author: zhaochichi github: https://github.com/chichizhao/clone_jzg22
# function: count the varinats distribution in the each sample
# here we use the following structure to calculate the mutation site distribution
# samples
# snp_0/0
# snp_0/1
# snp_1/1
# snp_./.
# snp_total
# indel_0/0
# indel_0/1
# indel_1/1
# indel_./.
# indel_total
# mutation_total
#usage: python variants_count.py -vcf vcf_path -o output_prefix


import os
import argparse
import csv
from collections import defaultdict
import numpy as np
import gzip

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-vcf', '--vcf_path', help='each sample vcf path', required=True)
    parser.add_argument('-o', '--output', help='output prefix', required=True)
    args = parser.parse_args()


def calculate_mutation_site(vcfs_path):
    sample = []
    snp = defaultdict(list)
    snp['0/0'] = []
    snp['0/1'] = []
    snp['1/1'] = []
    snp['./.'] = []
    snp['total'] = []
    indel = defaultdict(list)
    indel['0/0'] = []
    indel['0/1'] = []
    indel['1/1'] = []
    indel['./.'] = []
    indel['total'] = []
    mutation_total = []
    i = 0
    for file in file_path:
        with gzip.open(file, 'rt') as f:
            sample.append(file.split('/')[-2].split('_')[1])
            snp['0/0'].append(0)
            snp['0/1'].append(0)
            snp['1/1'].append(0)
            snp['./.'].append(0)
            snp['total'].append(0)
            indel['0/0'].append(0)
            indel['0/1'].append(0)
            indel['1/1'].append(0)
            indel['./.'].append(0)
            indel['total'].append(0)
            for line in f:
                if line.startswith('#'):
                    continue
                else:
                    line = line.strip().split('\t')
                    if len(str(line[3])) == 1 and len(str(line[4])) == 1:
                        #print(str(line[4]))
                        snp[line[9].split(':')[0]][i] += 1
                        snp['total'][i] += 1
                    else:
                        print(str(line[4]))
                        indel[line[9].split(':')[0]][i] += 1
                        indel['total'][i] += 1
        mutation_total.append(snp['total'][i] + indel['total'][i])
        i += 1
    return sample, snp, indel, mutation_total

def write_to_csv(sample, snp, indel, mutation_total, output):
    with open(output, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['sample'] + sample)
        snp_0_0 = ['snp_0/0'] + snp['0/0']
        writer.writerow(snp_0_0)
        snp_0_1 = ['snp_0/1'] + snp['0/1']
        writer.writerow(snp_0_1)
        snp_1_1 = ['snp_1/1'] + snp['1/1']
        writer.writerow(snp_1_1)
        snp_dot_dot = ['snp_./.'] + snp['./.']
        writer.writerow(snp_dot_dot)
        snp_total = ['snp_total'] + snp['total']
        writer.writerow(snp_total)
        indel_0_0 = ['indel_0/0'] + indel['0/0']
        writer.writerow(indel_0_0)
        indel_0_1 = ['indel_0/1'] + indel['0/1']
        writer.writerow(indel_0_1)
        indel_1_1 = ['indel_1/1'] + indel['1/1']
        writer.writerow(indel_1_1)
        indel_dot_dot = ['indel_./.'] + indel['./.']
        writer.writerow(indel_dot_dot)
        indel_total = ['indel_total'] + indel['total']
        writer.writerow(indel_total)
        mutation_total = ['mutation_total'] + mutation_total
        writer.writerow(mutation_total)

def main():
    sample, snp, indel, mutation_total = calculate_mutation_site(args.vcf_file)
    write_to_csv(sample, snp, indel, mutation_total, output)

if __name__ == '__main__':
    args = get_args()
    main()

