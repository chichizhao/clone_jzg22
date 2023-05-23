#!/py3
# -- coding: utf-8 --
# author: zhaochichi github: https://github.com/chichizhao/clone_jzg22
# function: remove REJECT snps and indels in the each sample
# usage: python remove_reject.py -i input.vcf.gz -o output.vcf

import sys
import argparse
import gzip

def parse_args():
    parser = argparse.ArgumentParser(description='remove REJECT snps and indels in the each sample')
    parser.add_argument('-i', '--input', help='input vcf file', required=True)
    parser.add_argument('-o', '--output', help='output vcf file', required=True)
    args = parser.parse_args()
    return args
def remove_REJECT(input, output):
    for line in input:
        if line.startswith('#'):
            output.write(line)
        else:
            line = line.strip().split('\t')
            if 'REJECT' in line[6]:
                continue
            else:
                output.write('\t'.join(line) + '\n')
def main():
    args = parse_args()
    input = gzip.open(args.input, 'rt')
    output = open(args.output, 'w')
    remove_REJECT(input, output)
    input.close()
    output.close()

if __name__ == '__main__':
    main()
