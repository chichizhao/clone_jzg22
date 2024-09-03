#!/python
# -- coding: utf-8 --
# author: zhaochichi github: https://github.com/chichizhao/clone_jzg22
#function: calculate the somatic mutation fix number in the in all vcf files 
# usage: python check_fix_somatic_mutaions.py
groupA = [] # add the samples' name in the groupA
import gzip
average = []
num = 0
vcf = {}
for i in groupA:
    vcf1 = gzip.open('path/somaticseq_data/somatic_'+str(i)+'/indels.snvs_filter.vcf.gz', 'rt')
    def read_vcf(vcf):
        mutation = {}
        for line in vcf:
            if line.startswith('#'):
                continue
            elif line.startswith('Chr'):
                #print(line)
                line = line.strip().split('\t')
                pos = line[0]+'_'+line[1]
                genetype = line[9].split(':')[0]
                ref = line[3]
                alt = line[4]    
                af = float(line[7].split('AF=')[1])
                mutation[pos] = [ref, alt,af]   
        return mutation
    mutation = read_vcf(vcf1)
    vcf[i] = mutation
#print(vcf)

# get the farthest distance between two mutation samples of the fix number

max_fix = []
for i in groupA:
    fix = []
    for j in groupA:
        if i != j:
            fix_number = 0
            for k in vcf[i]:
                if vcf[i][k][2] > 0.5 :
                        fix_number += 1
            for k in vcf[j]:
                if vcf[j][k][2] > 0.5 :
                        fix_number += 1
            fix_number = fix_number/2            
            fix.append(fix_number)
    max_fix.append(max(fix))
    print(i, fix)
    print(i, max(fix))
print(max_fix)