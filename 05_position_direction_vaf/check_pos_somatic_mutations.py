
#/py3
# -- coding: utf-8 --
# author: zhaochichi github: https://github.com/chichizhao/clone_jzg22
# function: check the mutation position in gene feature in four categories: Promoter, gene, Transposable elements, Intergenic region
# usage: python3 path/check_pos_somatic_mutation.py
import pandas as pd
import gzip
# read the change head 
# path path/script/change_head.csv
change_header = pd.read_csv('path/script/change_head.csv', sep = "\t", header=0)

# read the CpGs in gff file 
path = "path/met/typha_cpgs.gff"
CpGs ={}
with open(path, "r") as f:
    for line in f:
        if line.startswith("#"):
            continue
        else:
            line = line.strip().split("\t")
            chrom = line[0]
            start = line[3]
            end = line[4]
            if chrom not in CpGs:
                CpGs[chrom] = [start, end]
            else:
                CpGs[chrom].append([start, end])
        
# read the chromesome length file
# path /home/chichi/finalresult/Tl/10.genedenisty/chrlen.txt
chrom_len = pd.read_csv("/home/chichi/finalresult/Tl/10.genedenisty/chrlen.txt", sep = "\t", header=None)

for name in range(1,73):
    old_name = change_header[change_header["change"] == name]["old"].values[0]
    
    # read the mutation file
    # path path/somaticseq_data/somatic_1
    sample = name

    #mutation_file ="path/somaticseq_data/somatic_1/indels.snvs_filter.vcf.gz"
    mutation_file = "path/somaticseq_data/somatic_" + str(old_name) + "/indels.snvs_filter.vcf.gz"
    mutation_site = {}
    mutation_site_num = 0
    with gzip.open(mutation_file, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            else:
                line = line.strip().split("\t")
                chrom = line[0]
                pos = line[1]
                ref = line[3]
                alt = line[4]
                start = int(pos)
                end = int(pos) + 1
                mutation_site_num += 1
                if chrom not in mutation_site:
                    mutation_site[chrom] =[[start, end, ref, alt]]
                else:
                    mutation_site[chrom].append([start, end, ref, alt])

    CpGs_mut = {}
    snp_type_in_CpGs = {}
    num = 0
    num1 = 0
                    
    # find the mutation in the CpGs
    for chrom in mutation_site:
        if chrom in CpGs:
            for i in range(0, len(mutation_site[chrom])):
                start = int(mutation_site[chrom][i][0])
                end = int(mutation_site[chrom][i][1])
                for j in range(0, len(CpGs[chrom])):
                    CpG_start = int(CpGs[chrom][j][0])
                    CpG_end = int(CpGs[chrom][j][1])
                    if start >= CpG_start and end <= CpG_end:
                        if chrom not in CpGs_mut:
                            CpGs_mut[chrom] = [[start, end, mutation_site[chrom][i][2], mutation_site[chrom][i][3]]]
                        else:
                            CpGs_mut[chrom].append([start, end, mutation_site[chrom][i][2], mutation_site[chrom][i][3]])
                        if len(mutation_site[chrom][i][2]) == 1 and len(mutation_site[chrom][i][3]) == 1:
                            num += 1
                            if str(mutation_site[chrom][i][2]) + str(mutation_site[chrom][i][3]) not in snp_type_in_CpGs:
                                snp_type_in_CpGs[str(mutation_site[chrom][i][2]) + str(mutation_site[chrom][i][3])] = 1
                            else:
                                snp_type_in_CpGs[str(mutation_site[chrom][i][2]) + str(mutation_site[chrom][i][3])] += 1
                        else:
                            continue
                        
                    elif CpG_start >= start :
                        break
                    else:
                        continue

    # check the CpGs_mut is valid
    num = 0
    num1 = 0
    for chrom in CpGs_mut:
        element = CpGs_mut[chrom]
        num += len(element)
        element = list(set([tuple(t) for t in element]))
        CpGs_mut[chrom] = element
        num1 += len(element)


    # we need check the mutation happened site if in specific gene structure, like exon, intron, utr, etc
    # read the gene structure file
    # path /home/chichi/finalresult/Tl/03.structure_annotation/pasa2.longest.filter.gff3
    gene_structure = {}
    gene_structure_file = "/home/chichi/finalresult/Tl/03.structure_annotation/pasa2.longest.filter.gff3"
    # we need frist sort the gene structure file by the chromesome then the start site
    gff = pd.read_csv(gene_structure_file, sep = "\t", header=None)

    sort_gff = gff.sort_values(by=[0,3])

    for i in range(0, len(sort_gff)):
        chrom = sort_gff.iloc[i,0]
        start = sort_gff.iloc[i,3]
        end = sort_gff.iloc[i,4]
        gene_type = sort_gff.iloc[i,2]
        if chrom not in gene_structure:
            gene_structure[chrom] = [[start, end, gene_type]]
        else:
            gene_structure[chrom].append([start, end, gene_type])


    num2 = 0
    num4 = 0
    mutation_in_gene_structure = {}

    for chrom in mutation_site:
        if chrom in gene_structure:
            for i in range(0, len(mutation_site[chrom])):
                start = int(mutation_site[chrom][i][0])
                end = int(mutation_site[chrom][i][1])
                ref = mutation_site[chrom][i][2]
                alt = mutation_site[chrom][i][3]
                for j in range(0, len(gene_structure[chrom])):
                    gene_start = int(gene_structure[chrom][j][0])
                    gene_end = int(gene_structure[chrom][j][1])
                    gene_type = gene_structure[chrom][j][2]
                    if start >= gene_start and end <= gene_end:
                        if gene_type == "gene" :
                            '''
                            if chrom not in mutation_in_gene_structure:
                                mutation_in_gene_structure[chrom] = [[start, end, gene_type]]
                                num2 += 1
                            else:
                                mutation_in_gene_structure[chrom].append([start, end, gene_type])
                                num2 += 1
                            '''
                            # as all the gene future are from the gene, so we can make the type according to the gene future
                            # or we the gene in the in tron, exon, utr, etc
                            # and also the gene probably have the overlap, so we also need to check the overlap
                            
                            # first, we need to check the gene feature
                            # get the gene feature
                            gene_feature = []
                            gene_number = 0
                            for k in range(j, len(gene_structure[chrom])):
                                gene_start = int(gene_structure[chrom][k][0])
                                gene_end = int(gene_structure[chrom][k][1])
                                gene_type = gene_structure[chrom][k][2]
                                gene_feature.append([gene_start, gene_end, gene_type])
                                if gene_type == "gene":
                                    gene_number += 1
                                    if gene_number == 2:
                                        gene_next_start = gene_start
                                        gene_next_end = gene_end
                                        gene_next_type = gene_type
                                        break

                            # check the gene feature
                            for k in range(0, len(gene_feature)):
                                gene_start = gene_feature[k][0]
                                gene_end = gene_feature[k][1]
                                gene_type = gene_feature[k][2]
                                if start >= gene_start and end <= gene_end and gene_type != "gene" and gene_type != "mRNA":
                                    if chrom not in mutation_in_gene_structure:
                                        mutation_in_gene_structure[chrom] = [[start, end, gene_type,ref, alt]]
                                        num2 += 1
                                    else:
                                        mutation_in_gene_structure[chrom].append([start, end, gene_type, ref, alt])
                                        num2 += 1
                                    break
                                elif gene_start >= end:
                                    num2 += 1
                                    if chrom not in mutation_in_gene_structure:
                                        mutation_in_gene_structure[chrom] = [[start, end, "intron", ref, alt]]
                                    else:
                                        mutation_in_gene_structure[chrom].append([start, end, "intron", ref, alt])
                                    break
                                else:
                                    continue
                            # check the gene overlap
                            if gene_next_start < start and gene_next_end > end:
                                num4 += 1
                    
                    elif gene_start >= start:
                        break
                    else:
                        continue

    num3 = 0
    for chrom in mutation_in_gene_structure:
        num3 += len(mutation_in_gene_structure[chrom])

    gene = num3

    # check the mutation in the gene is valid
    num = 0
    for chrom in mutation_in_gene_structure:
        element = mutation_in_gene_structure[chrom]
        element = list(set([tuple(t) for t in element]))
        mutation_in_gene_structure[chrom] = element
        num += len(element)

    # plot the mutation in the chromesome like way

    # check the mutation in specific region of TEprotein Transposon  and TRF
    # read the annotation gff file
    # path /home/chichi/finalresult/Tl/11.plots/all.gff
    annotation_gff = "/home/chichi/finalresult/Tl/11.plots/all.gff"
    TEprotein = {}
    TRF = {}
    Transposon = {}
    for line in open(annotation_gff, "r"):
        if line.startswith("#"):
            continue
        else:
            line = line.strip().split("\t")
            chrom = line[0]
            start = line[3]
            end = line[4]
            feature = line[2]
            if feature == "TEprotein":
                if chrom not in TEprotein:
                    TEprotein[chrom] = [[start, end]]
                else:
                    TEprotein[chrom].append([start, end])
            elif feature == "TandemRepeat":
                if chrom not in TRF:
                    TRF[chrom] = [[start, end]]
                else:
                    TRF[chrom].append([start, end])
            elif feature == "Transposon":
                if chrom not in Transposon:
                    Transposon[chrom] = [[start, end]]
                else:
                    Transposon[chrom].append([start, end])
            else:
                continue
            

    # check the mutation in the TEprotein, TRF, Transposon
    num5 = 0
    num6 = 0
    num7 = 0
    mutation_in_TEprotein = {}
    mutation_in_TRF = {}
    mutation_in_Transposon = {}
    for chrom in mutation_site:
        if chrom in TEprotein:
            for i in range(0, len(mutation_site[chrom])):
                start = int(mutation_site[chrom][i][0])
                end = int(mutation_site[chrom][i][1])
                ref = mutation_site[chrom][i][2]
                alt = mutation_site[chrom][i][3]

                for j in range(0, len(TEprotein[chrom])):
                    TEprotein_start = int(TEprotein[chrom][j][0])
                    TEprotein_end = int(TEprotein[chrom][j][1])
                    if start >= TEprotein_start and end <= TEprotein_end:
                        if chrom not in mutation_in_TEprotein:
                            mutation_in_TEprotein[chrom] = [[start, end, ref, alt]]
                            num5 += 1
                        else:
                            mutation_in_TEprotein[chrom].append([start, end, ref, alt])
                            num5 += 1
                        # as the TEproteins probably have the overlap
                        # but once we consider this site as in TEprotein, we stop the loop
                        break

                    elif TEprotein_start >= end:
                        break
                    else:
                        continue

        if chrom in TRF:
            for i in range(0, len(mutation_site[chrom])):
                start = int(mutation_site[chrom][i][0])
                end = int(mutation_site[chrom][i][1])
                ref = mutation_site[chrom][i][2]
                alt = mutation_site[chrom][i][3]
                for j in range(0, len(TRF[chrom])):
                    TRF_start = int(TRF[chrom][j][0])
                    TRF_end = int(TRF[chrom][j][1])
                    if start >= TRF_start and end <= TRF_end:
                        if chrom not in mutation_in_TRF:
                            mutation_in_TRF[chrom] = [[start, end, ref, alt]]
                            num6 += 1
                        else:
                            mutation_in_TRF[chrom].append([start, end, ref, alt])
                            num6 += 1
                        # same as the TEprotein, we need to stop the loop
                        break
                    elif TRF_start >= end:
                        break
                    else:
                        continue
        if chrom in Transposon:
            for i in range(0, len(mutation_site[chrom])):
                start = int(mutation_site[chrom][i][0])
                end = int(mutation_site[chrom][i][1])
                ref = mutation_site[chrom][i][2]
                alt = mutation_site[chrom][i][3]
                for j in range(0, len(Transposon[chrom])):
                    Transposon_start = int(Transposon[chrom][j][0])
                    Transposon_end = int(Transposon[chrom][j][1])
                    if start >= Transposon_start and end <= Transposon_end:
                        if chrom not in mutation_in_Transposon:
                            mutation_in_Transposon[chrom] = [[start, end, ref, alt]]
                            num7 += 1
                        else:
                            mutation_in_Transposon[chrom].append([start, end, ref, alt])
                            num7 += 1
                        # same as the TEprotein, we need to stop the loop
                        break
                    elif Transposon_start >= end:
                        break
                    else:
                        continue
    transposon = num7

    # check the mutation in the promoter region (2000bp upstream of the gene)

    num8 = 0
    mutation_in_promoter = {}
    for chrom in mutation_site:
        if chrom in gene_structure:
            for i in range(0, len(mutation_site[chrom])):
                start = int(mutation_site[chrom][i][0])
                end = int(mutation_site[chrom][i][1])
                ref = mutation_site[chrom][i][2]
                alt = mutation_site[chrom][i][3]
                for j in range(0, len(gene_structure[chrom])):
                    gene_start = int(gene_structure[chrom][j][0])
                    gene_end = int(gene_structure[chrom][j][1])
                    gene_type = gene_structure[chrom][j][2]
                    if gene_type == "gene":
                        if gene_start - 2000 <= start and gene_start >= start:
                            if chrom not in mutation_in_promoter:
                                mutation_in_promoter[chrom] = [[start, end, ref, alt]]
                                num8 += 1
                            else:
                                mutation_in_promoter[chrom].append([start, end, ref, alt])
                                num8 += 1
                                break
                            break
                        elif gene_start - 2000 >= start:
                            break
                        else:
                            continue
    promoter = num8
    # check the overlap of the mutation in the TEprotein and transposon
    num8 = 0
    num9 = 0
    for chrom in mutation_in_TEprotein:
        for i in range(0, len(mutation_in_TEprotein[chrom])):
            start = int(mutation_in_TEprotein[chrom][i][0])
            end = int(mutation_in_TEprotein[chrom][i][1])
            ref = mutation_in_TEprotein[chrom][i][2]
            alt = mutation_in_TEprotein[chrom][i][3]
            if chrom in mutation_in_Transposon:
                for j in range(0, len(mutation_in_Transposon[chrom])):
                    Transposon_start = int(mutation_in_Transposon[chrom][j][0])
                    Transposon_end = int(mutation_in_Transposon[chrom][j][1])
                    Transposon_ref = mutation_in_Transposon[chrom][j][2]
                    Transposon_alt = mutation_in_Transposon[chrom][j][3]
                    if start == Transposon_start and end == Transposon_end and ref == Transposon_ref and alt == Transposon_alt:
                        num8 += 1
                        break
            else:
                continue
                
    for chrom in mutation_in_TRF:
        for i in range(0, len(mutation_in_TRF[chrom])):
            start = int(mutation_in_TRF[chrom][i][0])
            end = int(mutation_in_TRF[chrom][i][1])
            ref = mutation_in_TRF[chrom][i][2]
            alt = mutation_in_TRF[chrom][i][3]
            if chrom in mutation_in_Transposon:
                for j in range(0, len(mutation_in_Transposon[chrom])):
                    Transposon_start = int(mutation_in_Transposon[chrom][j][0])
                    Transposon_end = int(mutation_in_Transposon[chrom][j][1])
                    Transposon_ref = mutation_in_Transposon[chrom][j][2]
                    Transposon_alt = mutation_in_Transposon[chrom][j][3]
                    if start == Transposon_start and end == Transposon_end and ref == Transposon_ref and alt == Transposon_alt:
                        num9 += 1
                        break
            else:
                continue

    # check the overlap of the mutation in the gene and TEprotein, TRF, Transposon
    num10 = 0

    for chrom in mutation_in_gene_structure:
        for i in range(0, len(mutation_in_gene_structure[chrom])):
            start = int(mutation_in_gene_structure[chrom][i][0])
            end = int(mutation_in_gene_structure[chrom][i][1])
            ref = mutation_in_gene_structure[chrom][i][3]
            alt = mutation_in_gene_structure[chrom][i][4]
            if chrom in mutation_in_TEprotein:
                for j in range(0, len(mutation_in_Transposon[chrom])):
                    Transposon_start = int(mutation_in_Transposon[chrom][j][0])
                    Transposon_end = int(mutation_in_Transposon[chrom][j][1])
                    Transposon_ref = mutation_in_Transposon[chrom][j][2]
                    Transposon_alt = mutation_in_Transposon[chrom][j][3]
                    if start == Transposon_start and end == Transposon_end and ref == Transposon_ref and alt == Transposon_alt:
                        num10 += 1
                        break
            else:
                continue
    gene_transposon = num10
    print(num10)

    # check the overlap of in promoter and Transposon
    num11 = 0

    for chrom in mutation_in_promoter:
        for i in range(0, len(mutation_in_promoter[chrom])):
            start = int(mutation_in_promoter[chrom][i][0])
            end = int(mutation_in_promoter[chrom][i][1])
            ref = mutation_in_promoter[chrom][i][2]
            alt = mutation_in_promoter[chrom][i][3]
            if chrom in mutation_in_Transposon:
                for j in range(0, len(mutation_in_Transposon[chrom])):
                    Transposon_start = int(mutation_in_Transposon[chrom][j][0])
                    Transposon_end = int(mutation_in_Transposon[chrom][j][1])
                    Transposon_ref = mutation_in_Transposon[chrom][j][2]
                    Transposon_alt = mutation_in_Transposon[chrom][j][3]
                    if start == Transposon_start and end == Transposon_end and ref== Transposon_ref and alt == Transposon_alt:
                        num11 += 1
                        break
            else:
                continue
    promoter_transposon = num11

    # check the overlap of the promoter and gene
    num12 = 0
    for chrom in mutation_in_promoter:
        for i in range(0, len(mutation_in_promoter[chrom])):
            start = int(mutation_in_promoter[chrom][i][0])
            end = int(mutation_in_promoter[chrom][i][1])
            ref = mutation_in_promoter[chrom][i][2]
            alt = mutation_in_promoter[chrom][i][3]
            if chrom in mutation_in_gene_structure:
                for j in range(0, len(mutation_in_gene_structure[chrom])):
                    gene_start = int(mutation_in_gene_structure[chrom][j][0])
                    gene_end = int(mutation_in_gene_structure[chrom][j][1])
                    gene_ref = mutation_in_gene_structure[chrom][j][3]
                    gene_alt = mutation_in_gene_structure[chrom][j][4]
                    if start == gene_start and end == gene_end and ref == gene_ref and alt == gene_alt:
                        num12 += 1
                        break
            else:
                continue
    promoter_gene = num12

    # check the overlap of the mutation in the promoter, gene, and transposon
    num13 = 0
    for chrom in mutation_in_promoter:
        for i in range(0, len(mutation_in_promoter[chrom])):
            start = int(mutation_in_promoter[chrom][i][0])
            end = int(mutation_in_promoter[chrom][i][1])
            ref = mutation_in_promoter[chrom][i][2]
            alt = mutation_in_promoter[chrom][i][3]
            if chrom in mutation_in_gene_structure:
                for j in range(0, len(mutation_in_gene_structure[chrom])):
                    gene_start = int(mutation_in_gene_structure[chrom][j][0])
                    gene_end = int(mutation_in_gene_structure[chrom][j][1])
                    gene_ref = mutation_in_gene_structure[chrom][j][3]
                    gene_alt = mutation_in_gene_structure[chrom][j][4]
                    if start == gene_start and end == gene_end and ref == gene_ref and alt == gene_alt:
                        if chrom in mutation_in_Transposon:
                            for k in range(0, len(mutation_in_Transposon[chrom])):
                                Transposon_start = int(mutation_in_Transposon[chrom][k][0])
                                Transposon_end = int(mutation_in_Transposon[chrom][k][1])
                                Transposon_ref = mutation_in_Transposon[chrom][k][2]
                                Transposon_alt = mutation_in_Transposon[chrom][k][3]
                                if start == Transposon_start and end == Transposon_end and ref == Transposon_ref and alt == Transposon_alt:
                                    num13 += 1
                                    break
                        else:
                            continue
                            
            else:
                continue
    promoter_gene_transposon = num13
    promoter_gene_out_transposon = promoter_gene - promoter_gene_transposon
    promoter_transposon_out_gene = promoter_transposon - promoter_gene_transposon
    promoter_out_gene_transposon = promoter - promoter_gene_transposon - promoter_gene_out_transposon - promoter_transposon_out_gene
    gene_transposon_out_promoter = gene_transposon - promoter_gene_transposon
    gene_out_promoter_transposon = gene- promoter_gene_transposon - gene_transposon_out_promoter- promoter_gene_out_transposon
    transposon_out_promoter_gene = transposon - promoter_gene_transposon-gene_transposon_out_promoter-promoter_transposon_out_gene


    # check the overlap of the TEprotein, Transposon, gene, promoter
    num14 = 0
    for chrom in mutation_in_TEprotein:
        for i in range(0, len(mutation_in_TEprotein[chrom])):
            start = int(mutation_in_TEprotein[chrom][i][0])
            end = int(mutation_in_TEprotein[chrom][i][1])
            ref = mutation_in_TEprotein[chrom][i][2]
            alt = mutation_in_TEprotein[chrom][i][3]
            control_element = 0
            if chrom in mutation_in_Transposon:
                for j in range(0, len(mutation_in_Transposon[chrom])):
                    Transposon_start = int(mutation_in_Transposon[chrom][j][0])
                    Transposon_end = int(mutation_in_Transposon[chrom][j][1])
                    Transposon_ref = mutation_in_Transposon[chrom][j][2]
                    Transposon_alt = mutation_in_Transposon[chrom][j][3]
                    if start == Transposon_start and end == Transposon_end and ref == Transposon_ref and alt == Transposon_alt:
                        control_element = 1
                        break
                    elif Transposon_start >= end:
                        break
                    else:
                        continue
            if control_element != 1:
                for j in range(0, len(mutation_in_gene_structure[chrom])):
                    gene_start = int(mutation_in_gene_structure[chrom][j][0])
                    gene_end = int(mutation_in_gene_structure[chrom][j][1])
                    gene_ref = mutation_in_gene_structure[chrom][j][2]
                    gene_alt = mutation_in_gene_structure[chrom][j][3]
                    if start == gene_start and end == gene_end and ref == gene_ref and alt == gene_alt:
                        control_element = 2
                        break
                    elif gene_start >= end:
                        break
                    else:
                        continue
            if control_element != 1 and control_element != 2:
                for j in range(0, len(mutation_in_promoter[chrom])):
                    promoter_start = int(mutation_in_promoter[chrom][j][0])
                    promoter_end = int(mutation_in_promoter[chrom][j][1])
                    promoter_ref = mutation_in_promoter[chrom][j][2]
                    promoter_alt = mutation_in_promoter[chrom][j][3]
                    if start == promoter_start and end == promoter_end and ref == promoter_ref and alt == promoter_alt:
                        control_element = 3
                        break
                    elif promoter_start >= end:
                        break
                    else:
                        continue
            if control_element != 1 and control_element != 2 and control_element != 3:
                num14 += 1

    TEprotein_out = num14


    # check the overlap of the TRF, Transposon, gene, promoter
    num15 = 0
    num15_1 = 0
    num15_2 = 0
    num15_3 = 0
    for chrom in mutation_in_TRF:
        for i in range(0, len(mutation_in_TRF[chrom])):
            start = int(mutation_in_TRF[chrom][i][0])
            end = int(mutation_in_TRF[chrom][i][1])
            ref = mutation_in_TRF[chrom][i][2]
            alt = mutation_in_TRF[chrom][i][3]
            control_element = 0
            if chrom in mutation_in_Transposon:
                for j in range(0, len(mutation_in_Transposon[chrom])):
                    Transposon_start = int(mutation_in_Transposon[chrom][j][0])
                    Transposon_end = int(mutation_in_Transposon[chrom][j][1])
                    Transposon_ref = mutation_in_Transposon[chrom][j][2]
                    Transposon_alt = mutation_in_Transposon[chrom][j][3]
                    if start == Transposon_start and end == Transposon_end and ref == Transposon_ref and alt == Transposon_alt:
                        control_element = 1
                        break
                    elif Transposon_start >= end:
                        break
                    else:
                        continue
            if control_element != 1:
                if chrom in mutation_in_gene_structure:
                    for j in range(0, len(mutation_in_gene_structure[chrom])):
                        gene_start = int(mutation_in_gene_structure[chrom][j][0])
                        gene_end = int(mutation_in_gene_structure[chrom][j][1])
                        gene_ref = mutation_in_gene_structure[chrom][j][2]
                        gene_alt = mutation_in_gene_structure[chrom][j][3]
                        if start == gene_start and end == gene_end and ref == gene_ref and alt == gene_alt:
                            control_element = 2
                            num15_2 += 1
                            break
                        elif gene_start >= end:
                            break
                        else:
                            continue
                else:
                    continue
            if control_element != 1 and control_element != 2:
                if chrom in mutation_in_promoter:
                    for j in range(0, len(mutation_in_promoter[chrom])):    
                        promoter_start = int(mutation_in_promoter[chrom][j][0])
                        promoter_end = int(mutation_in_promoter[chrom][j][1])
                        promoter_ref = mutation_in_promoter[chrom][j][2]
                        promoter_alt = mutation_in_promoter[chrom][j][3]
                        if start == promoter_start and end == promoter_end and ref == promoter_ref and alt == promoter_alt:
                            control_element = 3
                            num15_3 += 1
                            break
                        elif promoter_start >= end:
                            break
                        else:
                            continue
                else:
                    continue
            if control_element != 1 and control_element != 2 and control_element != 3:
                num15 += 1


    # check the overlap between CpGs, promoter, gene,Treansposon，TRF
    num16 = 0
    num16_1 = 0
    num16_2 = 0
    num16_3 = 0
    num16_4 = 0
    for chrom in CpGs_mut:
        for i in range(0, len(CpGs_mut[chrom])):
            start = int(CpGs_mut[chrom][i][0])
            end = int(CpGs_mut[chrom][i][1])
            ref = CpGs_mut[chrom][i][2]
            alt = CpGs_mut[chrom][i][3]
            control_element = 0
            if chrom in mutation_in_Transposon:
                for j in range(0, len(mutation_in_Transposon[chrom])):
                    Transposon_start = int(mutation_in_Transposon[chrom][j][0])
                    Transposon_end = int(mutation_in_Transposon[chrom][j][1])
                    Transposon_ref = mutation_in_Transposon[chrom][j][2]
                    Transposon_alt = mutation_in_Transposon[chrom][j][3]
                    if start == Transposon_start and end == Transposon_end and ref == Transposon_ref and alt == Transposon_alt:
                        control_element = 1
                        num16_1 += 1
                        break
                    elif Transposon_start >= end:
                        break
                    else:
                        continue
            if control_element != 1:
                if chrom in mutation_in_gene_structure:
                    for j in range(0, len(mutation_in_gene_structure[chrom])):
                        gene_start = int(mutation_in_gene_structure[chrom][j][0])
                        gene_end = int(mutation_in_gene_structure[chrom][j][1])
                        gene_ref = mutation_in_gene_structure[chrom][j][2]
                        gene_alt = mutation_in_gene_structure[chrom][j][3]
                        if start == gene_start and end == gene_end and ref == gene_ref and alt == gene_alt:
                            control_element = 2
                            num16_2 += 1
                            break
                        elif gene_start >= end:
                            break
                        else:
                            continue
                else:
                    continue
            if control_element != 1 and control_element != 2:
                if chrom in mutation_in_promoter:
                    for j in range(0, len(mutation_in_promoter[chrom])):
                        promoter_start = int(mutation_in_promoter[chrom][j][0])
                        promoter_end = int(mutation_in_promoter[chrom][j][1])
                        promoter_ref = mutation_in_promoter[chrom][j][2]
                        promoter_alt = mutation_in_promoter[chrom][j][3]
                        if start == promoter_start and end == promoter_end and ref == promoter_ref and alt == promoter_alt:
                            control_element = 3
                            num16_3 += 1
                            break
                        elif promoter_start >= end:
                            break
                        else:
                            continue
                else:
                    continue
            if control_element != 1 and control_element != 2 and control_element != 3:
                if chrom in mutation_in_TRF:
                    for j in range(0, len(mutation_in_TRF[chrom])):
                        TRF_start = int(mutation_in_TRF[chrom][j][0])
                        TRF_end = int(mutation_in_TRF[chrom][j][1])
                        TRF_ref = mutation_in_TRF[chrom][j][2]
                        TRF_alt = mutation_in_TRF[chrom][j][3]
                        if start == TRF_start and end == TRF_end and ref == TRF_ref and alt == TRF_alt:
                            control_element = 4
                            num16_4 += 1
                            break
                        elif TRF_start >= end:
                            break
                        else:
                            continue
                else:
                    continue
            if control_element != 1 and control_element != 2 and control_element != 3 and control_element != 4:

                num16 += 1



    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    import numpy as np
    import random

    fig, ax = plt.subplots()
    fig.set_size_inches(16, 10)
    ax.set_xlim(0, 24800000)
    ax.set_ylim(0, 15.2*1000000)
    # plot the chromesome
    ax.axis('off')
    num10 = 0
    #for i in range(len(chrom_len)):
    for i in range(15):

        chrom = chrom_len.iloc[i,0]
        length = chrom_len.iloc[i,1]
        ax.add_patch(patches.Rectangle((0, i*1000000), length, 200000, facecolor = "none", edgecolor = "black", linewidth=0.5))
        ax.add_patch(patches.Rectangle((0, (i+0.3)*1000000), length, 200000, facecolor = "none", edgecolor = "black", linewidth=0.5))
        ax.add_patch(patches.Rectangle((0, (i+0.6)*1000000), length, 150000, facecolor = "none", edgecolor = "black", linewidth=0.5))
        ax.text(-600000, (i+0.4)*1000000, chrom, ha='center', va='center', fontsize=12)

    
        # plot the TEs(TEprotein, Transposon)
        if chrom in TEprotein:
            for j in range(0, len(TEprotein[chrom])):
                start = TEprotein[chrom][j][0]
                end = TEprotein[chrom][j][1]
                ax.add_patch(patches.Rectangle((int(start), i*1000000), int(end)-int(start), 200000, facecolor = "green", edgecolor = "green",linewidth=0.01,alpha=0.5))
        if chrom in Transposon:
            for j in range(0, len(Transposon[chrom])):
                start = Transposon[chrom][j][0]
                end = Transposon[chrom][j][1]
                ax.add_patch(patches.Rectangle((int(start), i*1000000), int(end)-int(start), 200000, facecolor = "red", edgecolor = "red",linewidth=0.01,alpha=0.5))
        
        
        
        # plot the gene structure
        if chrom in gene_structure:
            for j in range(0, len(gene_structure[chrom])):
                start = gene_structure[chrom][j][0]
                end = gene_structure[chrom][j][1]
                type = gene_structure[chrom][j][2]
                if type == "gene":
                    ax.add_patch(patches.Rectangle((int(start), (i+0.3)*1000000), int(end)-int(start), 200000, facecolor = "blue", edgecolor = "blue",linewidth=0.01,alpha=0.1))
                elif type == "exon":
                    ax.add_patch(patches.Rectangle((int(start), (i+0.3)*1000000), int(end)-int(start), 200000, facecolor = "green", edgecolor = "green",linewidth=0.01))
                elif type == "intron":
                    ax.add_patch(patches.Rectangle((int(start), (i+0.3)*1000000), int(end)-int(start), 200000, facecolor = "yellow", edgecolor = "yellow",linewidth=0.01))
                elif type == "five_prime_UTR":
                    ax.add_patch(patches.Rectangle((int(start), (i+0.3)*1000000), int(end)-int(start), 200000, facecolor = "black", edgecolor = "black",linewidth=0.1))
                #else:
                #    ax.add_patch(patches.Rectangle((int(start), (i+0.3)*1000000), int(end)-int(start), 200000, facecolor = "red", edgecolor = "red",linewidth=0.1))
        
        # add the snp in the sample
        # 1. add the mutation in all
        if chrom in mutation_site:
            for j in range(0, len(mutation_site[chrom])):
                start = mutation_site[chrom][j][0]
                end = mutation_site[chrom][j][1]
                ax.add_patch(patches.Rectangle((int(start), (i+0.6)*1000000), int(end)-int(start), 150000, facecolor = "grey", edgecolor = "grey",linewidth=0.5, alpha=0.5))
                
        # 2. add the mutation in the promoter
        
        if chrom in mutation_in_promoter:
            for j in range(0, len(mutation_in_promoter[chrom])):
                start = mutation_in_promoter[chrom][j][0]
                end = mutation_in_promoter[chrom][j][1]
                ax.add_patch(patches.Rectangle((int(start), (i+0.6)*1000000), int(end)-int(start), 150000, facecolor = "yellow", edgecolor = "yellow",linewidth=0.5, alpha=0.5))
                
        # 3. add the mutation in the gene
        if chrom in mutation_in_gene_structure:
            for j in range(0, len(mutation_in_gene_structure[chrom])):
            
                if mutation_in_gene_structure[chrom][j][2] == "exon":
                    ax.add_patch(patches.Rectangle((int(start), (i+0.6)*1000000), int(end)-int(start), 150000, facecolor = "green", edgecolor = "green",linewidth=0.5, alpha=0.5))   
                elif mutation_in_gene_structure[chrom][j][2] == "intron":
                    ax.add_patch(patches.Rectangle((int(start), (i+0.6)*1000000), int(end)-int(start), 150000, facecolor = "purple", edgecolor = "purple",linewidth=0.5, alpha=0.5))
                elif mutation_in_gene_structure[chrom][j][2] == "five_prime_UTR":
                    ax.add_patch(patches.Rectangle((int(start), (i+0.6)*1000000), int(end)-int(start), 150000, facecolor = "black", edgecolor = "black",linewidth=0.5, alpha=0.5))
        
        # 4. add the mutation in the TEprotein
        if chrom in mutation_in_TEprotein:
            for j in range(0, len(mutation_in_TEprotein[chrom])):
                start = mutation_in_TEprotein[chrom][j][0]
                end = mutation_in_TEprotein[chrom][j][1]
                ax.add_patch(patches.Rectangle((int(start), (i+0.6)*1000000), int(end)-int(start), 150000, facecolor = "blue", edgecolor = "blue",linewidth=0.5, alpha=0.5))
        
        # finally, add the mutation in the transposon
        if chrom in mutation_in_Transposon:
            for j in range(0, len(mutation_in_Transposon[chrom])):
                start = mutation_in_Transposon[chrom][j][0]
                end = mutation_in_Transposon[chrom][j][1]
                ax.add_patch(patches.Rectangle((int(start), (i+0.6)*1000000), int(end)-int(start), 150000, facecolor = "red", edgecolor = "red",linewidth=0.5, alpha=0.5))
            
        # no we would mark the muatations that in gene 

        if chrom in mutation_in_gene_structure:
            for j in range(0, len(mutation_in_gene_structure[chrom])):
                start = mutation_in_gene_structure[chrom][j][0]
                end = mutation_in_gene_structure[chrom][j][1]
                type = mutation_in_gene_structure[chrom][j][2]
                #ax.add_patch(patches.Rectangle((int(start), i*1000000), int(end)-int(start), 0.5, facecolor = "blue", edgecolor = "blue",linewidth=0.1))
                # we need give the type different mark for position
                if type == "exon":
                    ax.add_patch(patches.Rectangle((int(start), (i+0.6)*1000000), int(end)-int(start), 200000, facecolor = "green", edgecolor = "green",linewidth=0.1))
                    # add the mark box 
                    ax.add_patch(patches.Rectangle((int(start)-25000, (i+0.8)*1000000), 50000, 100000, facecolor = "green", edgecolor = "green",linewidth=0.5))
                    
                #elif type == "intron":
                    #ax.add_patch(patches.Rectangle((int(start), (i+0.6)*1000000), int(end)-int(start), 200000, facecolor = "yellow", edgecolor = "yellow",linewidth=0.1))
                elif type == "five_prime_UTR":
                    ax.add_patch(patches.Rectangle((int(start), (i+0.6)*1000000), int(end)-int(start), 200000, facecolor = "black", edgecolor = "black",linewidth=0.1))
                    # add the mark box
                    ax.add_patch(patches.Rectangle((int(start)-25000, (i+0.8)*1000000), 50000, 100000, facecolor = "black", edgecolor = "black",linewidth=0.5))
                    num10 += 1
                #else:
                    #ax.add_patch(patches.Rectangle((int(start), (i+0.6)*1000000), int(end)-int(start), 200000, facecolor = "red", edgecolor = "red",linewidth=0.1))
                    # add the mark box
                    #ax.add_patch(patches.Rectangle((int(start)-50000, (i+0.8)*1000000), 100000, 100000, facecolor = "red", edgecolor = "red",linewidth=0.5))
        
        ax.add_patch(patches.Rectangle((0, i*1000000), length, 200000, facecolor = "none", edgecolor = "black",linewidth=0.5))
        ax.add_patch(patches.Rectangle((0, (i+0.3)*1000000), length, 200000, facecolor = "none", edgecolor = "black",linewidth=0.5))
        ax.add_patch(patches.Rectangle((0, (i+0.6)*1000000), length, 150000, facecolor = "none", edgecolor = "black",linewidth=0.5))
        if i == 2:
            plt.plot([length,19600000],[2750000,8100000],linewidth = 1, color = "black", alpha = 0.5)
            plt.plot([length,19600000],[2600000,3830000],linewidth = 1, color = "black", alpha = 0.5)
            plt.plot([length,19600000],[2500000,3700000],linewidth = 1, color = "black", alpha = 0.5)
            plt.plot([length,19600000],[2300000,2530000],linewidth = 1, color = "black", alpha = 0.5)
            plt.plot([length,19600000],[2200000,2320000],linewidth = 1, color = "black", alpha = 0.5)
            plt.plot([length,19600000],[2000000,1550000],linewidth = 1, color = "black", alpha = 0.5)

        
    # add the veen plot of the mutation location distribution
    # add the box of total mutation # chrom_len = 17
    ax.add_patch(patches.Rectangle((17200000,9500000), 7000000, 5500000, facecolor = "white", edgecolor = "black",linewidth=0.5))
    ax.add_patch(patches.Rectangle((17200000,9500000), 7000000, 5500000, facecolor = "grey", edgecolor = "black",linewidth=0.5,alpha=0.05))
    # add the three circle of venn plot
    ax.add_patch(patches.Circle((18500000, 13000000), 1200000, facecolor = "yellow", edgecolor = "yellow",linewidth=0.5,alpha=0.2))
    ax.add_patch(patches.Circle((18500000, 13000000), 1200000, facecolor = "none", edgecolor = "black",linewidth=0.5))
    ax.add_patch(patches.Circle((19900000, 13000000), 1200000, facecolor = "green", edgecolor = "green",linewidth=0.5,alpha=0.2))
    ax.add_patch(patches.Circle((19900000, 13000000), 1200000, facecolor = "none", edgecolor = "green",linewidth=0.5))
    ax.add_patch(patches.Circle((19200000, 12000000), 1200000, facecolor = "red", edgecolor = "red",linewidth=0.5,alpha=0.2))
    ax.add_patch(patches.Circle((19200000, 12000000), 1200000, facecolor = "none", edgecolor = "red",linewidth=0.5))
    # add the number of the mutation in the venn plot
    ax.text(18200000, 13100000, str(promoter_out_gene_transposon), ha='center', va='center', fontsize=10)
    ax.text(18100000, 13400000, str("Promoter"), ha='center', va='center', fontsize=10)
    ax.text(20200000, 13100000, str(gene_out_promoter_transposon), ha='center', va='center', fontsize=10)
    ax.text(20200000, 13400000, str("Gene"), ha='center', va='center', fontsize=10)
    ax.text(19200000, 11500000, str(transposon_out_promoter_gene + TEprotein_out), ha='center', va='center', fontsize=10)
    ax.text(19200000, 11200000, str("TEs"), ha='center', va='center', fontsize=10)
    ax.text(19200000, 12600000, str(promoter_gene_transposon), ha='center', va='center', fontsize=10)
    ax.text(19200000, 13350000, str(promoter_gene_out_transposon), ha='center', va='center', fontsize=10)
    ax.text(18500000, 12300000, str(promoter_transposon_out_gene), ha='center', va='center', fontsize=10)
    ax.text(19900000, 12300000, str(gene_transposon_out_promoter), ha='center', va='center', fontsize=10)
    # add one circle of the TEprotein
    #ax.add_patch(patches.Circle((22000000, 11500000), 800000, facecolor = "pink", edgecolor = "pink",linewidth=0.5,alpha=0.2))
    #ax.add_patch(patches.Circle((22000000, 11500000), 800000, facecolor = "none", edgecolor = "black",linewidth=0.5))
    #ax.text(22000000, 11500000, str(TEprotein_out), ha='center', va='center', fontsize=10)
    #ax.text(22000000, 11200000, str("TE protein"), ha='center', va='center', fontsize=10)
    IGR = mutation_site_num - promoter - gene_transposon_out_promoter - gene_out_promoter_transposon - transposon_out_promoter_gene - TEprotein_out 
    ax.text(21800000, 11500000, str("IGR:"), ha='center', va='center', fontsize=10)
    ax.text(22400000, 11500000, str(IGR), ha='center', va='center', fontsize=10)
    ax.text(21400000, 10000000, str("Sample"+str(sample)+ ":"), ha='center', va='center', fontsize=12)
    ax.text(23000000, 10000000, str(mutation_site_num), ha='center', va='center', fontsize=12)

    # add the legend box    
    ax.add_patch(patches.Rectangle((19500000, 1500000), 5500000, 7000000, facecolor = "white", edgecolor = "white",linewidth=0.5))

    ax.add_patch(patches.Rectangle((19700000, 7400000), 900000, 200000, facecolor = "none", edgecolor = "black",linewidth=0.5))
    ax.add_patch(patches.Rectangle((20010000, 7400000), 20000, 400000, facecolor = "green", edgecolor = "green",linewidth=0.5))
    ax.add_patch(patches.Rectangle((19920000, 7650000), 200000, 250000, facecolor = "green", edgecolor = "green",linewidth=0.5))


    ax.add_patch(patches.Rectangle((19700000, 6800000), 900000, 200000, facecolor = "none", edgecolor = "black",linewidth=0.5))
    ax.add_patch(patches.Rectangle((20310000, 6800000), 20000, 400000, facecolor = "black", edgecolor = "black",linewidth=0.5))
    ax.add_patch(patches.Rectangle((20220000, 7050000), 200000, 250000, facecolor = "black", edgecolor = "black",linewidth=0.5))


    ax.add_patch(patches.Rectangle((19700000, 6300000), 900000, 200000, facecolor = "none", edgecolor = "black",linewidth=0.5))
    ax.add_patch(patches.Rectangle((20200000, 6300000), 50000, 200000, facecolor = "#CCCC00", edgecolor = "#CCCC00",linewidth=0.5,alpha=0.5))
    ax.add_patch(patches.Rectangle((19700000, 6300000), 900000, 200000, facecolor = "none", edgecolor = "black",linewidth=0.5))

    ax.add_patch(patches.Rectangle((19700000, 5900000), 900000, 200000, facecolor = "none", edgecolor = "black",linewidth=0.5))
    ax.add_patch(patches.Rectangle((20200000, 5900000), 50000, 200000, facecolor = "green", edgecolor = "green",linewidth=0.5,alpha=0.5))
    ax.add_patch(patches.Rectangle((19700000, 5900000), 900000, 200000, facecolor = "none", edgecolor = "black",linewidth=0.5))

    ax.add_patch(patches.Rectangle((19700000, 5500000), 900000, 200000, facecolor = "none", edgecolor = "black",linewidth=0.5))
    ax.add_patch(patches.Rectangle((20200000, 5500000), 50000, 200000, facecolor = "blue", edgecolor = "blue",linewidth=0.5,alpha=0.5))
    ax.add_patch(patches.Rectangle((19700000, 5500000), 900000, 200000, facecolor = "none", edgecolor = "black",linewidth=0.5))

    ax.add_patch(patches.Rectangle((19700000, 5100000), 900000, 200000, facecolor = "none", edgecolor = "black",linewidth=0.5))
    ax.add_patch(patches.Rectangle((20200000, 5100000), 50000, 200000, facecolor = "black", edgecolor = "black",linewidth=0.5,alpha=0.5))
    ax.add_patch(patches.Rectangle((19700000, 5100000), 900000, 200000, facecolor = "none", edgecolor = "black",linewidth=0.5))

    ax.add_patch(patches.Rectangle((19700000, 4700000), 900000, 200000, facecolor = "none", edgecolor = "black",linewidth=0.5))
    ax.add_patch(patches.Rectangle((20200000, 4700000), 50000, 200000, facecolor = "pink", edgecolor = "pink",linewidth=0.5,alpha=0.5))
    ax.add_patch(patches.Rectangle((19700000, 4700000), 900000, 200000, facecolor = "none", edgecolor = "black",linewidth=0.5))

    ax.add_patch(patches.Rectangle((19700000, 4300000), 900000, 200000, facecolor = "none", edgecolor = "black",linewidth=0.5))
    ax.add_patch(patches.Rectangle((20200000, 4300000), 50000, 200000, facecolor = "red", edgecolor = "red",linewidth=0.5,alpha=0.5))
    ax.add_patch(patches.Rectangle((19700000, 4300000), 900000, 200000, facecolor = "none", edgecolor = "black",linewidth=0.5))

    ax.add_patch(patches.Rectangle((19700000, 3900000), 900000, 200000, facecolor = "none", edgecolor = "black",linewidth=0.5))
    ax.add_patch(patches.Rectangle((20200000, 3900000), 50000, 200000, facecolor = "grey", edgecolor = "grey",linewidth=0.5,alpha=0.5))
    ax.add_patch(patches.Rectangle((19700000, 3900000), 900000, 200000, facecolor = "none", edgecolor = "black",linewidth=0.5))

    ax.add_patch(patches.Rectangle((19700000, 3400000), 900000, 200000, facecolor = "none", edgecolor = "black",linewidth=0.5))
    ax.add_patch(patches.Rectangle((19800000, 3400000), 500000, 200000, facecolor = "blue", edgecolor = "blue",linewidth=0.5,alpha=0.3))

    ax.add_patch(patches.Rectangle((19700000, 3000000), 900000, 200000, facecolor = "none", edgecolor = "black",linewidth=0.5))
    # we need generate a series site to show the exon
    distance = 0
    length = 0
    position = 19700000
    for j in range(0, 50):
        distance = random.randint(3,8)*7000
        length = random.randint(2,5)*4000
        position += distance
        ax.add_patch(patches.Rectangle((position, 3000000), length, 200000, facecolor = "green", edgecolor = "green",linewidth=0.5,alpha=0.5))
        position += length
        if position >= 20590000:
            break
    #ax.add_patch(patches.Rectangle((19800000, 3200000), 200000, 200000, facecolor = "green", edgecolor = "green",linewidth=0.5,alpha=0.1))
    ax.add_patch(patches.Rectangle((19700000, 3000000), 900000, 200000, facecolor = "none", edgecolor = "black",linewidth=0.5))

    ax.add_patch(patches.Rectangle((19700000, 2600000), 900000, 200000, facecolor = "none", edgecolor = "black",linewidth=0.5))
    ax.add_patch(patches.Rectangle((19800000, 2600000), 50000, 200000, facecolor = "black", edgecolor = "black",linewidth=0.5,alpha=0.6))

    ax.add_patch(patches.Rectangle((19700000, 2000000), 900000, 200000, facecolor = "none", edgecolor = "black",linewidth=0.5))
    distance = 0
    length = 0
    position = 19700000
    for j in range(0, 50):
        distance = random.randint(2,5)*9000
        length = random.randint(2,5)*5000
        position += distance
        ax.add_patch(patches.Rectangle((position, 2000000), length, 200000, facecolor = "green", edgecolor = "green",linewidth=0.5,alpha=0.5))
        position += length
        if position >= 20590000:
            break
    ax.add_patch(patches.Rectangle((19700000, 2000000), 900000, 200000, facecolor = "none", edgecolor = "black",linewidth=0.5))

    ax.add_patch(patches.Rectangle((19700000, 1600000), 900000, 200000, facecolor = "none", edgecolor = "black",linewidth=0.5))
    distance = 0
    length = 0
    position = 19700000
    for j in range(0, 50):
        distance = random.randint(3,8)*7000
        length = random.randint(2,5)*4000
        position += distance
        ax.add_patch(patches.Rectangle((position, 1600000), length, 200000, facecolor = "red", edgecolor = "red",linewidth=0.5,alpha=0.5))
        position += length
        if position >= 20590000:
            break
    ax.add_patch(patches.Rectangle((19700000, 1600000), 900000, 200000, facecolor = "none", edgecolor = "black",linewidth=0.5))

    #ax.add_patch(patches.Rectangle((19800000, 2000000), 400000, 200000, facecolor = "red", edgecolor = "red",linewidth=0.5,alpha=0.1))

    ax.text(20800000, 7500000, str("Exon"), ha='left', va='center', fontsize=10)
    ax.text(20800000, 6900000, str("UTR"), ha='left', va='center', fontsize=10)

    ax.text(20800000, 6400000, str("Promoter"), ha='left', va='center', fontsize=10)
    ax.text(20800000, 6000000, str("Exon"), ha='left', va='center', fontsize=10)
    ax.text(20800000, 5600000, str("Intron"), ha='left', va='center', fontsize=10)
    ax.text(20800000, 5200000, str("UTR"), ha='left', va='center', fontsize=10)
    ax.text(20800000, 4800000, str("TE protein"), ha='left', va='center', fontsize=10)
    ax.text(20800000, 4400000, str("Transposon"), ha='left', va='center', fontsize=10)
    ax.text(20800000, 4000000, str("IGR"), ha='left', va='center', fontsize=10)

    ax.text(20800000, 3500000, str("mRNA"), ha='left', va='center', fontsize=10)
    ax.text(20800000, 3100000, str("Exons"), ha='left', va='center', fontsize=10)
    ax.text(20800000, 2700000, str("UTR"), ha='left', va='center', fontsize=10)

    ax.text(20800000, 2100000, str("TE proteins"), ha='left', va='center', fontsize=10)
    ax.text(20800000, 1700000, str("Transposons"), ha='left', va='center', fontsize=10)
    # add the box to be continued

    ax.add_patch(patches.Rectangle((19600000, 1550000), 4500000, 770000, facecolor = "none", edgecolor = "black",linewidth=0.5))
    ax.add_patch(patches.Rectangle((23000000, 1600000), 20000, 700000, facecolor = "black", edgecolor = "black",linewidth=0.5))
    ax.text(23500000, 1930000, str("TEs"), ha='center', va='center', fontsize=10,rotation=90)
    ax.add_patch(patches.Rectangle((19600000, 2530000), 4500000, 1170000, facecolor = "none", edgecolor = "black",linewidth=0.5))
    ax.add_patch(patches.Rectangle((23000000, 2580000), 20000, 1070000, facecolor = "black", edgecolor = "black",linewidth=0.5))
    ax.text(23500000, 3020000, str("Gene"), ha='center', va='center', fontsize=10,rotation=90)
    ax.add_patch(patches.Rectangle((19600000, 3830000), 4500000, 4270000, facecolor = "none", edgecolor = "black",linewidth=0.5))
    ax.add_patch(patches.Rectangle((23000000, 3880000), 20000, 4170000, facecolor = "black", edgecolor = "black",linewidth=0.5))
    ax.text(23500000, 5900000, str("Somatic mutations"), ha='center', va='center', fontsize=10, rotation=90)
    #ax.add_patch(patches.Rectangle((19600000, 6800000), 4500000, 1400000, facecolor = "none", edgecolor = "black",linewidth=0.5))

    #ax.add_patch(patches.Rectangle((17000000, 7500000), 2000000, 200000, facecolor = "none", edgecolor = "black",linewidth=0.5))
    #ax.add_patch(patches.Rectangle((17000000, 7100000), 2000000, 200000, facecolor = "none", edgecolor = "black",linewidth=0.5))
    #ax.add_patch(patches.Rectangle((17000000, 6700000), 2000000, 200000, facecolor = "none", edgecolor = "black",linewidth=0.5))


    # bar plot the promoter
    ax.add_patch(patches.Rectangle((115000000, 9000000), 7000000, 6000000, facecolor = "white", edgecolor = "none",linewidth=0.5))
    ax.add_patch(patches.Rectangle((13000000, 9500000), 400000, promoter*2500, facecolor = "yellow", edgecolor = "black",linewidth=0.5,alpha=0.2))
    ax.add_patch(patches.Rectangle((13000000, 9500000), 400000, promoter_gene*2500, facecolor = "green", edgecolor = "none",linewidth=0.5,alpha=0.2))
    ax.add_patch(patches.Rectangle((13000000, 9500000), 400000, promoter_gene_transposon*2500, facecolor = "red", edgecolor = "none",linewidth=0.5,alpha=0.2))
    ax.add_patch(patches.Rectangle((13000000, 9500000+promoter_gene*2500), 400000, promoter_transposon_out_gene*2500, facecolor = "red", edgecolor = "none",linewidth=0.5,alpha=0.2))
                
    ax.add_patch(patches.Rectangle((14000000, 9500000), 400000, gene*2500, facecolor = "green", edgecolor = "black",linewidth=0.5,alpha=0.2))
    ax.add_patch(patches.Rectangle((14000000, 9500000), 400000, promoter_gene*2500, facecolor = "yellow", edgecolor = "none",linewidth=0.5,alpha=0.2))
    ax.add_patch(patches.Rectangle((14000000, 9500000), 400000, promoter_gene_transposon*2500, facecolor = "red", edgecolor = "none",linewidth=0.5,alpha=0.2))
    ax.add_patch(patches.Rectangle((14000000, 9500000+promoter_gene*2500), 400000, gene_transposon_out_promoter*2500, facecolor = "red", edgecolor = "none",linewidth=0.5,alpha=0.2))
    
    ax.add_patch(patches.Rectangle((15000000, 9500000), 400000, promoter_gene_transposon*2500, facecolor = "green", edgecolor = "none",linewidth=0.5,alpha=0.2))
    ax.add_patch(patches.Rectangle((15000000, 9500000), 400000, promoter_transposon*2500, facecolor = "yellow", edgecolor = "none",linewidth=0.5,alpha=0.2))
    ax.add_patch(patches.Rectangle((15000000, 9500000), 400000, (transposon+TEprotein_out)*2500, facecolor = "red", edgecolor = "black",linewidth=0.5,alpha=0.2))
    ax.add_patch(patches.Rectangle((15000000, 9500000+promoter_transposon*2500), 400000, gene_transposon_out_promoter*2500, facecolor = "green", edgecolor = "none",linewidth=0.5,alpha=0.2))
    ax.add_patch(patches.Rectangle((16000000, 9500000), 400000, IGR*2500, facecolor = "grey", edgecolor = "black",linewidth=0.1))

    # add the text
    ax.text(13400000, 9300000, "Promoter", fontsize=8, ha='center',va='center')
    ax.text(13200000, 9500000 +300000, str(format(promoter/mutation_site_num*100,'.0f'))+'%', ha='center', va='center', fontsize=8, rotation=90)
    ax.text(13200000, 9500000+(promoter + 50)*2500, str(promoter), fontsize=10, ha='center',va='center')
    ax.text(14400000, 9300000, "Gene", fontsize=10, ha='center',va='center')
    ax.text(14200000, 9500000+ 300000, str(format(gene/mutation_site_num*100,'.0f'))+'%', ha='center', va='center', fontsize=8, rotation=90)
    ax.text(14200000, 9500000+(gene + 50)*2500, str(gene), fontsize=10, ha='center',va='center')
    ax.text(15400000, 9300000, "TEs", fontsize=10, ha='center',va='center')
    ax.text(15200000, 9500000+ 300000, str(format((transposon + TEprotein_out)/mutation_site_num*100,'.0f'))+'%', ha='center', va='center', fontsize=8, rotation=90)
    ax.text(15200000, 9500000+(transposon + TEprotein_out + 50)*2500, str(transposon+TEprotein_out), fontsize=10, ha='center',va='center')
    ax.text(16400000, 9300000, "IGR", fontsize=10, ha='center',va='center')
    ax.text(16200000, 9500000+ 300000, str(format(IGR/mutation_site_num*100,'.0f'))+'%', ha='center', va='center', fontsize=8, rotation=90)
    ax.text(16200000, 9500000+(IGR + 50)*2500, str(IGR), fontsize=10, ha='center',va='center')

    # add the axies marks
    ax.add_patch(patches.Rectangle((12800000, 9500000), 20000, 2100*2500, facecolor = "black", edgecolor = "black",linewidth=0.5))
    for i in range(0, 2100, 500):
        ax.add_patch(patches.Rectangle((12700000, i*2500+9500000), 100000, 5*2500, facecolor = "black", edgecolor = "black",linewidth=0.5))
        ax.text(12670000, i*2500+9500000, str(i), fontsize=8, ha='right',va='center')
    ax.text(12000000, 1000*2500 + 9500000, "Number of mutation sites", fontsize=10, ha='center',va='center',rotation=90)
    
    
    genome_size = 216955151
    gene = 99495975
    promoter = 38963587
    TEs = 73155645
    gene_TE = 12776936
    promoter_TE = 5524050
    promoter_gene = 10782633
    promoter_gene_TE = 350126
    IGRS = genome_size - gene - promoter - TEs + gene_TE + promoter_TE + promoter_gene - promoter_gene_TE
    ax.add_patch(patches.Rectangle((13500000, 9500000), 400000, mutation_site_num*2500 * promoter/genome_size, facecolor = "yellow", edgecolor = "black",linewidth=0.5,alpha=0.2))
    ax.text(13700000, 9500000+300000, str(format(promoter/genome_size*100, '.0f'))+'%', ha='center', va='center', fontsize=8, rotation=90)
    ax.add_patch(patches.Rectangle((14500000, 9500000), 400000, mutation_site_num*2500 * gene/genome_size, facecolor = "green", edgecolor = "black",linewidth=0.5,alpha=0.2))
    ax.text(14700000, 9500000+300000, str(format(gene/genome_size*100, '.0f'))+'%', ha='center', va='center', fontsize=8, rotation=90)
    ax.add_patch(patches.Rectangle((15500000, 9500000), 400000, mutation_site_num*2500 * TEs/genome_size, facecolor = "red", edgecolor = "black",linewidth=0.5,alpha=0.2))
    ax.text(15700000, 9500000+300000, str(format(TEs/genome_size*100, '.0f'))+'%', ha='center', va='center', fontsize=8, rotation=90)
    ax.add_patch(patches.Rectangle((16500000, 9500000), 400000, mutation_site_num*2500 * IGRS/genome_size, facecolor = "grey", edgecolor = "black",linewidth=0.5,alpha=0.2))
    ax.text(16700000, 9500000+300000, str(format(IGRS/genome_size*100, '.0f'))+'%', ha='center', va='center', fontsize=8, rotation=90)
    
    ax.text(10000000, 14500000, "a", fontsize=24, ha='center',va='center')
    ax.text(13300000, 14500000, "b", fontsize=24, ha='center',va='center')
    ax.text(17500000, 14500000, "c", fontsize=24, ha='center',va='center')
    

    # plot the mutation

    plt.savefig("path/met/somatic_mutation_location_sample_" + str(sample) + ".pdf", dpi=600, format='pdf',bbox_inches='tight')
    #close the file
    plt.close()
