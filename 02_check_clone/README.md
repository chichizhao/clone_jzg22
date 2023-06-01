#  check the clone population
## 1 check the clone population by GATK haplotypecaller
### material
    reference genome: typha.fa
    raw NGS reads: see data.list in data directory
### build the index of reference genome
    BWA: https://github.com/lh3/bwa
    Samtools: https://github.com/samtools/samtools
    GATK4: https://github.com/broadgsa/gatk
    bwa index -p typha typha.fa
    samtools faidx typha.fa
    java -Djava.io.tmpdir=/tmp -Xmx128g -jar path/gatk-package-4.2.4.1-local.jar CreateSequenceDictionary -R typha.fa -O typha.dict
### filter the raw reads
    fastp: https://github.com/OpenGene/fastp
    for sample in $(cat sample_list.txt);
    do 
    fastp -i ${sample}_1.fq.gz -I ${sample}_2.fq.gz -o ${sample}_1.fq.gz -O ${sample}_2.fq.gz -h ${sample}.html -j ${sample}.json -w 8 -q 20 -u 20 -n 5 -l 50 -y 20 -x --umi --umi_loc read1 --umi_len 10;
    done

### mapping the reads to reference genome
    BWA: https://github.com/lh3/bwa
    Samtools: https://github.com/samtools/samtools
    for sample in $(cat sample_list.txt);
    do bwa mem -t 8 -R "@RG\tID:${sample}\tSM:${sample}\tPL:illumina" typha ${sample}_1.fq.gz ${sample}_2.fq.gz | samtools view -bS | samtools sort -o ${sample}.bam;
    samtools index ${sample}.bam;
    java -Djava.io.tmpdir=/tmp -Xmx128g -jar path/gatk-package-4.2.4.1-local.jar MarkDuplicates -I ${sample}.bam -O ${sample}_dedup.bam -M ${sample}_dedup.metrics.txt;
    samtools index ${sample}_dedup.bam;
    done
### call variants
    GATK4:  https://github.com/broadgsa/gatk
    for sample in $(cat sample_list.txt);
    do java -Djava.io.tmpdir=/tmp -Xmx128g -jar path/gatk-package-4.2.4.1-local.jar HaplotypeCaller --native-pair-hmm-threads 16 -R typha.fa -I ${sample}_dedup.bam -ERC GVCF -O ${sample}.g.vcf.gz;
    done
### join genotyping
    GATK4:  https://github.com/broadgsa/gatk
    for file in $(cat list);do echo -en --variant ${file}.g.vcf' '; done > mergevcf
    echo  >> mergevcf
    cat mergevcf | while read line;do
    java -Djava.io.tmpdir=/tmp -Xmx128g -jar path/gatk-package-4.2.4.1-local.jar CombineGVCFs -R typha.fa $line -O merge.g.vcf.gz
    java -Djava.io.tmpdir=/tmp -Xmx128g -jar path/gatk-package-4.2.4.1-local.jar GenotypeGVCFs -R typha.fa -V merge.g.vcf.gz -O typha_raw.vcf.gz
    done
#### here we got the raw variants vcf file
    typha_raw.vcf.gz


## 2 confirm the variants by Sanger sequencing
### prepare files
    genome: typha.fa
    ab1 files: see data.list in data directory
|Chromosome|primer_start|primer_end|
|:---:|:---:|:---:|
|Chr03|6328725|6329308|
|Chr04|13002465|13003003|
|Chr05|13386809|13387563|
|Chr07|3574237|3574712|
|Chr09|5593369|5593746|
|Chr11|677683|678158|
|Chr12|3432665|3433448|
|Chr13|9904283|9904920|
|Chr14_1|1654473|1655303|
|Chr14_2|4407470|4408007|
|Chr15|2720291|2721212|
    seqkit: https://github.com/shenwei356/seqkit
    mafft: https://github.com/GSLBiotech/mafft
### chr03
    seqkit subseq --chr Chr03 -r 6328705:6329308 path/typha.fa| awk -F [' '] '{print $1}' | seqkit seq -p -r > chr3.fa
    cat $(ls jzg_Chr03_*.fasta) > chr3.fa
    mafft --threads 40 --auto chr3.fa > chr3.aln
    # manually check the alignment
    python plot_alignment.py --alignment chr3.aln --sanger jzg_Chr03_*.ab1_path --output chr3
### chr04
    seqkit subseq --chr Chr04 -r 13002465:13003003 path/typha.fa| awk -F [' '] '{print $1}' > chr4.fa
    cat $(ls jzg_Chr04_*.fasta) > chr4.fa
    mafft --threads 40 --auto chr4.fa > chr4.aln
    # manually check the alignment
    python plot_alignment.py --alignment chr4.aln --sanger jzg_Chr04_*.ab1_path --output chr4
### chr05
    seqkit subseq --chr Chr05 -r 13386809:13387563 path/typha.fa| awk -F [' '] '{print $1}'| seqkit seq -p -r > chr5.fa
    cat $(ls jzg_Chr05_*.fasta) > chr5.fa
    mafft --threads 40 --auto chr5.fa > chr5.aln
    # manually check the alignment
    python plot_alignment.py --alignment chr5.aln --sanger jzg_Chr05_*.ab1_path --output chr5
### chr07
    seqkit subseq --chr Chr07 -r 3574237:3574712 path/typha.fa| awk -F [' '] '{print $1}'> chr7.fa
    cat $(ls jzg_Chr07_*.fasta) > chr7.fa
    mafft --threads 40 --auto chr7.fa > chr7.aln
    # manually check the alignment
    python plot_alignment.py --alignment chr7.aln --sanger jzg_Chr07_*.ab1_path --output chr7
### chr09
    seqkit subseq --chr Chr09 -r 5593369:5593746 path/typha.fa| awk -F [' '] '{print $1}' | seqkit seq -p -r > chr9.fa
    cat $(ls jzg_Chr09_*.fasta) > chr9.fa
    mafft --threads 40 --auto chr9.fa > chr9.aln
    # manually check the alignment
    python plot_alignment.py --alignment chr9.aln --sanger jzg_Chr09_*.ab1_path --output chr9
### chr11
    seqkit subseq --chr Chr11 -r 677683:678158 path/typha.fa| awk -F [' '] '{print $1}' | seqkit seq -p -r > chr11.fa
    cat $(ls jzg_Chr11_*.fasta) > chr11.fa
    mafft --threads 40 --auto chr11.fa > chr11.aln
    # manually check the alignment
    python plot_alignment.py --alignment chr11.aln --sanger jzg_Chr11_*.ab1_path --output chr11
### chr12
    seqkit subseq --chr Chr12 -r 3432665:3433448 path/typha.fa| awk -F [' '] '{print $1}' > chr12.fa
    cat $(ls jzg_Chr12_*.fasta) > chr12.fa
    mafft --threads 40 --auto chr12.fa > chr12.aln
    # manually check the alignment
    python plot_alignment.py --alignment chr12.aln --sanger jzg_Chr12_*.ab1_path --output chr12
### chr13
    seqkit subseq --chr Chr13 -r 9904283:9904920 path/typha.fa| awk -F [' '] '{print $1}' | seqkit seq -p -r > chr13.fa
    cat $(ls jzg_Chr13_*.fasta) > chr13.fa
    mafft --threads 40 --auto chr13.fa > chr13.aln
    # manually check the alignment
    python plot_alignment.py --alignment chr13.aln --sanger jzg_Chr13_*.ab1_path --output chr13
### chr14_1
    seqkit subseq --chr Chr14 -r 1654473:1655303 path/typha.fa| awk -F [' '] '{print $1}' > chr14_1.fa
    cat $(ls jzg_Chr14_1_*.fasta) > chr14_1.fa
    mafft --threads 40 --auto chr14_1.fa > chr14_1.aln
    # manually check the alignment
    python plot_alignment.py --alignment chr14_1.aln --sanger jzg_Chr14_1_*.ab1_path --output chr14_1
### chr14_2
    seqkit subseq --chr Chr14 -r 4407470:4408007 path/typha.fa| awk -F [' '] '{print $1}' > chr14_2.fa
    cat $(ls jzg_Chr14_2_*.fasta) > chr14_2.fa
    mafft --threads 40 --auto chr14_2.fa > chr14_2.aln
    # manually check the alignment
    python plot_alignment.py --alignment chr14_2.aln --sanger jzg_Chr14_2_*.ab1_path --output chr14_2
### chr15
    seqkit subseq --chr Chr15 -r 2720291:2721212 path/typha.fa| awk -F [' '] '{print $1}' > chr15.fa
    cat $(ls jzg_Chr15_*.fasta) > chr15.fa
    mafft --threads 40 --auto chr15.fa > chr15.aln
    # manually check the alignment
    python plot_alignment.py --alignment chr15.aln --sanger jzg_Chr15_*.ab1_path --output chr15
#### after the plot we got the variants rawdata of each primer
    chr03.txt
    chr04.txt
    chr05.txt
    chr07.txt
    chr09.txt
    chr11.txt
    chr12.txt
    chr13.txt
    chr14_1.txt
    chr14_2.txt
    chr15.txt
### plot the all 11 primers variants
    python plot_total_sanger.py --sanger txt_path

## 3 check the clone population by somatic mutation calling tools
### 3.1 mutect2 tumor-only mode
    GATK4:https://github.com/broadgsa/gatk
    for sample in $(cat sample_list.txt);
    java -Djava.io.tmpdir=/tmp -Xmx128g -jar path/gatk-package-4.2.4.1-local.jar Mutect2 -R typha.fa -I ${sample}_dedup.bam -tumor ${sample} -O ${sample}_mutect2.vcf.gz
    done
### 3.2 strelka2 tumor-only mode
    strelka2: https://github.com/Illumina/strelka
    for sample in $(cat sample_list.txt);
    do 
    configureStrelkaSomaticWorkflow.py --bam ${sample}_dedup.bam --referenceFasta typha.fa --runDir ${sample}_strelka2
    done
    for sample in $(cat sample_list.txt);
    do
    cd ${sample}_strelka2
    ./runWorkflow.py -m local -j 8 --quiet
    cd ..
    mv ${sample}_strelka2/results/variants/somatic.snvs.vcf.gz ${sample}_strelka2.vcf.gz
    done
### 3.3 varscan2 tumor-only mode
    Varscan2: https://github.com/Jeltje/varscan2
    Samtools: https://github.com/samtools/samtools
    for sample in $(cat sample_list.txt);
    do
    samtools mpileup -f typha.fa ${sample}_dedup.bam | varscan mpileup2snp --output-vcf 1 | varscan filter --output-vcf 1 > ${sample}_varscan2.vcf
    done
### 3.4 lofreq tumor-only mode
    lofreq: https://github.com/CSB5/lofreq
    for sample in $(cat sample_list.txt);
    do
    lofreq call-parallel --pp-threads 8 -f typha.fa -o ${sample}_lofreq.vcf ${sample}_dedup.bam
    done
### 3.5 summarize the somatic mutation calling results with SomaticSeq tumor-only mode
    SomaticSeq: https://github.com/bioinform/somaticseq
    for sample in $(cat sample_list.txt);
    do
    somaticseq_parallel.py --output-directory somatic_'${sample}' --genome-reference typha.fa single --bam-file '${sample}'_dedup.bam --mutect2-vcf '${sample}'_mutect2.vcf.gz --varscan-vcf '${sample}'_varscan2.vcf --lofreq-vcf '${sample}'_lofreq.vcf --strelka-vcf '${sample}'_strelka2.vcf.gz
    done

### 3.6 final vcfs merge
#### remove the REJECT variants
    bcftools: https://github.com/samtools/bcftools
    for file in $(cat list);do
    python remove_reject.py -i somatic_${file}/SomaticSeq.vcf -o somatic_${file}/SomaticSeq_remove_reject.vcf
    done
    for file in $(cat list);do 
    bgzip somatic_${file}/SomaticSeq.vcf && tabix -p vcf somatic_${file}/SomaticSeq.vcf.gz
    done
    bcftools merge -m none -O v -o typha_somatic.vcf.gz somatic_*.vcf.gz
#### here we got the final merge.vcf.gz
    typha_somatic.vcf.gz

## 4 scipts we used for plot the variants distribution

### For the germline variants call method, we use the following script to count and plot the figs.
### First we count the variants distribution of each samples, here we use the python scripts to achieve it.
    python variants_count.py -vcf vcf_path -o output_prefix
### Then we use the R scripts to plot the variants distribution of each samples.
    Rscript plot_variants_count.r
### Next we count the each type variants distribution in samples, here we use the python scripts to achieve it.
    python variants_type_count.py -vcf vcf_in_all_samples -o output_prefix
### After that we got the output_prefix_snps.csv and output_prefix_indels.csv, then we use the matlab scripts to plot the variants distribution of each type.
### we need load the snps and indels csv files in matlab, then run the following scripts.
    matlab plot_variants_type_count.m
### For the Sanger sequencing variants call method, we use the following script to count and plot the figs and it also used in the step 2.
    python plot_total_sanger.py --sanger txt_path
### similarly, for the solid somatic mutation calling method, we use the following script which same above to count and plot the figs.
    python variants_count.py -vcf vcf_solids_path -o output_prefix
    Rscript plot_variants_count.r
    python variants_type_count.py -vcf vcf_solids_in_all_samples -o output_prefix
    matlab plot_variants_type_count.m
### for the venne plot, we use the following script to plot the figs. we plot it manually.







