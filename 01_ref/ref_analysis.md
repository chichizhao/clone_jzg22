# 01 Flow of build reference genome
## prepared files
    jzg_nanopore.fq.gz
    jzg_NGS_1.fq.gz
    jzg_NGS_2.fq.gz
    jzg_Hic_1.fq.gz
    jzg_Hic_2.fq.gz
    typha.flRNA.fastq.gz
    jzg_leef1_1.fq.gz
    jzg_leef1_2.fq.gz
    jzg_leef2_1.fq.gz
    jzg_leef2_2.fq.gz
    jzg_shoot1_1.fq.gz
    jzg_shoot1_2.fq.gz
    jzg_shoot2_1.fq.gz
    jzg_shoot2_2.fq.gz
    jzg_shoot2_2.fq.gz
## estimatate the genome size
### Filter NGS reads
    Fastp: https://github.com/OpenGene/fastp
    fastp -i jzg_NGS_1.fq.gz -I jzg_NGS_2.fq.gz -o jzg_NGS_1.clean.fq.gz -O jzg_NGS_2.clean.fq.gz -h jzg_NGS.html -j jzg_NGS.json
### estimate the genome size 
    GenomeScope: https://github.com/schatzlab/genomescope
    jellyfish count -m 21 -s 1000000000 -t 40 -C jzg_NGS_1.clean.fq.gz jzg_NGS_2.clean.fq.gz -o jzg_NGS.jf
    jellyfish histo -t 40 jzg_NGS.jf > jzg_NGS.histo
    Rscript path/genomescope.R jzg_NGS.histo 21 150 jzg_NGS
## Assembly
### Filter nanopore reads
    NanoFilt: https://github.com/wdecoster/nanofilt
    NanoFilt -q 10 -l 500 -t 40 -o jzg_nanopore.clean.fq.gz jzg_nanopore.fq.gz
### Assembly
    NextDenovo: https://github.com/Nextomics/NextDenovo
    nextDenovo Nextdenovo.cfg
### polish the assembly
#### Polish the assembly both by NGS and nanopore reads twice.
    NextPolish: https://github.com/Nextomics/NextPolish
    nextPolish NextPolish.cfg
## estimate draft assembly
### estimate the draft assembly by BUSCO
    BUSCO: https://busco.ezlab.org/
    busco -i typha_draft.fa -o typha_draft -l path/embryophyta_odb10 -m genome -c 40
### estimate the draft's snp and indel by NGS reads
    BWA: https://github.com/lh3/bwa
    Samtools: https://github.com/samtools/samtools
    bcftools: https://github.com/samtools/bcftools
    bwa index typha_draft.fa
    bwa mem -t 40 typha_draft.fa jzg_NGS_1.clean.fq.gz jzg_NGS_2.clean.fq.gz | samtools view -bS - | samtools sort -@ 40 -o jzg_NGS.bam
    samtools index jzg_NGS.bam
    samtools mpileup -uf typha_draft.fa jzg_NGS.bam | bcftools call -mv -Oz -o jzg_NGS.vcf.gz
## align the Hi-C reads to the draft assembly
### filter the Hi-C reads
    Fastp: https://github.com/OpenGene/fastp
    fastp -i jzg_Hic_1.fq.gz -I jzg_Hic_2.fq.gz -o jzg_Hic_1.clean.fq.gz -O jzg_Hic_2.clean.fq.gz -h jzg_Hic.html -j jzg_Hic.json
### align the Hi-C reads to the draft assembly
    3d-dna: https://github.com/aidenlab/3d-dna
    juicer: https://github.com/aidenlab/juicer
    BWA: https://github.com/lh3/bwa
    bwa index typha_draft.fa
    python path/juicer/misc/generate_site_positions.py DpnII typha typha_draft.fa
    awk 'BEGIN{OFS="\t"}{print $1, $NF}' typha_DpnII.txt > typha_DpnII.bed
    juicer.sh -g tyha -s DnpII -z typha_draft.fa -D path/juicer/CPU/common/ -d path/juicer/CPU/ -y typha_DnpII.txt -p typha_DpnII.bed -t 8
    run-asm-pipeline.sh path/typha_draft.fa path/aligned/merged_nodups.txt
### After the 3d-dna, we get the draft assembly with chromosome level. and we check the final result by Juicebox.
    Juicebox:https://github.com/aidenlab/Juicebox
### rerun the 3d-dna with the chromosome level draft assembly
    run-asm-pipeline-review.sh -r path/typha_draft.fa path/aligned/merged_nodups.txt

# 02 Annotation
## prepare files
    typha_draft.fa
    
## Repeat annotation
    RepeatModeler: https://github.com/Dfam-consortium/RepeatModeler
    RepeatMasker: https://github.com/rmhubley/RepeatMasker
    BuildDatabase -name typha -engine ncbi typha_draft.fa
    RepeatModeler -database typha -engine ncbi -pa 40
    RepeatMasker -pa 40 -s -lib path/typha-families.fa typha_draft.fa
## Gene structure annotation
### prepare files
    typha_draft.fa.masked
### Ab initio gene prediction
    BRAKER: https://github.com/Gaius-Augustus/BRAKER
    GeneMark-ES: http://exon.gatech.edu/GeneMark/
    AUGUSTUS: https://github.com/Gaius-Augustus/Augustus
    braker.pl --species=typha --genome=typha_draft.fa.masked --prot_seq=proteins.fasta --softmasking --gff3 --cores=48 --workingdir=ab_initio --min_contig=1000
#### after the BRAKER, we get the ab initio gene annotation file 
    mv augustus.hints.gff3 typha_ab_initio.gff3
    
### Homology-based gene prediction
    GenomeThreader: https://github.com/genometools/genomethreader
    Homology-based species: Oryza sativa, Arabidopsis thaliana, Zea mays
#### prepare files
    Oryza sativa protein: Oryza_sativa.IRGSP-1.0.pep.all.fa
    Arabidopsis thaliana protein: Arabidopsis_thaliana.TAIR10.pep.all.fa
    Zea mays protein: Zea_mays.B73_RefGen_v4.pep.all.fa
    typha_draft.fa.masked
    
    cat Oryza_sativa.IRGSP-1.0.pep.all.fa Arabidopsis_thaliana.TAIR10.pep.all.fa Zea_mays.B73_RefGen_v4.pep.all.fa > homology_species.pep.fa
    makeblastdb -in homology_species.pep.fa -dbtype prot -out homology_species.pep.fa
    GenomeThreader -gff3out -genomic typha_draft.fa.masked -protein homology_species.pep.fa -skipalignmentout -skipalignme -o typha_homology.gff3
#### after the GenomeThreader, we get the homology-based gene annotation file
    typha_homology.gff3

### RNA-seq based gene prediction
    StringTie: https://github.com/gpertea/stringtie
    HISAT2: https://github.com/DaehwanKimLab/hisat2
    Samtools: https://github.com/samtools/samtools
    TransDecoder: https://github.com/TransDecoder/TransDecoder
#### prepare files
    Full-length RNA-seq reads: typha.flRNA.fastq.gz
    RNA-seq reads: 
        jzg_leef1_1.fq.gz
        jzg_leef1_2.fq.gz
        jzg_leef2_1.fq.gz
        jzg_leef2_2.fq.gz
        jzg_shoot1_1.fq.gz
        jzg_shoot1_2.fq.gz
        jzg_shoot2_1.fq.gz
        jzg_shoot2_2.fq.gz
        jzg_shoot2_2.fq.gz
#### align the RNA-seq reads to the draft assembly
    hisat2-build typha_draft.fa.masked typha_draft
    hisat2 -p 40 -x typha_draft -1 jzg_leef1_1.fq.gz -2 jzg_leef1_2.fq.gz | samtools sort -@ 40 -o jzg_leef1.bam
    hisat2 -p 40 -x typha_draft -1 jzg_leef2_1.fq.gz -2 jzg_leef2_2.fq.gz | samtools sort -@ 40 -o jzg_leef2.bam
    hisat2 -p 40 -x typha_draft -1 jzg_shoot1_1.fq.gz -2 jzg_shoot1_2.fq.gz| samtools sort -@ 40 -o jzg_shoot1.bam
    hisat2 -p 40 -x typha_draft -1 jzg_shoot2_1.fq.gz -2 jzg_shoot2_2.fq.gz | samtools sort -@ 40 -o jzg_shoot2.bam
    hisat2 -p 40 -x typha_draft -U typha.flRNA.fastq.gz | samtools sort -@ 40 -o typha.flRNA.bam
    samtools merge -@ 40 -f bam -O bam jzg_leef1.bam jzg_leef2.bam jzg_shoot1.bam jzg_shoot2.bam > typha_rna.bam
    stringtie -p 40 --mix -o typha_rna.gtf typha_rna.bam typha.flRNA.bam
    gtf_genome_to_cdna_fasta.pl typha_rna.gtf typha_draft.fa.masked > typha_rna.fa
    gtf_to_alignment_gff3.pl typha_rna.gtf > typha_rna.gff3
    TransDecoder.LongOrfs -t typha_rna.fa
#### After the TransDecoder, we get the RNA-seq based gene annotation file
    typha_rna.gff3
    
## Merge the three gene annotation files
#### we use the EVidenceModeler to merge the three gene annotation files
    EVidenceModeler: https://github.com/EVidenceModeler/EVidenceModeler
    
    weight.txt
        Ab_initio 2
        Homology 2
        RNA-seq 6
        
    prepare files
        typha_ab_initio.gff3
        typha_homology.gff3
        typha_rna.gff3
    partition_EVM_inputs.pl --genome typha_draft.fa.masked --gene_predictions typha_ab_initio.gff3 --protein_alignments typha_homology.gff3 --transcript_alignments typha_rna.gff3 --segmentSize 100000 --overlapSize 10000 --partition_listing partitions_list.out
    write_EVM_commands.pl --genome typha_draft.fa.masked --weights weight.txt --gene_predictions typha_ab_initio.gff3 --protein_alignments typha_homology.gff3 --transcript_alignments typha_rna.gff3 --output_file_name evm.out --partitions partitions_list.out >  commands.list
    bash commands.list
    recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out
    convert_EVM_outputs_to_GFF3.pl  --partitions partitions_list.out --output evm.out  --genome typha_draft.fa.masked
#### after the EVidenceModeler, we get the final gene annotation file
    evm.out.gff3
### update the gene annotation file with the PASA
    PASApipeline: https://github.com/PASApipeline/PASApipeline
#### prepare files
    typha_draft.fa.masked
    typha_rna.fa
    typha_rna.gff3
    evm.out.gff3
    
    gtf_to_cdna_fasta.pl typha_rna.gff3 typha_rna.fa > transcripts.fasta
    seqclean transcripts.fasta
    Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g typha_draft.fa.masked -t transcripts.fasta.clean -T -u transcripts.fasta.clean --ALIGNERS gmap --CPU 16
#### here we need run it twice.
    pasa_gff3_validator.pl pasa_assemblies.gff3 > pasa_assemblies.gff3.valid
    Load_Current_Gene_Annotations.dbi -c pasa_conf.txt -g typha_draft.fa.masked -P pasa_assemblies.gff3.valid
    
#### after the PASA, we get the final gene annotation file
    pasa_assemblies.gff3.valid
    mv pasa_assemblies.gff3.valid typha_final.gff3

## functional annotation
### Swiss-Prot, pfam, GO, KEGG, NR, InterPro
    annotate wit blast2go

# 03 species phylogenetic tree
## Species
| Assembly  | Species   | Source |
| :---: | :---: | :---: |
|GCF_004353265.1    |*Vutis ripapria* |https://www.ncbi.nlm.nih.gov/assembly/GCF_004353265.1/|
|GCA_902729315.2	|*Spirodela intermedia*   |https://www.ncbi.nlm.nih.gov/assembly/GCA_902729315.2/|
|GCF_000001605.2	|*Sorghun bicolor*    |https://www.ncbi.nlm.nih.gov/assembly/GCF_000003195.3/|
|GCF_001263595.1	|*Phalaenopsis equestris*   |https://www.ncbi.nlm.nih.gov/assembly/GCF_001263595.1/|
|GCF_001433935.1	|*Oryza sativa*   |https://www.ncbi.nlm.nih.gov/assembly/GCF_001433935.1/|
|GCF_000313855.2	|*Musa acuminata*  |https://www.ncbi.nlm.nih.gov/assembly/GCF_000313855.2/|
|GCF_000442705.1	|*Elaeis guineensis*  |https://www.ncbi.nlm.nih.gov/assembly/GCF_000442705.1/|
|GCF_001876935.1	|*Asparagus officinalis*  |https://www.ncbi.nlm.nih.gov/assembly/GCF_001876935.1/|
|GCF_001540865.1	|*Ananas comosus* |https://www.ncbi.nlm.nih.gov/assembly/GCF_001540865.1/|
|GCF_000130695.1	|*Amborella trichopoda*   |https://www.ncbi.nlm.nih.gov/assembly/GCF_000130695.1/|

## build the tree with OrthoFinder
    OrthoFinder:https://github.com/davidemms/OrthoFinder
    mafft: https://mafft.cbrc.jp/alignment/software/
    mcmctree: https://github.com/abacus-gene/paml
#### prepare files
    typha.pep
    Vutis.pep
    Spirodela.pep
    Sorghun.pep
    Phalaenopsis.pep
    Oryza.pep
    Musa.pep
    Elaeis.pep
    Asparagus.pep
    Ananas.pep
    Amborella.pep
    orthofinder -f ./ -t 40 -S diamond
#### tree file is in the OrthoFinder/Results_data/WorkingDirectory/Species_Tree/ directory
    cp OrthoFinder/Results_data/Singe_Copy_Orthologue_Sequences/* ./data/
    for file in $(ls data/*.fa);do
        mafft --auto $file > $file.mafft.fas
    done
    python3 merge_species.py data merged.fas
    mafft --auto merged.fas > merged.fas.mafft.fas
    trimal -in merged.fas.mafft.fas -out merged.fas.mafft.fas.trimal -automated1
    mcmctree mcmctree.ctl
#### after the mcmctree, we get the species tree

# 04 WGD analysis
    MCScanX:   https://github.com/wyp1125/MCScanX
    ParaAT: https://github.com/wonaya/ParaAT
    Kaks_Calculator: https://sourceforge.net/projects/kakscalculator2/
## prepare files
    typha.pep
    typha.gff3
    typha_final.gff3
    typha_mcscanx.gff
    axt_cal.py 
    
    makeblastdb -in typha.pep -dbtype prot -out typha.pep
    blastp -query typha.pep -db typha.pep -outfmt 6 -evalue 1e-10 -num_threads 40 -out typha.pep
    cp typha_mcscanx.gff typha.gff
    MCScanX typha.gff
    cat typha.gff.collinearity | rep "Chr" | awk '{print $3"\t"$4}'> typha.homolog
    cat 40 > threads
    perl path/ParaAT.pl -h typha.homolog -n typha.pep -a typha.pep -p threads -m clustalw -f axt -o typha_out
    for i in `ls *.axt`;do KaKs_Calculator -i $i -o ${i}.kaks -m YN;done
    for i in `ls *.axt`;do python axt_cal.py $i ${i}.one-line;done
    for i in `ls *.kaks`;do awk 'NR>1{print $1"\t"$3"\t"$4"\t"$5}' $i >>all-kaks.txt;done
    sort all-kaks.txt|uniq >all-kaks.results
#### here we got the all-kaks.results file, by this method, we can get the WGD events in *Oryza sativa*, *Sorghun bicolor*, and *Ananas comosus*


# 05 plots and calculates
## plot the reads distrubution here we use Matlab file to plot it. 
    plot_reads.m
## HiC heatmap we plot it by python package plotHicGenome
    plotHicGenome: https://github.com/Atvar2/plotHicGenome
## Genome circos plot with python packages pycircos in plot_genome.py
    pycircos: https://github.com/ponnhide/pyCircos
    plot_genome.py
## phylogentic tree plot with itol
    itol: https://itol.embl.de
## *Ananas comosus* and *Typha latifolia* linkages is plot by python packages jcvi
    jcvi: https://github.com/tanghaibao/jcvi
## the curve fit of *Oryza sativa*, *Sorghun bicolor*,  *Ananas comosus* and *Typha latifolia* are use Matlab 2022b curve fitter tools. And the plot we plot it with kaks_plot.m
    kaks_plot.m




    
    
    
    

    

