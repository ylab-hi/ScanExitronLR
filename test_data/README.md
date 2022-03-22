# Aradopsis Test Data

In this example, you will run ScanExitronLR on two small bam files from _Aradopsis thaliana_ and find a specific transcript enriched in exitrons. They are two technical replicates from the floral bud ([NCBI accession PRJNA605023](https://www.ncbi.nlm.nih.gov/bioproject?LinkName=biosample_bioproject&from_uid=14008366)). 

Please download the example bam files [here](https://drive.google.com/drive/folders/1JHY2dqf6O5QgOAJpALhBh2iAquf3EY7O?usp=sharing) or from this github page.

## Running ScanExitron

First, it is necessary to download a specifically formated version of TAIR10 reference genome and annotation, [which can be found here.](https://drive.google.com/drive/folders/1FNZ5HRJOvGeiMxMObXBPgTGC2E0l3yeE?usp=sharing) This includes a sorted GFF file and a GFF database (with the '.db' extension).  ScanExitron will create such a database if none is found, but downloading it here will save time. We require that the GFF annotation file is sorted for speed and memory efficiency. This example will assume you have downloaded these files and have them in the current working directory folder.

**STEP 1:** Enter the directory you have downloaded the BAM and TAIR files and run ScanExitronLR with the following command:

    selr extract -i ara_example_1.bam -g TAIR10.fas -r TAIR10_GFF3_sorted.gff.gz -sa -arabidopsis

The `-sa` flag tells ScanExitronLR to save the abundances for downstream analysis.  The `-aradopsis` runs ScanExitronLR in *Aradopsis* mode. If the flag is not set and *Aradopsis* alignments are used, you will get a `contig not found` error.

Run the same command for the second example:

    selr extract -i ara_example_2.bam -g TAIR10.fas -r TAIR10_GFF3_sorted.gff.gz -sa -arabidopsis

**STEP 2:** You will now have two exitron files, `ara_example_1.exitron` and  `ara_example_2.exitron`. In addition, because we saved the abundances, you should also have the files

    ara_example_1.isoform.exitrons
    ara_example_1.isoform.normals
    ara_example_2.isoform.exitrons
    ara_example_2.isoform.normals

In order to run LIQA's differential analysis command we need to save the path to these files in two lists:

    ls ara_example_*.isoform.exitrons > exitron_list
    ls ara_example_*.isoform.normals > normal_list

We can then run LIQA:

    liqa -task diff -condition_1 exitron_list -condition_2 normal_list -out exitron_das

Inspect the output: 

    cat exitron_das
    
    AT1G80245	0.647158117957978
    AT2G26030	0.000188791935434085

The second column is a p-value. We see that isoforms in the gene `AT2G26030`  are differentially expressed in exitron reads compared to normal reads. To see which transcript specifically is enriched with exitrons, we can use a simple `grep`: 

    grep "AT2G26030" *.isoform.*
    
    ara_example_1.isoform.exitrons:AT2G26030	AT2G26030.2	2.0	1.0	1.0
    ara_example_1.isoform.exitrons:AT2G26030	AT2G26030.3	0.0	0.0	1.0
    ara_example_1.isoform.exitrons:AT2G26030	AT2G26030.1	0.0	0.0	1.0
    ara_example_1.isoform.normals:AT2G26030	AT2G26030.2	2.200010756671194	0.20000097787919946	0.9090909090909091
    ara_example_1.isoform.normals:AT2G26030	AT2G26030.3	1.100015669145584	0.10000142446778036	0.9090909090909091
    ara_example_1.isoform.normals:AT2G26030	AT2G26030.1	7.699973574183222	0.6999975976530202	0.9090909090909091
    ara_example_2.isoform.exitrons:AT2G26030	AT2G26030.2	2.9994203350635424	0.999806778354514	0.3333333333333333
    ara_example_2.isoform.exitrons:AT2G26030	AT2G26030.3	0.000297190852227	9.906361740900907e-05	0.3333333333333333
    ara_example_2.isoform.exitrons:AT2G26030	AT2G26030.1	0.0002824740842305	9.415802807685568e-05	0.3333333333333333
    ara_example_2.isoform.normals:AT2G26030	AT2G26030.2	1.5715240176549718	0.14286581978681562	0.6363636363636364
    ara_example_2.isoform.normals:AT2G26030	AT2G26030.3	0.0001464890324780173	1.3317184770728845e-05	0.6363636363636364
    ara_example_2.isoform.normals:AT2G26030	AT2G26030.1	9.42832949331255	0.8571208630284136	0.6363636363636364

We see that exitrons occur predominantly in the AT2G26030.2 transcript in both replicates, while AT2G26030.1 is the predominant transcript in reads without exitrons--the fourth column reports relative abundance. 

## Annotation

	selr annotate -i ara_example_1.exitron -g TAIR10.fas -r TAIR10_GFF3_sorted.gff.gz -b ara_example_1.bam -fasta -arabidopsis

	selr annotate -i ara_example_2.exitron -g TAIR10.fas -r TAIR10_GFF3_sorted.gff.gz -b ara_example_2.bam -fasta -arabidopsis

