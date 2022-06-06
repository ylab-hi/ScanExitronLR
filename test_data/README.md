# Aradopsis Test Data

In this example, you will run ScanExitronLR on two small bam files from _Aradopsis thaliana_ and find a specific transcript enriched in exitrons. They are two technical replicates from the floral bud ([NCBI accession PRJNA605023](https://www.ncbi.nlm.nih.gov/bioproject?LinkName=biosample_bioproject&from_uid=14008366)).

Please download the example bam files [here](https://drive.google.com/drive/folders/1JHY2dqf6O5QgOAJpALhBh2iAquf3EY7O?usp=sharing) or from this github page.

There is also a human example below.

## Running ScanExitronLR

First, it is necessary to download a specifically formated version of TAIR10 reference genome and annotation, [which can be found here.](https://drive.google.com/drive/folders/1FNZ5HRJOvGeiMxMObXBPgTGC2E0l3yeE?usp=sharing) This includes a sorted GFF file and a GFF database (with the '.db' extension). ScanExitron will create such a database if none is found, but downloading it here will save time. We require that the GFF annotation file is sorted for speed and memory efficiency. This example will assume you have downloaded these files and have them in the current working directory folder.

**STEP 1:** Enter the directory you have downloaded the BAM and TAIR files and run ScanExitronLR with the following command:

    selr extract -i ara_example_1.bam -g TAIR10.fas -r TAIR10_sorted.gtf.gz -a 1 -p 0.05 -sa

The `-sa` flag tells ScanExitronLR to save the abundances for downstream analysis. 

Run the same command for the second example:

    selr extract -i ara_example_2.bam -g TAIR10.fas -r TAIR10_sorted.gtf.gz -a 1 -p 0.05 -sa

**STEP 2:** You will now have two exitron files, `ara_example_1.exitron` and `ara_example_2.exitron`. In addition, because we saved the abundances, you should also have the files

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

The second column is a p-value. We see that isoforms in the gene `AT2G26030` are differentially expressed in exitron reads compared to normal reads. To see which transcript specifically is enriched with exitrons, we can use a simple `grep`:

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

You can also annotate the exitrons with the following commands:

    selr annotate -i ara_example_1.exitron -g TAIR10.fas -r TAIR10_GFF3_sorted.gff.gz -b ara_example_1.bam -fasta -arabidopsis

    selr annotate -i ara_example_2.exitron -g TAIR10.fas -r TAIR10_GFF3_sorted.gff.gz -b ara_example_2.bam -fasta -arabidopsis

If you inspect the second example, you will see that the 90 length exitron in AT2G26030 is a truncation plus substitution, though no amino acid was changed because the same amino acid happened to be substituted (G->G).

    Chr2	11092582	11092673	AT2G26030d11092582-11092673	CDS	2	-	AT2G26030	AT2G26030	90	GT-AG	AT2G26030.2,0.9998	0.1176	17	0.5	36	truncated+substitution	G->G	.	.	.	.	.	,SRR11031291.368547,SRR11031291.947060
    

# Human Test Data

In this example, we provide a small number of direct RNA alignments from the A549 cell line from the [SG-NEx sequencing project](https://github.com/GoekeLab/sg-nex-data). You can download the test data [here](https://drive.google.com/drive/folders/1pujKF4pb03gVcEWkJZtO0AcX8_6RL6M9?usp=sharing). 

Run ScanExitronLR using your GTF annotation:

    selr extract -i a549_example.bam -g hg38.fa -r your_annotation.gtf.gz

The output, `a549_example.exitron` should have the following exitrons:

    chrom	start	end	name	region	ao	strand	gene_name	gene_id	length	splice_site	transcript_id	pso	dp	cluster_purity	reads
    chr1	89712935	89713026	LRRC8Cd89712935-89713026	CDS	3	+	LRRC8C	ENSG00000171488.15	90	GT-AG	ENST00000370454.9,1.0	0.1	30	1.0	,2c6887a4-d77c-4ac8-a8fb-c3bde4e56183,8722a180-6687-4662-bb9e-67ca71f99009,bfce6d0f-0cd0-416e-8d9e-be42b5c1ec98
    chr1	11923115	11923191	KIAA2013d11923115-11923191	CDS	6	-	KIAA2013	ENSG00000116685.16	75	GT-AG	ENST00000376572.8,1.0	0.0295	203	0.83	,8d85452b-84c7-41bb-9ac5-d542b36f1e51,035f42c1-23ba-42e5-af9b-5c70e17e979e,e635263a-c821-43a8-baa1-db8a1cc14866,1acc9269-babf-4015-9b35-001e9d135c87,e81bd395-1047-4156-a8f6-a6ee77e00b9d,6c3b2381-c739-4c41-ab61-db6ed2b38805
    chr2	181915414	181915922	ITPRID2d181915414-181915922	CDS	8	+	ITPRID2	ENSG00000138434.17	507	GT-AG	ENST00000320370.11,0.5176;ENST00000431877.7,0.4824	0.06741573033707865	118	1.0	,acf1976a-9aa8-4b77-b110-ebf2860baf0b,5fe479d8-69fc-4b04-bf3e-eeeaec3d8b18,7694f254-b279-4246-a013-293960cb1e92,3b25ba1a-251c-4d87-afec-b9909189dcdf,REALIGNED_ac325add-9bc6-4483-a28d-a3757282aeb2,REALIGNED_d7326ea2-cc3d-4aee-bda8-7ef12f82ef11,REALIGNED_7dce71d4-c76d-458b-b4a8-1a57c8d4ddd4,REALIGNED_0b8405bb-0332-4fce-aa12-cd4cea46bd67
    chr20	47635980	47636314	NCOA3d47635980-47636314	CDS	7	+	NCOA3	ENSG00000124151.19	333	GT-AG	ENST00000372004.7,1.0	0.06231454005934718	112	1.0	,8e61ccfc-d0fd-43be-a059-095eb176e0ef,909821ec-91a8-4180-a683-db791a1b68b6,786a06fc-55e4-4458-8f3e-9dc047da13cf,098f5359-724f-4ecf-b312-8e1f6fedd382,155c1da6-7135-46d4-b928-1dcf2df9fb9c,REALIGNED_61fc108d-feba-449b-90cf-1496dd25e376,REALIGNED_82e8eb05-729c-46c0-83e8-e5d86b4a0839
    chr20	47636178	47636630	NCOA3d47636178-47636630	CDS	2	+	NCOA3	ENSG00000124151.19	451	GT-AG	ENST00000372004.7,1.0	0.0172	116	0.5	,ed4bcc6d-0d41-4dfb-a95f-d8b3896caea9,9691ee69-db9d-4328-9f8d-db66a3220971
