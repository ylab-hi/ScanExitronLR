# ScanExitronLR


A computational workflow for exitron splicing identification in short-read RNAseq data
<img align="right" width="500" src="https://github.com/ylab-hi/ScanExitronLR/blob/main/splice_type.png#gh-light-mode-only">
<img align="right" width="500" src="https://github.com/ylab-hi/ScanExitronLR/blob/main/splice_type_dark.png#gh-dark-mode-only">

# Installation
The recommended way to install `ScanExitronLR` is using [pip](https://pip.pypa.io/en/stable/):

```bash
pip install selr
```
This will pull and install the latest stable release from [PyPi](https://pypi.org/). ScanExitronLR requires Python 3.7+. Thus you need to make sure that the `pip` is for python3 using e.g. `which pip` or using: 
```bash 
pip3 install selr
```

To test your installation, run:

```bash
selr
```

You should see the version number, usage instructions and commands.

# Usage
ScanExitronLR has two modes, `extract` and `annotate`. Use `extract` when calling exitrons in an alignment and `annotate` when annotating exitrons called using `extract`.



## Extract
`extract` requires three inputs: (1) a BAM alignment file of long-reads containing the ts:A flag (provided by default by Minimap2), (2) a reference genome and (3) a sorted and bgzip'd gene annotation file. Currently only gtf files, e.g. Gencode v38, and the TAIR _Arabidopsis thaliana_ GFF3 annotation file are supported. You can download the TAIR10 annotation and reference genome [here](https://drive.google.com/drive/folders/1FNZ5HRJOvGeiMxMObXBPgTGC2E0l3yeE?usp=sharing).

To sort your gtf file, use the command:

```
awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k4,4n -k5,5n"}' in.gtf > out_sorted.gtf
```

To bgzip your gene annotation file, use:

```
bgzip in.gtf
```

`bgzip` is part of the [htslib](http://www.htslib.org/), which you most likely already have installed if you care about BAM files. Otherwise, you can get it [here](http://www.htslib.org/). It is important to note that if you download the [latest GENCODE release](https://www.gencodegenes.org/human/) it will be in the gzip form, not bgzip. You will need to run `gzip -d` and then `bgzip`. 

ScanExitronLR utilizes the `gffutils` package, which requires an SQL-lite database of the annotation file. You do not need to provide such a file, as ScanExitronLR will create one if one is not found, though it may take ~20 minutes to build.  It will be saved as `your_annotation.gtf.gz.db` in the same location as your annotation. In addition, we require a tabix index, and it will be created if one is not found. This should only take seconds.  It will be saved as `your_annotation.gtf.gz.tbi`. 

Thus, if you are running ScanExitronLR on a shared server and using a shared annotation database, you may not have writing privelages in the shared space. You will need to copy your annotation file to your local directory.

To run ScanExitron2 in extract mode, simply run

```bash
selr extract ...
```
with the following parameters

```
usage:
    		-i STR		REQUIRED: Input BAM file
    		-g STR		REQUIRED: Input genome reference (e.g. hg38.fa)
    		-r STR		REQUIRED: Input *sorted* and *bgzip'd* annotation reference (e.g. gencode_v38_sorted.gtf.gz).
    		-o STR		Output filename (e.g. bam_filename.exitron <- this is default)
    		-a/--ao INT	Reports only exitrons with AO of INT or above (default: 1).
    		-p/--pso FLOAT	Reports only exitrons with PSO of FLOAT or above (default: 0.01).
    		-c/--cores INT	Use INT cores (default: 1). Use as many as you can spare. Even large BAM files only use 4GB total memory on 10 cores.
    		-arabidopsis	Use this flag if using alignments from Arabidopsis. See github page for annotation file/genome assumptions.
    		-m/--mapq INT	Only considers reads with mapq score >= INT (default: 50)
    		-j/--jitter INT	Treat splice-sites with fuzzy boundry of +/- INT (default: 10).
    		-sr		Use this flag to skip the realignment step.
    		-sa		Use this flag to save isoform abundance files for downstream differential isoform usage analysis with LIQA.
    				Files are of the form: input.isoform.exitrons, input.isoform.normals
```

## Annotate

To run ScanExitron2 in annotate mode, simply run

```bash
selr annotate ...
```
with the following parameters

```
usage:
    		-i STR		REQUIRED: Input exitron file, generated from selr extract
    		-g STR		REQUIRED: Input genome reference (e.g. hg38.fa)
    		-r STR		REQUIRED: Input *sorted* and *gzip'd* annotation reference (e.g. gencode_v38_sorted.gtf.gz).
    		-o STR		Output filename (e.g. bam_filename.exitron.annotation <- this is default)
    		-b/--bam-file STR		If specified, annotation includes read supported NMD status directly from alignments.
    		-arabidopsis	Use this flag if using alignments from Arabidopsis. See github page for annotation file/genome assumptions.
```

The output is a tab-separated file.

# Example

See [here](https://github.com/ylab-hi/ScanExitronLR/tree/main/test_data) for an example. 

# Contact

Please feel free to post any issues here on github.

# Citation

TBD

