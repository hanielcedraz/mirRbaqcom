# mirBaqcom – A Pipeline for Quality Control, Mapping, Identification, and Counting of miRNA Sequences

**mirBaqcom** is a user-friendly pipeline that implements automated processes for small RNA sequences quality control, mapping, prediction, identification, and counting. The pipeline includes:

- **Quality Control**: Uses Cutadapt to trim adapters and low-quality sequences (Martin, 2011), and UMI-tools for handling Unique Molecular Identifiers (Smith et al., 2017).
- **Mapping**: Utilizes Bowtie for sequence alignment (Langmead et al., 2009).
- **Identification and Counting**: Employs mirDeep2 to identify, predict, and count RNA sequences (Friedländer et al., 2012).

The main advantage of using mirBaqcom is its ability to analyze multiple samples simultaneously, generating grouped report outputs.

## Step 1: Download the Repository

Clone the repository to your preferred path:
```bash
$ git clone https://github.com/hanielcedraz/mirBaqcom.git
$ cd mirBaqcom
$ chmod +x mirBaqcom*.R
```

Alternatively, download the YAML file that contains the Conda environment:
```bash
$ conda env create -f baqcom_env.yaml
```

## Step 2: Preparing Sequence and Sample Files
Create a 00-Fastq folder and copy all your FASTQ files into this folder. Ensure your FASTQ files are named in a similar pattern:
SAMPLENAME_any_necessary_information_SE_001.fastq.gz

Create a file with sample and FASTQ file information using:
```bash
$ createSample.sh -f <single or paired>
```

## Step 3: Running mirBaqcom
### For UMI Removal
```bash
$ miRbaqcomQC-umiTools -h
```
### For Quality Control
```bash
$ miRbaqcomQC-CutAdapter -h
# Use -h to verify the trimming options.
```

### For Mapping with Bowtie
```bash
$ miRbaqcomMapping-Bowtie -h
```

### For Identifying, Predicting, and Counting miRNAs
```bash
$ miRbaqcomMapping-MirDeep2.pl -h
```

## Output Files
### For UMI Removal
A FASTQ file with UMIs in the header.
### For Quality Control
A trimmed FASTQ file for each sample (if paired-end, there will be two files: R1 and R2).
A .txt file with quality control information: total reads number, trimmed reads number, etc.
### For Mapping with Bowtie
A .sam file with mapped reads.
### For Identifying, Predicting, and Counting miRNAs
A known.counts file with known miRNAs.
A predict.counts file with predicted miRNAs.
A total.counts file with both predicted and known miRNAs.



### References
- Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, 17(1), pp. 10-12. doi:10.14806/ej.17.1.200
- Smith T, Heger A, Sudbery I. UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy. Genome Res. 2017 Mar;27(3):491-499. doi: 10.1101/gr.209601.116
- Langmead, B., Trapnell, C., Pop, M. et al. Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome Biol 10, R25 (2009). doi:10.1186/gb-2009-10-3-r25
- Friedländer MR, Mackowiak SD, Li N, Chen W, Rajewsky N. miRDeep2 accurately identifies known and hundreds of novel microRNA genes in seven animal clades. Nucleic Acids Res. 2012 Jan;40(1):37-52. doi: 10.1093/nar/gkr688








<!-- miRbaqcom
================

Download repository
``` bash
git clone https://github.com/hanielcedraz/miRbaqcom.git
```
Open miRbaqcom folder
``` bash
cd miRbaqcom
```
make pipelines executable
``` bash
chmod +x miRbaqcom*.R
```

Download the yaml file that contains the conda env: [baqcom.yaml](https://github.com/hanielcedraz/miRbaqcom/blob/master/mirBaqcom_env.yaml)

Create baqcom env from yaml
``` bash
conda env create -f baqcom_env.yaml
```

## GitHub Documents

This pipeline for miRNA-Seq analysis was based on the paper [A
bioinformatics approach to microRNA-sequencing
analysis](https://www.sciencedirect.com/science/article/pii/S266591312030131X),
as following.

## 

![Pipeline Figure](README_figs/README-pipeline.png)

## UMI Tools

### Instalation:

``` bash
$ pip install umi_tools
```

``` html
https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1601-4

https://www.encodeproject.org/microrna/microrna-seq/

https://academic.oup.com/nargab/article/3/3/lqab068/6325159

https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/miRNA_Pipeline/

https://www.sciencedirect.com/science/article/pii/S266591312030131X#fig1

https://cutadapt.readthedocs.io/en/stable/installation.html

https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC3245920/

https://github.com/OpenGene/fastp

https://cutadapt.readthedocs.io/en/stable/

https://rpubs.com/hanielcedraz/macs3pipeline

https://pip.pypa.io/en/stable/installation/

https://github.com/CGATOxford/UMI-tools

https://github.com/CGATOxford/UMI-tools/blob/master/doc/QUICK_START.md#step-4–mapping

https://github.com/salasouza/ecosystem

https://hub.docker.com/_/centos
``` -->

