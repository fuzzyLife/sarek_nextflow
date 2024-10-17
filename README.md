# [RNA-seq](https://rnabio.org) analysis with [nextflow](https://nf-co.re/rnaseq/3.16.1/)


## Create directory or simply clone this [repo](https://github.com/animesh/rnaseq_nextflow)
```
git clone https://github.com/animesh/rnaseq_nextflow
cd rnaseq_nextflow/
```

## Download raw [data](https://rnabio.org/module-01-inputs/0001/05/01/RNAseq_Data/)
```
wget http://genomedata.org/rnaseq-tutorial/HBR_UHR_ERCC_ds_5pc.tar
tar -xvf HBR_UHR_ERCC_ds_5pc.tar
```

## Create [sample-sheet](https://nf-co.re/rnaseq/usage#samplesheet-input) from the downloaded Data 
```
echo "sample,fastq_1,fastq_2,strandedness" > samples.csv
ls -1 *1.fastq.gz | awk -F "_" '{print $1 $2}' > c0
ls -1 $PWD/*1.fastq.gz > c1
ls -1 $PWD/*2.fastq.gz > c2
printf 'auto\n%.0s' {1..`ls *1.fastq.gz`} > c3
paste -d "," c? >> samples.csv
cat samples.csv
```

## Download [Reference Genome](https://rnabio.org/module-01-inputs/0001/02/01/Reference_Genomes/)
```
wget http://genomedata.org/rnaseq-tutorial/fasta/GRCh38/chr22_with_ERCC92.fa
```
## And the [annotation for the Reference Genome](https://rnabio.org/module-01-inputs/0001/03/01/Annotations/)
```
wget http://genomedata.org/rnaseq-tutorial/annotations/GRCh38/chr22_with_ERCC92.gtf
sed "s/exon_number \"1\";$/exon_number \"1\";gene_biotype "protein_coding";/g" chr22_with_ERCC92.gtf > chr22_with_ERCC92.biotype.gtf
```

## Install [nextflow](https://www.nextflow.io/docs/latest/install.html)
Pre-reqs: [curl](https://ec.haxx.se/install/linux.html), [Java](https://askubuntu.com/questions/1492571/install-openjdk-21) and [Docker](https://docs.docker.com/engine/install/ubuntu/) or [Singularity](https://www.youtube.com/watch?v=OB1qs56g0VA) on Linux, Windows users can run via [WSL](https://www.nextflow.io/blog/2021/setup-nextflow-on-windows.html)
```
curl -s https://get.nextflow.io | bash
```
```

      N E X T F L O W
      version 24.04.4 build 5917
      created 01-08-2024 07:05 UTC (09:05 CEST)
      cite doi:10.1038/nbt.3820
      http://nextflow.io


Nextflow installation completed. Please note:
- the executable file `nextflow` has been created in the folder: /home/ash022/rnaseq_nextflow
- you may complete the installation by moving it to a directory in your $PATH
```
## [test](https://www.nf-test.com/)
```
./nextflow run hello
```
```

 N E X T F L O W   ~  version 24.04.4

Launching `https://github.com/nextflow-io/hello` [trusting_cori] DSL2 - revision: afff16a9b4 [master]

executor >  local (4)
executor >  local (4)
[4d/db3896] sayHello (1) [100%] 4 of 4 ✔
Ciao world!

Hello world!

Hola world!

Bonjour world!
```

## can finally test nextflow-rnaseq-pipeline 
```
./nextflow run nf-core/rnaseq -profile docker,test --outdir test
```

## if it goes well, run it over the downloaded data above using the created [sample-sheet](./samples.csv), [reference genome](./chr22_with_ERCC92.fa) and the [corresponding annotation](./chr22_with_ERCC92.biotype.gtf)
```
./nextflow run nf-core/rnaseq --max_memory '16.GB' --max_cpus 6  --input samples.csv --outdir results --gtf chr22_with_ERCC92.biotype.gtf --fasta chr22_with_ERCC92.fa -profile docker 
```

## AND in case above fails, check the *.nextflow* logs, lets say [cosmic rays flipped the bits](https://radiolab.org/podcast/bit-flip), one can try to [resume](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html) from last successfull step/checkpoint just with a -resume switch!
```
./nextflow run nf-core/rnaseq --max_memory '16.GB' --max_cpus 6  --input samples.csv --outdir results --gtf chr22_with_ERCC92.biotype.gtf --fasta chr22_with_ERCC92.fa -profile docker -resume
```
but if nothing seems to be working, try https://github.com/nf-core/rnaseq/issues 🤞


## [QC](https://rnabio.org/module-02-alignment/0002/06/01/Alignment_QC/#:~:text=%24RNA_ALIGN_DIR%0Amultiqc%20./-,MultiQC%20screenshot,-View%20a%20pre) 
```
wget https://rnabio.org/assets/module_2/multiqc.png
```
looks [similar](./Screenshot%202024-10-17%20174136.png) to nextflow [multiqc_report](results/multiqc/star_salmon/multiqc_report.html#rseqc)



## [Result](https://rnabio.org/module-09-appendix/0009/09/03/POSIT_Setup/) 
compare [original](http://genomedata.org/rnaseq-tutorial/results/cshl2022/rnaseq/gene_read_counts_table_all_final.tsv) 
```
wget http://genomedata.org/rnaseq-tutorial/results/cshl2022/rnaseq/gene_read_counts_table_all_final.tsv	
```
with resulting [nextflow-gene-counts](./results/star_salmon/salmon.merged.gene_counts.tsv), spearman rank-correlation is ~ 0.99 calculated using [Perseus](https://www.nature.com/articles/nmeth.3901) shown in blue in [scatter-plot](./Screenshot%202024-10-17%20170444.png) and Euclidean distance [clusters](./Screenshot%202024-10-17%20165609.png)) accordingly, [more](./session1.sps) on this in blog/[post](https://open.substack.com/pub/fuzzylife/p/rna-seq-analysis-with-nextflow?r=a55q5&utm_campaign=post&utm_medium=web&showWelcomeOnShare=true) https://fuzzylife.substack.com/p/rna-seq-analysis-with-nextflow

#### For mapping to the latest whole Human-genome assembly, need to download and run nextflow the appropriate Reference Genome and its corresponding Annotation, for example
```
wget http://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
wget http://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz
./nextflow run nf-core/rnaseq --max_memory '62.GB' --max_cpus 16  --input samples.csv --outdir results --gtf Homo_sapiens.GRCh38.112.gtf.gz --fasta Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz -profile docker
```

## cleanup 
```
rm -rf work .nextflow c? HBR_UHR_ERCC_ds_5pc.tar
```

## issues faced (so far...)

#### docker daemon not running?
```
sudo touch nohup.out
sudo nohup dockerd &
```

#### crash? maybe due to RAM, usually due to STAR genome-indexing, try to increase allocated RAM
```
./nextflow run nf-core/rnaseq --max_memory '92.GB' --max_cpus 16  --input samples.csv --outdir results --gtf Homo_sapiens.GRCh38.112.gtf.gz --fasta Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz -profile docker
```

