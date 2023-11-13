# UNAP Nextflow

## Overview

This NextFlow workflow is that latest iteration of NVIDIA’s [UNAP workflow](https://developer.nvidia.com/blog/boosting-ultra-rapid-nanopore-sequencing-analysis-on-nvidia-dgx-a100/), designed for high speed end-to-end Oxford Nanopore germline sequencing analysis on GPU.

This workflow includes all the steps that are needed to run :  
1. Basecalling and integrated alignment with Dorado/Minimap2
2. Variant calling with [DeepVariant](https://docs.nvidia.com/clara/parabricks/latest/documentation/tooldocs/man_deepvariant.html) in Parabricks
3. Structural variant calling with Sniffles2.
  
For more information on the NVIDIA Parabricks portion of this workflow, you can go to our [webpage](https://www.nvidia.com/en-gb/clara/genomics/), our latest [docs page](https://docs.nvidia.com/clara/parabricks/latest/index.html), or read the latest [developer blogs on Parabricks](https://developer.nvidia.com/blog/search-posts/?q=parabricks).
  
As input you will need an ONT Pod5, ONT R10.4 data file.
You will be able to output 
* The Bam file with index after the base calling and sorting. 
* The variants as a VCF file by Parabricks Deepvariant
* The structural variants as a VCF file by Sniffles2

All the outputs are stored in *params.out_dir* (default to *$PWD/output*)

## The Software 

### Nextflow 

First you need to install Nextflow, version 23.04.1 or above

https://nextflow.io/index.html#GetStarted


### Docker containers

We will need three containers for this project :

* NVIDIA Clara Parabricks

```
docker pull nvcr.io/nvidia/clara/clara-parabricks:4.2.0-1
```

* ONT Dorado 

```
docker pull nanoporetech/dorado:sha06c3f67311887d7976509748ae9357788fafce94
```

* Sniffles2
Then you will build the sniffles2 container
```
cd docker
docker build -t sniffles2:latest .
```

## Running the pipleine 

This workflow works with POD5 file as input

### Help 
You can print the help using the option  :
`run nextflow/main.nf --help` 

```
$ nextflow run nextflow/main.nf --help
N E X T F L O W  ~  version 23.04.1
Launching `nextflow/main.nf` [sad_poitras] DSL2 - revision: e4e4a44464

    Usage: nextflow run nextflow/main.nf  [options]

    Mandatory options:
        --input            Path to the directory containing Pod5 file for Dorado basecalling and minimap2 alligement
        --ref              Path to reference tar file with the same name as the FASTA file.
                           For example, for a 'Homo_sapiens_assembly38.fasta' file with indexes the tar file should
                           be 'Homo_sapiens_assembly38.fasta.tar'

    Options:
        --out_dir          Directory path to store all the outputs. [default : output]
        --tmp_dir          Directory to store temporary files. If it does not exists it will be created. [default : tempDir]
        --sort_threads     Number of threads to use for samtools sort [default : 32]
        --sniffles_threads Numbers of threads to run sniffles2 [default : 32]
        --reads            samtools addreplacerg -r option. It should be specified as this example :  --reads '-r "SM:GM24385" -r "ID:GM24385"'

    Example:
      nextflow run main.nf  --input "/raid/pod5" --ref "/raid/parabricks_sample/Ref/Homo_sapiens_assembly38.fasta"

```

### Run with a sample

```
# Example command :
nextflow run -c config/local.nf.conf \
    main.nf \
    --input "Data/Pod5" \
    --ref "parabricks_sample/Homo_sapiens_assembly38.fasta.tar"  \
    --reads '-r "SM:GM24385" -r "ID:GM24385"'
```
