#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/**
 ********************************** UNAP NextFlow ******************************************
 * 1 - Base calling, alligment and data preparation
 *  a. Base calling + alignment : ONT Dorado and minimap2
 *  b. Sorting, adding read group information and creating index : samtools
 * 2 - Long read variant calling : Clara Parabricks Deepvariant 1.5
 * 3 - Structural variants : Sniffles2
 *******************************************************************************************
*/


/**
* Inputs
*/

params.input = null
params.ref = null
params.tmp_dir = "tempDir"
params.out_dir = "output"
params.sort_threads = 32
params.sniffles_threads = 32
params.rParams = [] // Initialize an empty list to store -r parameters

// DeepVariant mode. Here we use ONT data as input so we fix it at "ont"
params.pbDVMode = "ont"
params.pbPATH = "pbrun"

params.help = null
params.test = null
params.reads = null

// Show help message
if (params.help) {
   log.info """
    Usage: nextflow run nextflow/main.nf  [options]

    Mandatory options:
        --input            Path to the directory containing Pod5 file for Dorado basecalling and minimap2 alligement
        --ref              Path to reference tar file with the same name as the FASTA file.
               For example, for a 'Homo_sapiens_assembly38.fasta' file with indexes the tar file should
               be 'Homo_sapiens_assembly38.fasta.tar'

    Options:
        --out_dir          Directory path to store all the outputs. [default : ${params.out_dir}]
        --tmp_dir          Directory to store temporary files. If it does not exists it will be created. [default : ${params.tmp_dir}]
        --sort_threads     Number of threads to use for samtools sort [default : ${params.sort_threads}]
    --sniffles_threads Numbers of threads to run sniffles2 [default : ${params.sniffles_threads}]
    --reads            samtools addreplacerg -r option. It should be specified as this example :  --reads '-r "SM:GM24385" -r "ID:GM24385"'

    Example:
      nextflow run main.nf  --input "/raid/pod5" --ref "/raid/parabricks_sample/Ref/Homo_sapiens_assembly38.fasta"
    """
    exit 0
}

// Verify that the mandatory parameters are provided
if (params.ref == null) error "The reference genome file is mandatory. Please specify it with --ref"
if (params.input == null) error "The path to the input POD5 files is mandatory, please specify it with --input_dir"

include { doradoMinimap } from './nextflow/basecalling.nf'
include { deepVariant } from './nextflow/pbDeepVariant.nf'
include { sniffles2 } from './nextflow/sniffles2.nf'

workflow {

    doradoMinimap(inputFile=params.input,
            inputRef=params.ref
        )

    deepVariant( inputBam=doradoMinimap.out.bam_file,
            inputBai=doradoMinimap.out.bai_file,
            inputRef=params.ref,
            pbDVMode=params.pbDVMode,
            pbPATH=params.pbPATH,
            tmpDir=params.tmp_dir
        )

    sniffles2( inputBam=doradoMinimap.out.bam_file,
            inputBai=doradoMinimap.out.bai_file,
            inputRef=params.ref
        )

}

workflow.onComplete {
    if(workflow.success) {
    println ( "This UNAP workflow is now complete!\n Your outpus are located in : " + params.out_dir )
    }
    else {
    println ( "Oops .. something went wrong, please look into the log file, and error messages into " + workDir )
    }
}

