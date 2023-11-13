process deepVariant {
    label 'localGPU'
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"

    //stageInMode "copy"

    input:
    path inputBam
    path inputBai
    path inputRef
    val pbDVMode
    val pbPATH
    val tmpDir

    output:
      path "${inputBam.baseName}.dv.vcf"

    script:
    """
    tar xf ${inputRef.Name} && \
    mkdir -p ${tmpDir} && \
    time ${pbPATH} deepvariant \
    --tmp-dir ${tmpDir} \
    --in-bam ${inputBam} \
    --ref ${inputRef.baseName} \
    --out-variants ${inputBam.baseName}.dv.vcf \
    --mode ${pbDVMode} \
    --run-partition --norealign-reads
    """

}
