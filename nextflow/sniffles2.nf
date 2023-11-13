process sniffles2 {
    label 'localGPU'
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"

    input:
        path inputBam
        path inputBai
        path inputRef

    output:
        path "${inputBam.baseName}.sniffles2.vcf"

    script:
        """
	tar xf ${inputRef.Name} && \
        sniffles --threads ${params.sniffles_threads} --allow-overwrite \
                    --reference ${inputRef.baseName}  \
                    --input ${inputBam} \
                    --vcf ${inputBam.baseName}.sniffles2.vcf
        """

}
