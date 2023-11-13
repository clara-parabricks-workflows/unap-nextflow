process doradoMinimap {
    label 'localGPU'
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"

    input:
        path inputFile
        path inputRef
	

    output:
        path "${inputFile.baseName}.sorted.index.bam", emit: bam_file
        path "${inputFile.baseName}.sorted.index.bam.bai", emit: bai_file

    script:
	def r_args = params.reads ?: ''
	
        """
	tar xf ${inputRef.Name} && \
        dorado basecaller /models/dna_r10.4.1_e8.2_400bps_hac@v4.1.0 ${inputFile} --reference ${inputRef.baseName} | samtools sort --threads ${params.sort_threads} -O BAM -o ${inputFile.baseName}.sorted.bam -

        samtools addreplacerg ${r_args} \
                -@10 -o ${inputFile.baseName}.sorted.index.bam \
                ${inputFile.baseName}.sorted.bam

        samtools index -@10 -b ${inputFile.baseName}.sorted.index.bam
        """

}
