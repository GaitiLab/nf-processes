#!/usr/bin/env nextflow

nextflow.enable.dsl=2


include { use_introns; concat_pattern_dir } from "../../utils/utils.nf"

process cellranger_count {

       label 'cellranger'

       publishDir path: "${output_dir}/cellranger_count/", mode: "copy"

       input: 
       val sample_name
       path input_dir
       path output_dir
       val ref
       val expected_cells
       val introns
    
       output: 
       path "${sample_name}/*"
       path "${sample_name}/outs/*"
       tuple val(sample_name), path("${sample_name}/outs/filtered_feature_bc_matrix/matrix.mtx.gz"), emit: cellranger_filtered_matrix

       script: 
       """
       cellranger count --id=${sample_name} \
                   --transcriptome=${ref} \
                   --fastqs=${input_dir} \
                   --sample=${sample_name} \
                   --expect-cells=${expected_cells} \
                   --localcores=8 \
                   --localmem=64 \
                   ${introns}
                   
       """
}

workflow cellranger {

       main:


       if ( params.cellranger.sample_sheet == '' ) {

       fastqs = Channel.fromFilePairs( concat_pattern_dir(params.input_dir, params.cellranger.fastq_pattern) )
       .ifEmpty { exit 1, "Cannot find any reads matching: ${params.input_dir}\n" }

       //strip the fastq names from the pairs input channel if no sample sheet is provided, split at _S
       stripped = fastqs
                .flatMap { tuple ->
                tuple[0] = tuple[0].split('_S')[0]
                }
                .unique()
          
       cellranger_count(stripped, params.input_dir, params.output_dir, params.cellranger.ref,
       params.cellranger.expected_cells, use_introns())

}
               

}


workflow {

     main: 

      cellranger()

}



     
