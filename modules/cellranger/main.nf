#!/usr/bin/env nextflow

nextflow.enable.dsl=2


include { use_introns; concat_pattern_dir; addRecursiveSearch } from "../../utils/utils.nf"

process cellranger_count {

       label 'cellranger'

       publishDir path: "${output_dir}/cellranger_count/", mode: "copy"

       input: 
       tuple val(sample_name), path(input_dir)
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

       fastqs = Channel.fromFilePairs( params.input_dir + '/' + addRecursiveSearch(params.recursive_search) + params.cellranger.fastq_pattern )
       .ifEmpty { exit 1, "Cannot find any reads matching: ${params.input_dir}\n" }

       //strip the fastq names from the pairs input channel if no sample sheet is provided, split at _S
       // combine the sample names with the input directory
       samples_with_path = fastqs
                .flatMap { tuple ->
                tuple[0] = tuple[0].split(params.cellranger.sample_name_split)[0]
                }
                .unique()
                .combine(Channel.fromPath(params.input_dir))
       

} else {
       // if using a sample sheet, create a channel identical to above
       // requires two columns: one with the sample name, other with the path to the fastq files
       samples_with_path = Channel.fromPath( params.cellranger.sample_sheet )
       .ifEmpty { exit 1, "Cannot find any sample sheet at: ${params.cellranger.sample_sheet}\n" }
       .splitCsv(header:true)
       .map{ row-> tuple(row.sample_name, file(row.sample_path)) }  
}
               
       cellranger_count(samples_with_path, params.output_dir, params.cellranger.ref,
       params.cellranger.expected_cells, use_introns(params.cellranger.include_introns))   
}


workflow {

     main: 

      cellranger()

}



     
