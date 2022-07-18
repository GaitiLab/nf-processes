#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { addRecursiveSearch; setSingleforKallisto; addKallistoSingleParams } from "../../utils/utils.nf"


process kallisto_quant {
     
     label 'kallisto'

     publishDir path: "${params.output_dir}/kallisto_quant/${name}/", mode: "copy"
   
     input:
     tuple val(name), path(reads)
     path transcriptome_index
     val single_reads
     val single_read_params

     output: 
     path "abundance.h5", emit: cell_gene_matrix
     path "abundance.tsv", emit: cell_gene_tsv
     path "run_info.json", emit: run_information


     script: 
     """
     kallisto quant -i ${transcriptome_index} -o . ${single_reads} ${single_read_params} ${reads} 
     """

}


workflow kallisto_qulk {

       main:

       fastqs = Channel.fromFilePairs( params.input_dir + '/' + addRecursiveSearch(params.recursive_search) + params.fastq_pattern )
       .ifEmpty { exit 1, "Cannot find any reads matching: ${params.input_dir} with recursive set to: ${params.recursive_search}\n" }

       kallisto_quant(fastqs, params.transcriptome_index, setSingleforKallisto(params.single_read), addKallistoSingleParams(
          params.single_read, params.fragment_length, params.fragment_deviation))
}


workflow {

     main: 

      kallisto_qulk()

}



      
