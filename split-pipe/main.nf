#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process merge_fastqs {


       publishDir path: "${params.output_dir}/merged_fastqs/", mode: "symlink"

       input: 
       val sublibrary
    
       output: 
       path "${sublibrary}_R1.fastq.gz", emit: read_1_merged
       path "${sublibrary}_R2.fastq.gz", emit: read_2_merged

       script: 
       """
       cat ${params.input_dir}*R1*fastq.gz > ${sublibrary}_R1.fastq.gz
       cat ${params.input_dir}*R2*fastq.gz > ${sublibrary}_R2.fastq.gz        
       """
}


process splitpipe_all {

       publishDir path: "${params.output_dir}/split_pipe/sublibrary", mode: "copy"

       input: 
       path read_1
       path read_2
       val sublibrary

       output: 
       path ("split_pipe/*"), emit: splitpipe_all_outputs

       script:
       """ 
       split-pipe --mode ${params.mode} \
       --kit ${params.kit} \
       --genome_dir ${params.ref} \
       --fq1 ${read_1} \
       --fq2 ${read_2} \
       --output_dir split_pipe/ \
       --samp_list ${sample_list}
       """

}


process splitpipe_combine {

       publishDir path: "${params.output_dir}/split_pipe/combined", mode: "copy"

       input: 
       
       val sample_name
       path sublibrary_path_list

       output: 
       path("split_pipe/${sample_name}"), emit: splitpipe_combined_by_sample


       script: 
       """
       split-pipe
       --mode comb \
       --sublib_list ${sublibrary_path_list} \
       --output_dir split_pipe/${sample_name}
       """
}



workflow pb_splitpipe {

       main:

       if ( params.merge_fastqs ) {

       samples_sublibraries = Channel.fromList(params.sublibrary) 
       
       merge_fastqs(samples_sublibraries)
       // splitpipe_all(params.sample_list, merge_fastqs.out.read_1_merged, merge_fastqs.out.read_2_merged)
}

}
           


workflow {

     main: 

      pb_splitpipe()

}




