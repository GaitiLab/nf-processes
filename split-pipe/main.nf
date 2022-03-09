#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process merge_fastqs {


       publishDir path: "${params.output_dir}/merged_fastqs/", mode: "symlink"

       input: 
       val sublibrary
    
       output: 
       tuple val(sublibrary), path("${sublibrary}_R1.fastq.gz"), emit: read_1_merged
       tuple val(sublibrary), path("${sublibrary}_R2.fastq.gz"), emit: read_2_merged

       script: 
       """
       cat ${params.input_dir}/${sublibrary}*R1*fastq.gz > ${sublibrary}_R1.fastq.gz
       cat ${params.input_dir}/${sublibrary}*R2*fastq.gz > ${sublibrary}_R2.fastq.gz        
       """
}


process splitpipe_all {

       publishDir path: "${params.output_dir}/split_pipe/all/", mode: "copy"

       input: 
       tuple val(sublibrary), path(read_1), path(read_2)

       output: 
       path ("${sublibrary}/split_pipe/*"), emit: splitpipe_all_outputs

       script:
       """ 
       split-pipe --mode all \
       --kit ${params.kit} \
       --genome_dir ${params.ref} \
       --fq1 ${read_1} \
       --fq2 ${read_2} \
       --output_dir ${sublibrary}/split_pipe/ \
       --samp_list ${params.sample_list}
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
       --sublibraries ${sublibrary_path_list} \
       --output_dir split_pipe/${sample_name}
       """
}



workflow pb_splitpipe {

       main:

       if ( params.merge_fastqs ) {

       samples_sublibraries = Channel.fromList( params.sublibrary) 
       
       merge_fastqs(samples_sublibraries)
       splitpipe_all(merge_fastqs.out.read_1_merged.join(merge_fastqs.out.read_2_merged))
       // splitpipe_combine(NEED TO INSERT HERE CHANNEL FOR SAMPLE NAMES, splitpipe_all.out.splitpipe_all_outputs.collect())
}


}
           


workflow {

     main: 

      pb_splitpipe()

}




