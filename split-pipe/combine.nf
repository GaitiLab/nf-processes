#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process splitpipe_combine {

       publishDir path: "${params.output_dir}/split_pipe/", mode: "copy"

       input: 

       path sublibrary_path_list

       output: 
       path("split_pipe/${sample_name}"), emit: splitpipe_combined_by_sample


       script: 
       """
       split-pipe
       --mode comb \
       --sublibraries ${sublibrary_path_list} \
       --output_dir combined/
       """
}


workflow {

     main: 

     
     sublibrary_paths = Channel.from( params.test ).flatMap { n -> [ n + ' ' + '\\' ] }
     splitpipe_combine(sublibrary_paths)



}
