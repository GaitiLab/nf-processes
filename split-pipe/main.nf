#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process merge_fastqs {


       publishDir path: "${params.output_dir}/split_pipe/merged_fastqs/", mode: "symlink"

       input: 
       val sublibrary
    
       output: 
       tuple val(sublibrary), path("${sublibrary}_R1.fastq.gz"), path("${sublibrary}_R2.fastq.gz"), emit: sublibrary_read_pairs

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
       path ("${sublibrary}/*"), emit: splitpipe_all_outputs
       path ("${sublibrary}/"), emit: splitpipe_all_dir

       script:
       """ 
       split-pipe --mode all \
       --kit ${params.kit} \
       --genome_dir ${params.ref} \
       --fq1 ${read_1} \
       --fq2 ${read_2} \
       --output_dir ${sublibrary}/ \
       --samp_list ${params.sample_list}
       """

}


process splitpipe_combine {

       publishDir path: "${params.output_dir}/split_pipe/", mode: "copy"

       input: 
       
       path sublibrary_path_list

       output: 
       path("combined/"), emit: splitpipe_combined_by_sample


       script: 
       """
       split-pipe --mode comb \
       --sublibraries ${sublibrary_path_list} \
       --output_dir combined/
       """
}



workflow pb_splitpipe {

       take: 
       sub_libraries

       main:

       if ( sub_libraries instanceof List ) {
       samples_sublibraries = Channel.fromList(sub_libraries)
       } else if ( sub_libraries.contains('.txt') & sub_libraries instanceof String) {
       samples_sublibraries = Channel.fromList(file(sub_libraries).readLines())
       } else {
       println("sublibraries must be either a Groovy list of strings of a .txt file with one sublibrary per line.")
       System.exit(1)
       }
       
       if ( params.merge_fastqs ) {
       
       merge_fastqs(samples_sublibraries)
       splitpipe_all(merge_fastqs.out.sublibrary_read_pairs)
 
       } else {
       splitpipe_all(samples_sublibraries)
       }
       if ( params.combine ) {

       splitpipe_combine(splitpipe_all.out.splitpipe_all_dir.collect())
}
      
}
           


workflow {

     main: 

     
     pb_splitpipe(params.sublibrary)

}




