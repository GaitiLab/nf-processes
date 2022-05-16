#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process merge_fastqs {


       publishDir path: "${output_dir}/split_pipe/merged_fastqs/", mode: "symlink"

       input: 
       path input_dir
       path output_dir
       val sublibrary
    
       output: 
       tuple val(sublibrary), path("${sublibrary}_R1.fastq.gz"), path("${sublibrary}_R2.fastq.gz"), emit: sublibrary_read_pairs

       script: 
       """
       cat ${input_dir}/${sublibrary}*R1*fastq.gz > ${sublibrary}_R1.fastq.gz
       cat ${input_dir}/${sublibrary}*R2*fastq.gz > ${sublibrary}_R2.fastq.gz        
       """
}


process splitpipe_all {

       publishDir path: "${output_dir}/split_pipe/all/", mode: "copy"

       input: 
       tuple val(sublibrary), path(read_1), path(read_2)
       val kit
       val reference
       path sample_list
       path output_dir

       output: 
       path ("${sublibrary}/*"), emit: splitpipe_all_outputs
       path ("${sublibrary}/"), emit: splitpipe_all_dir

       script:
       """ 
       split-pipe --mode all \
       --kit ${kit} \
       --genome_dir ${reference} \
       --fq1 ${read_1} \
       --fq2 ${read_2} \
       --output_dir ${sublibrary}/ \
       --samp_list ${sample_list}
       """

}


process splitpipe_combine {

       publishDir path: "${output_dir}/split_pipe/", mode: "copy"

       input: 
       
       path sublibrary_path_list
       path output_dir
       val sample_name

       output: 
       path ("combined/"), emit: splitpipe_combined_by_sample
       tuple val(sample_name), path("combined/${sample_name}/DGE_filtered/DGE.mtx"), emit: sample_filtered_matrix


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
       
       if ( params.merge_fastqs ) {

        if ( sub_libraries instanceof List ) {
       samples_sublibraries = Channel.fromList(sub_libraries)
       } else if ( sub_libraries.contains('.txt') & sub_libraries instanceof String) {
       samples_sublibraries = Channel.fromList(file(sub_libraries).readLines())
       } else {
       println("If merging FASTQs, sublibraries must be either a list of strings of a .txt file with one sublibrary per line.")
       System.exit(1)
       }
       
       merge_fastqs(params.input_dir, params.output_dir, samples_sublibraries)
       splitpipe_all(merge_fastqs.out.sublibrary_read_pairs, params.kit, params.ref, params.sample_list, params.output_dir)
 
       } else {
       println("Not merging FASTQ files. Detecting sublibrary file pairs using the input directory and FASTQ pattern.")
       samples_sublibraries = Channel.fromFilePairs( params.input_dir + '/' + params.fastq_pattern, flat: true )
       splitpipe_all(samples_sublibraries, params.kit, params.ref, params.sample_list, params.output_dir)

       }
       if ( params.combine ) {


       sample_names = Channel.fromList(file(params.sample_list).readLines()).map { i -> spaceSplit(i)[0] }.view()
       splitpipe_combine(splitpipe_all.out.splitpipe_all_dir.collect(), params.output_dir, sample_name)
}
      
}
           


workflow {

     main: 

     pb_splitpipe(params.sublibrary)

}