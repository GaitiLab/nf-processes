#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { addRecursiveSearch } from "../../utils/utils.nf"


process merge_fastqs {

       label 'splitpipe'

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

       label 'splitpipe'

       publishDir path: "${output_dir}/split_pipe/all/", mode: "copy"

       input: 
       tuple val(sublibrary), path(read_1), path(read_2)
       val chemistry
       val reference
       path sample_list
       path output_dir

       output: 
       path ("${sublibrary}/*"), emit: splitpipe_all_outputs
       path ("${sublibrary}/"), emit: splitpipe_all_dir

       script:
       """ 
       split-pipe --mode all \
       --chemistry ${chemistry} \
       --genome_dir ${reference} \
       --fq1 ${read_1} \
       --fq2 ${read_2} \
       --output_dir ${sublibrary}/ \
       --samp_list ${sample_list}
       """

}


process splitpipe_combine {

       label 'splitpipe'

       publishDir path: "${output_dir}/split_pipe/", mode: "copy"

       input: 
       
       path sublibrary_path_list
       path output_dir

       output: 
       path ("combined/"), emit: splitpipe_combined_by_sample


       script: 
       """
       split-pipe --mode comb \
       --sublibraries ${sublibrary_path_list} \
       --output_dir combined/
       """
}



workflow pb_splitpipe {

       main:
       
       if ( params.merge_fastqs ) {

       if ( params.sublibrary instanceof List ) {
       samples_sublibraries = Channel.fromList(params.sublibrary)
       } else if ( params.sublibrary.contains('.txt') & params.sublibrary instanceof String) {
       samples_sublibraries = Channel.fromList(file(params.sublibrary).readLines())
       } else {
       println("If merging FASTQs, sublibraries must be either a list of strings of a .txt file with one sublibrary per line.")
       System.exit(1)
       }
       
       merge_fastqs(params.input_dir, params.output_dir, samples_sublibraries)
       splitpipe_all(merge_fastqs.out.sublibrary_read_pairs, params.chemistry, params.ref, params.sample_list, params.output_dir)

       fastqs_out = merge_fastqs.out.sublibrary_read_pairs
       
       } else {
       println("Not merging FASTQ files. Detecting sublibrary file pairs using the input directory and FASTQ pattern.")
       samples_sublibraries = Channel.fromFilePairs( params.input_dir + '/' + addRecursiveSearch(params.recursive_search) + params.fastq_pattern, flat: true )
       .ifEmpty { exit 1, "Cannot find any reads matching: ${params.input_dir} with recursive set to: ${params.recursive_search}\n" }
       
       splitpipe_all(samples_sublibraries, params.chemistry, params.ref, params.sample_list, params.output_dir)

       fastqs_out = samples_sublibraries

       }
       if ( params.combine ) {

       splitpipe_combine(splitpipe_all.out.splitpipe_all_dir.collect(), params.output_dir)
}

       emit: 
       samples = fastqs_out
       paths = splitpipe_combine.out.splitpipe_combined_by_sample
      
}

workflow {

     main: 

     pb_splitpipe()

}
