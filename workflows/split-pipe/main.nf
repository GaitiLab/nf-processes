#!/usr/bin/env nextflow

nextflow.enable.dsl=2

import java.nio.file.Path
import java.nio.file.Paths


include { fastqc; multiqc; } from "../../modules/fastqc_multiqc/main.nf"

include { merge_fastqs; splitpipe_all; splitpipe_combine } from "../../modules/split-pipe/main.nf"

include { scrublet; } from "../../modules/scrublet/main.nf"


def spaceSplit(string) { 
     string.split(" ")
     }

workflow splitpipe_pb_extended {
     
     main: 

     if ( params.merge_fastqs ) {

     if ( sub_libraries instanceof List ) {
       samples_sublibraries = Channel.fromList(sub_libraries)
       } else if ( sub_libraries.contains('.txt') & sub_libraries instanceof String) {
       samples_sublibraries = Channel.fromList(file(sub_libraries).readLines())
       } else {
       println("If merging FASTQs, sublibraries must be either a list of strings or a .txt file with one sublibrary per line.")
       System.exit(1)
       }
       
       merge_fastqs(params.input_dir, params.output_dir, samples_sublibraries)
       split_pipe_input = merge_fastqs.out.sublibrary_read_pairs

       fastqc(merge_fastqs.out.sublibrary_read_pairs)
 
       } else {
       println("Not merging FASTQ files. Detecting sublibrary file pairs using the input directory and FASTQ pattern.")
       samples_sublibraries = Channel.fromFilePairs( params.input_dir + '/' + params.fastq_pattern, flat: true )
       fastqc(samples_sublibraries)
       split_pipe_input = samples_sublibraries
       }
       
       splitpipe_all(split_pipe_input, params.kit, params.ref, params.sample_list, params.output_dir)

       if ( params.multiqc ) {
       multiqc(fastqc.out.fastqc_outputs.collect())

       }

       if ( params.combine ) {


     // create a channel with sample name and the path of the combined output matrix (cast as string)

     sample_names = Channel.fromList(file(params.sample_list).readLines()).map { i -> spaceSplit(i)[0] }


     splitpipe_combine(splitpipe_all.out.splitpipe_all_dir.collect(), params.output_dir, sample_names)
       

     scrublet(splitpipe_combine.out.sample_filtered_matrix)

}

}

workflow {
     main:
     splitpipe_pb_extended()
}

