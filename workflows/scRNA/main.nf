#!/usr/bin/env nextflow

nextflow.enable.dsl=2

import java.nio.file.Path
import java.nio.file.Paths


include { fastqc; multiqc; } from "../../modules/fastqc_multiqc/main.nf"

include { merge_fastqs; splitpipe_all; splitpipe_combine } from "../../modules/split-pipe/main.nf"

include { scrublet; } from "../../modules/scrublet/main.nf"

include {concat_pattern_dir; use_introns; cellranger_count } from "../../modules/cellranger/main.nf"

binDir = Paths.get(workflow.projectDir.toString(), "bin/")


def spaceSplit(string) { 
     string.split(" ")
     }


workflow scRNA {
     
     main: 

     
     if ( params.method != "cellranger" || params.method != "split-pipe") {
          println("Error: incorrect scRNA method specified. Please select one of cellranger or split-pipe.")
          System.exit(1)
     }

     if (params.method == "cellranger") {

          if ( params.sample_sheet == '' ) {

       fastqs = Channel.fromFilePairs( concat_pattern_dir(params.input_dir, params.fastq_pattern) )
       .ifEmpty { exit 1, "Cannot find any reads matching: ${params.fastq_pattern} in ${params.input_dir}\n" }

       //strip the fastq names from the pairs input channel if no sample sheet is provided, split at _S
       stripped = fastqs
                .flatMap { tuple ->
                tuple[0] = tuple[0].split('_S')[0]
                }
                .unique()
          
       cellranger_count(stripped, params.input_dir, params.output_dir, params.ref, params.expected_cells, use_introns())

       scrublet_input = cellranger_count.out.cellranger_filtered_matrix

}
     } else if ( params.method == "split-pipe" ) {

     if ( params.merge_fastqs ) {

     if ( sub_libraries instanceof List ) {
       samples_sublibraries = Channel.fromList(sub_libraries)
       .ifEmpty { exit 1, "Cannot find any sample names in : ${sublibraries}\n" }
       } else if ( sub_libraries.contains('.txt') & sub_libraries instanceof String) {
       samples_sublibraries = Channel.fromList(file(sub_libraries).readLines())
       .ifEmpty { exit 1, "Cannot find any sample names in : ${sub_libraries}\n" }
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
       .ifEmpty { exit 1, "Cannot find any reads matching: ${params.input_dir}\n" }
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
     .ifEmpty { exit 1, "Could not parse sample names from : ${params.sample_list}\n" }


     splitpipe_combine(splitpipe_all.out.splitpipe_all_dir.collect(), params.output_dir)

     scrublet_input = splitpipe_combine.out.sample_filtered_matrix

     }


     scrublet(scrublet_input)

}

}

workflow {
     main:
     scRNA()
}

