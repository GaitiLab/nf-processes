#!/usr/bin/env nextflow

nextflow.enable.dsl=2

import java.nio.file.Path
import java.nio.file.Paths

binDir = Paths.get(workflow.projectDir.toString(), "bin/")


include { fastqc; multiqc; } from "../../modules/fastqc_multiqc/main.nf"

include { merge_fastqs; splitpipe_all; splitpipe_combine } from "../../modules/split-pipe/main.nf"

include { scrublet } from "../../modules/scrublet/main.nf"

include { spaceSplit; makeCombinedMatrixPath; concat_pattern_dir; use_introns; toTranspose } from "../../utils/utils.nf"


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

       fastqc(merge_fastqs.out.sublibrary_read_pairs, params.output_dir)
 
       } else {
       println("Not merging FASTQ files. Detecting sublibrary file pairs using the input directory and FASTQ pattern.")
       samples_sublibraries = Channel.fromFilePairs( concat_pattern_dir(params.input_dir, params.fastq_pattern), flat: true )
       .ifEmpty { exit 1, "Cannot find any reads matching: ${params.fastq_pattern} in ${params.input_dir}\n" }
       fastqc(samples_sublibraries, params.output_dir)
       split_pipe_input = samples_sublibraries
       }
       
       splitpipe_all(split_pipe_input, params.kit, params.ref, params.sample_list, params.output_dir)

       if ( params.multiqc ) {
       multiqc(fastqc.out.fastqc_outputs.collect(), params.output_dir, params.multiqc_title)

       }

       if ( params.combine ) {

     // create a channel with sample name and the path of the combined output matrix (cast as string)

     splitpipe_combine(splitpipe_all.out.splitpipe_all_dir.collect(), params.output_dir)

     combine_out = splitpipe_combine.out.splitpipe_combined_by_sample

     sample_names = Channel.fromList(file(params.sample_list).readLines()).map { i -> spaceSplit(i)[0] }

     combined_output_matrices = sample_names.combine(combine_out).map { i, j -> [ i,
      makeCombinedMatrixPath(j, i)] }

     scrublet(combined_output_matrices, params.output_dir, params.expected_rate, params.min_counts, params.min_cells,
       params.gene_variability, params.princ_components, toTranspose(params.transpose))

}

}

workflow {
     main:
     splitpipe_pb_extended()
}

