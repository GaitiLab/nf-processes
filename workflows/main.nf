#!/usr/bin/env nextflow

nextflow.enable.dsl=2

import java.nio.file.Paths


include { fastqc; multiqc } from "../modules/fastqc_multiqc/main.nf"

include { pb_splitpipe } from "../modules/splitpipe/main.nf"

include { scrublet } from "../modules/scrublet/main.nf"


def spaceSplit(string) { 
     string.split(" ")
     }

          
def makeCombinedMatrixPath(string_1, string_2) {
     string_1 + "/" + string_2  + "/" + "DGE_filtered/DGE.mtx"
     }



workflow splitpipe_pb {
     
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
       
       merge_fastqs(samples_sublibraries)
       fastqc(merge_fastqs.out.sublibrary_read_pairs)
       splitpipe_all(merge_fastqs.out.sublibrary_read_pairs)
 
       } else {
       println("Not merging FASTQ files. Detecting sublibrary file pairs using the input directory and FASTQ pattern.")
       samples_sublibraries = Channel.fromFilePairs( params.input_dir + '/' + params.fastq_pattern, flat: true )
       fastqc(samples_sublibraries)
       splitpipe_all(samples_sublibraries)
       }
       
       multiqc(fastqc.out.fastqc_outputs.collect())
       if ( params.combine ) {

       splitpipe_combine(splitpipe_all.out.splitpipe_all_dir.collect())
}
       
       // read in the sample names to find them in the combined output directory for scrublet
       
       sample_names = Channel.fromList(file(params.sample_list).readLines()).map { i -> spaceSplit(i)[0] }

       combined_output_matrices = sample_names.map { i -> [ i, makeCombinedMatrixPath(output_path, i)] }

       scrublet(combined_output_matrices)

}

