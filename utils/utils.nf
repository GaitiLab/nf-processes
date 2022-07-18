#!/usr/bin/env nextflow

nextflow.enable.dsl=2

import java.nio.file.Path
import java.nio.file.Paths


def spaceSplit(string) { 
     string.split(" ")
     }

def makeCombinedMatrixPath(string_1, string_2) {
    Paths.get(string_1.toString() + "/" + string_2.toString()  + "/" + "DGE_filtered/DGE.mtx")
     }

def concat_pattern_dir(dir, pattern) { dir + '/' + pattern }


def use_introns (intron_param) { 
       if ( intron_param ) {

       introns = '--include-introns true'
       } else {
       introns = '--include-introns false'
}
       introns

}

def toTranspose(transpose) { 
     if (transpose) {
        transpose_addition = "-t" 
     } else {
         transpose_addition = ""
}
    transpose_addition
     }


def addRecursiveSearch(arg) {
     if (arg) {
          recursive = "*"
     } else {
          recursive = ""
     }
     recursive
}


/* create a channel from flat read pairs with a placeholder name and 
collected fastq names from a glob pattern
Steps: 
     - get the tail of the input channel to get just FASTQ names
     - flatten and collect all from the input channel
     - map into new channel with a placeholder name*/

def formatFASTQInputForFastQC(input_channel) {
     Channel.of("fastq_files").combine(input_channel.map( tuple -> tuple.tail()).flatten().collect()).
     map{ tuple -> [tuple[0], tuple.tail()]}
}

def setSingleforKallisto(single_read) {
     if (single_read) {
          single = "--single"
     } else {
          single = ""
     }
     single
}

/* If kallisto is set to single read, add the fragment length and deviation arguments provided by the params
if the mode is set paired, do not add (allow kallisto to calculate the length and deviation from the paired reads) 
*/


def addKallistoSingleParams(single_read, fragment_length, fragment_deviation) {

     if (single_read) {
          single_params = "-l ${fragment_length} -s ${fragment_deviation}"
     } else {
          single_params = ""
     }
     single_params

}

