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


def use_introns () { 
       if ( params.include_introns ) {

       introns = '--include-introns'
       } else {
       introns = ''
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




