#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process fastqc {

     publishDir path: "${params.output_dir}/fastqc/", mode: "copy"
     
     memory '4 GB'
   
     input:
     tuple val(name), path(reads)
     
     output: 
     path ("*.zip"), emit: fastqc_outputs
     path ("*.html"), emit: fastqc_htmls
     


     script: 
     """
     fastqc --outdir . \
     -t ${task.cpus} \
     ${reads.get(0)} \
     ${reads.get(1)}
     """

}


process multiqc {

     publishDir path: "${params.output_dir}/multiqc/", mode: "copy"

     input: 
     path fastqc_outputs


     output: 
     path('multiqc_report.html')
     path('multiqc_report_data/multiqc_general_stats.txt')

     script: 
     """
     mkdir -p ${params.output_dir}/multiqc/
     multiqc --force --interactive \
     --title "MultiQC_Output" \
     --filename "multiqc_report.html" \
     ${fastqc_outputs}
     """

}


workflow qc_workflow {

       main:

       fastqs = Channel.fromFilePairs( params.input_pattern )
       .ifEmpty { exit 1, "Cannot find any reads matching: ${params.input_pattern}\n" }

       fastqc(fastqs)

       multiqc(fastqc.out.fastqc_outputs.collect())

}


workflow {

     main: 

      qc_workflow()

}



      
