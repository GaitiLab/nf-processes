#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process fastqc {

     publishDir path: "${params.output_dir}/fastqc/", mode: "copy"
   
     input:
     tuple val(name), path(read_1), path(read_2)
     
     output: 
     path ("*.zip"), emit: fastqc_outputs
     path ("*.html"), emit: fastqc_htmls
     


     script: 
     """
     fastqc --outdir . \
     -t ${task.cpus} \
     ${read_1} \
     ${read_2}
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
     --title ${params.multiqc_title} \
     --filename "multiqc_report.html" \
     ${fastqc_outputs}
     """

}


workflow qc_workflow {

       main:

       fastqs = Channel.fromFilePairs( params.input_pattern, flat:true )
       .ifEmpty { exit 1, "Cannot find any reads matching: ${params.input_pattern}\n" }

       fastqc(fastqs)
 
       if ( params.multiqc ) {

       multiqc(fastqc.out.fastqc_outputs.collect())

}      

}


workflow {

     main: 

      qc_workflow()

}



      
