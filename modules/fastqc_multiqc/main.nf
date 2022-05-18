#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process fastqc {

     publishDir path: "${output_dir}/fastqc/", mode: "copy"
   
     input:
     tuple val(name), path(read_1), path(read_2)
     path output_dir
     
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

     publishDir path: "${output_dir}/multiqc/", mode: "copy"

     input: 
     path fastqc_outputs
     path output_dir
     val multiqc_title



     output: 
     path('multiqc_report.html')
     path('multiqc_report_data/multiqc_general_stats.txt')

     script: 
     """
     mkdir -p ${output_dir}/multiqc/
     multiqc --force --interactive \
     --title ${multiqc_title} \
     --filename "multiqc_report.html" \
     ${fastqc_outputs}
     """

}


workflow qc_workflow {

       main:

       fastqs = Channel.fromFilePairs( params.input_pattern, flat:true )
       .ifEmpty { exit 1, "Cannot find any reads matching: ${params.input_pattern}\n" }

       fastqc(fastqs, params.fastqc.output_dir)
 
       if ( params.multiqc ) {

       multiqc(fastqc.out.fastqc_outputs.collect(), params.fastqc.output_dir,params.multiqc.multiqc_title)

}      

}


workflow {

     main: 

      qc_workflow()

}



      
