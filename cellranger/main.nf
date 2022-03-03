#!/usr/bin/env nextflow

nextflow.enable.dsl=2


def concat_pattern_dir(dir, pattern) { dir + '/' + pattern }


def use_introns () { 
       if ( params.include_introns ) {

       introns = '--include-introns'
       } else {
       introns = ''
}
       introns

} 

intron_include = use_introns()

process cellranger_count {


       publishDir path: "${params.output_dir}/cellranger_count/", mode: "copy"

       input: 
       val sample_name
    
       output: 
       path "${sample_name}/*"
       path "${sample_name}/outs/*"

       script: 
       """
       cellranger count --id=${sample_name} \
                   --transcriptome=${params.transcriptome_ref} \
                   --fastqs=${params.input_dir} \
                   --sample=${sample_name} \
                   --expect-cells=${params.expected_cells} \
                   --localcores=8 \
                   --localmem=64 \
                   ${intron_include}
                   
       """
}

workflow cellranger {

       main:


       if ( params.sample_sheet == '' ) {

       fastqs = Channel.fromFilePairs( concat_pattern_dir(params.input_dir,params.fastq_pattern) )
       .ifEmpty { exit 1, "Cannot find any reads matching: ${params.input_dir}\n" }

       //strip the fastq names from the pairs input channel if no sample sheet is provided, split at _S
       stripped = fastqs
                .flatMap { tuple ->
                tuple[0] = tuple[0].split('_S')[0]
                }
                .unique()
          
       cellranger_count(stripped)

}
               

}


workflow {

     main: 

      cellranger()

}



     
