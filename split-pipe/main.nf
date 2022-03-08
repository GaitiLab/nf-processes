#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process merge_fastqs {


       publishDir path: "${params.input_dir}/merged_fastqs/", mode: "copy"

       input: 
       val sample_name
    
       output: 
       path "${params.input_dir}/merged_fastqs/${sample_name}_R1.fastq.gz", emit: read_1_merged
       path "${params.input_dir}/merged_fastqs/${sample_name}_R2.fastq.gz", emit: read_2_merged

       script: 
       """
       cat ${params.input_dir}/*R1*fastq* > ${sample_name}_R1.fastq.gz
       cat ${params.input_dir}/*R2*fastq* > ${sample_name}_R2.fastq.gz        
       """
}


process splitpipe_all {

       publishDir path: "${params.output_dir}/splitpipe/", mode: "copy"

       input: 
       path read_1
       path read_2
       path sample_list

       output: 
       
       path("${params.output_dir}/splitpipe/*")


       script: 
       """
       split-pipe \ 
       --mode all \
       --kit WT \
       --genome_dir ${params.ref} \ 
       --fq1 ${read_1} \
       --fq2 ${read_2} \
       --output_dir split pipe/ \ 
       --samp_list ${sample_list}
       """

{

workflow pb_splitpipe {

       main:

       if ( params.merge_fastqs ) {
       
       merge_fastqs(params.sample_name)
       // splitpipe_all(params.sample_list, merge_fastqs.out.read_1_merged, merge_fastqs.out.read_2_merged)

}


}
               

}


workflow {

     main: 

      cellranger()

}




