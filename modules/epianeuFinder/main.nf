#!/usr/bin/env nextflow

nextflow.enable.dsl=2

import java.nio.file.Paths

binDir = Paths.get(workflow.projectDir.toString(), "bin/")


process epianeuFinder {


       publishDir path: "${params.output_dir}/", mode: "copy"

       input: 
       tuple val(sample_name), path(input_fragments)
       each min_fragment_number  
         

       output: 
       tuple val(sample_name), path("${sample_name}_${min_fragment_number}/"), emit: epianeufinder_dir
       tuple val(sample_name), path("${sample_name}_${min_fragment_number}/results_table.tsv"), emit: epianeuFinder_results

       script: 
       """
       Rscript $binDir/run_epianeuFinder_HPC.R -i ${input_fragments} \
       -o ${sample_name}_${min_fragment_number}/ \
       -b ${params.blacklist} \
       -m ${min_fragment_number} \
       -g ${params.ref_genome} \
       -w ${params.window} \
       -n ${params.numcores} \
       -t ${params.title}                              
       """
}

workflow run_epianeuFinder_HPC {

       main: 

       fragment_paths = Channel.fromPath( params.input_fragments )
       .ifEmpty { exit 1, "Cannot find any matrices matching: ${params.input_pattern}\n" }
       .splitCsv(header:true)
       .map{ row-> tuple(row.sample_name, file(row.input_fragments)) }       

       possible_min_frags = Channel.fromList( params.minfrags )       

       epianeuFinder(fragment_paths, possible_min_frags)


}

workflow {

     main: 

     run_epianeuFinder_HPC()

}




