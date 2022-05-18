#!/usr/bin/env nextflow

import java.nio.file.Paths

nextflow.enable.dsl=2


def toTranspose(transpose) { 
     if (transpose) {
        transpose_addition = "-t" 
     } else {
         transpose_addition = ""
}
    transpose_addition
     }



binDir = Paths.get(workflow.projectDir.toString(), "bin/")

process scrublet {

        publishDir path: "${params.output_dir}/scrublet/${sample_name}/", mode: "copy"

        input: 
            tuple val(sample_name), path(matrix)
            path output_dir
            val expected_rate
            val min_counts
            val min_cells
            val gene_variability
            val princ_components
            val transpose
  
        output: 
            tuple val(sample_name), path("${sample_name}_scrublet_detection.csv"), emit: detections


        script:
        """
        python $binDir/scrublet_cell_prediction_CI.py -i ${matrix} -o ${sample_name}_scrublet_detection.csv \
        -e ${expected_rate} -mu ${min_counts} -mc ${min_cells} -g ${gene_variability} \
        -pc ${princ_components} ${transpose}
        """       

}


workflow scrublet_remove {

       main: 

       matrices = Channel.fromPath( params.input_pattern )
       .ifEmpty { exit 1, "Cannot find any matrices matching: ${params.input_pattern}\n" }
       .map { file -> tuple(file.baseName, file) }

       scrublet(matrices, params.output_dir, params.expected_rate, params.min_counts, params.min_cells,
       params.gene_variability, params.princ_components, toTranspose(params.transpose))


}


workflow {

     main: 

      scrublet_remove()

}
