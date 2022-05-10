#!/usr/bin/env nextflow

import java.nio.file.Paths

nextflow.enable.dsl=2

if (params.transpose) {
   transpose_addition = "-t"
} else {
   transpose_addition = ""
}

binDir = Paths.get(workflow.projectDir.toString(), "bin/")

process scrublet {

        publishDir path: "${params.output_dir}/scrublet/${sample_name}/", mode: "copy"

        input: 
            tuple val(sample_name), val(matrix)

        output: 
            tuple val(sample_name), path("${sample_name}_scrublet_detection.csv"), emit: detections


        script:
        """
        python $binDir/scrublet_cell_prediction_CI.py -i ${matrix} -o ${sample_name}_scrublet_detection.csv \
        -e ${params.expected_rate} -mu ${params.min_counts} -mc ${params.min_cells} -g ${params.gene_variability} \
        -pc ${params.princ_components} $transpose_addition
        """       

}


workflow scrublet_remove {

       main: 

       matrices = Channel.fromPath( params.input_pattern )
       .ifEmpty { exit 1, "Cannot find any matrices matching: ${params.input_pattern}\n" }
       .map { file -> tuple(file.baseName, file) }

       scrublet(matrices)


}


workflow {

     main: 

      scrublet_remove()

}
