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

        publishDir path: "${params.output_dir}/${matrix.baseName}/scrublet/", mode: "copy"

        input: 
            path(matrix)

        output: 
            path "${matrix.baseName}_scrublet_detection.csv"


        script:
        """
        python $binDir/scrublet_cell_prediction_CI.py -i ${matrix} -o ${matrix.baseName}_scrublet_detection.csv \
        -e ${params.expected_rate} -mu ${params.min_counts} -mc ${params.min_cells} -g ${params.gene_variability} \
        -pc ${params.princ_components} $transpose_addition
        """       

}


workflow scrublet_remove {

       main:

       matrices = Channel.fromPath( params.input_pattern )
       .ifEmpty { exit 1, "Cannot find any matrices matching: ${params.input_pattern}\n" }

       scrublet(matrices)


}


workflow {

     main: 

      scrublet_remove()

}
