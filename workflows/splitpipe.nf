#!/usr/bin/env nextflow

nextflow.enable.dsl=2

import java.nio.file.Paths


include { fastqc; multiqc } from "../modules/fastqc_multiqc/main.nf"

include { pb_splitpipe } from "../modules/splitpipe/main.nf"


workflow splitpipe_pb {

}

