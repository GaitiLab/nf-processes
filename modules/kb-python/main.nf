#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { addRecursiveSearch; getKBPythonParity } from "../../utils/utils.nf"


process kb_ref {

     label 'kb_python'
   
     publishDir path: "${params.output_dir}/kb_python_ref/", mode: "copy"

     input: 
     path cdna_fasta
     path primary_fasta
     path annotation_gtf


     output: 
     path "transcriptome.idx", emit: index
     path "transcripts_to_genes.txt", emit: gene_map


     script:
     """
     kb ref -i transcriptome.idx -g transcripts_to_genes.txt -f1 ${cdna_fasta} ${primary_fasta} ${annotation_gtf}
     """

}

process kb_count {

     label 'kb_python'

     publishDir path: "${params.output_dir}/kallisto_count/${name}/", mode: "copy"
   
     input:
     tuple val(name), path(reads)
     path transcriptome_index
     path gene_map
     val chemistry
     val max_memory
     val threads
     val parity
     val strand


     output: 
     path "inspect.json"
     path "matrix.ec"
     path "output.bus"
     path "output.unfiltered.bus"
     path "run_info.json"
     path "transcripts.txt"
     path "counts_unfiltered/cells_x_genes.barcodes.txt", emit: cell_gene_barcode
     path "counts_unfiltered/cells_x_genes.genes.txt", emit: gene_list
     path "counts_unfiltered/cells_x_genes.mtx", emit: cell_gene_matrix


     script: 
     """
     kb count -i ${transcriptome_index} -g ${gene_map} -x ${chemistry} -m ${max_memory} -t ${threads} --parity ${parity} --strand ${strand} ${reads}
     """

}


workflow kb_python {

       main:

       fastqs = Channel.fromFilePairs( params.input_dir + '/' + addRecursiveSearch(params.recursive_search) + params.fastq_pattern )
       .ifEmpty { exit 1, "Cannot find any reads matching: ${params.input_dir} with recursive set to: ${params.recursive_search}\n" }

       if (! params.gene_map || ! params.transcriptome_index) {

       kb_ref(params.transcript_fasta, params.primary_fasta, params.gtf)

       GENE_MAP = kb_ref.out.gene_map
       INDEX = kb_ref.out.index

       } else if ( params.gene_map && params.transcriptome_index ) {

       GENE_MAP = params.gene_map
       INDEX = params.transcriptome_index

       

} else {
     println("Error: If a pre-generated gene map and transcriptome index are not supplied, \n then a transcript FASTA, primary FASTA, and reference gtf file need to be supplied.")
     System.exit(1)
}

      kb_count(fastqs, INDEX, GENE_MAP, params.chemistry, params.max_memory, params.threads, getKBPythonParity(params.chemistry, params.read_parity),
      params.strand)
}


workflow {

     main: 

      kb_python()

}



      
