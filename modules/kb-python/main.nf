#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process kb_ref {

     module = 'python3'
   
     publishDir path: "${params.output_dir}/ref/", mode: "copy"

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

     publishDir path: "${params.output_dir}/kallisto_count/${name}/", mode: "copy"
     
     memory '4 GB'
   
     input:
     tuple val(name), path(reads)
     path transcriptome_index
     path gene_map


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
     kb count -i ${transcriptome_index} -g ${gene_map} -x ${params.chemistry} -m ${params.max_memory} -t ${params.threads} ${reads} 
     """

}


workflow kb_workflow {

       main:

       fastqs = Channel.fromFilePairs( params.input_pattern )
       .ifEmpty { exit 1, "Cannot find any reads matching: ${params.input_pattern}\n" }

       if (!params.gene_map || !params.transcriptome_index) {

       kb_ref(params.transcript_fasta, params.primary_fasta, params.gtf)

       GENE_MAP = kb_ref.out.gene_map
       INDEX = kb_ref.out.index

       } else {

       GENE_MAP = params.gene_map
       INDEX = params.transcriptome_index

       

}
      kb_count(fastqs, INDEX, GENE_MAP)
}


workflow {

     main: 

      kb_workflow()

}



      
