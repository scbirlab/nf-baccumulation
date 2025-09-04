#!/usr/bin/env nextflow

/*
========================================================================================
   Bacterial compound accumulation screen Nextflow Workflow
========================================================================================
   Github   : https://github.com/scbirlab/nf-baccumulation
   Contact  : Eachan Johnson <eachan.johnson@crick.ac.uk>
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2
pipeline_title = """\
                 S C B I R   C O M P O U N D   A C C U M U L A T I O N   P I P E L I N E
                 =======================================================================
                 """
                 .stripIndent()

/*
========================================================================================
   Help text
========================================================================================
*/
if ( params.help ) {
   log.info pipeline_title + """\
         Nextflow pipeline to .... 

         Usage:
            nextflow run sbcirlab/nf-baccumulation --sample_sheet <csv> [--inputs <dir>]
            nextflow run sbcirlab/nf-baccumulation -c <config-file>

         Required parameters:
            sample_sheet      Path to a CSV with information about the samples 
                                 to be processed

         Optional parameters (with defaults):
            inputs            Directory containing inputs. Default: "./inputs".
            outputs           Directory to contain outputs. Default: "./outputs".

         The parameters can be provided either in the `nextflow.config` file or on the `nextflow run` command.
   
   """.stripIndent()
   exit 0
}

/*
========================================================================================
   Check parameters
========================================================================================
*/
if ( !params.sample_sheet ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a path to sample_sheet")
}

log.info pipeline_title + """\
   inputs
      sample sheet   : ${params.sample_sheet}
      input dir.     : ${params.inputs}
   output            : ${params.outputs}
   """
   .stripIndent()

/*
========================================================================================
   MAIN Workflow
========================================================================================
*/

include { 
   merge_compounds_and_sample_sheet;
} from './modules/hts-tools.nf'
include { 
   annotate_compounds;
} from './modules/schemist.nf'
include { 
   parse_targetlynx_to_table;
   targetlynx_table_summary;
} from './modules/targetlynx.nf'

workflow {

   Channel.fromPath( 
      params.sample_sheet, 
      checkIfExists: true,
   )
      .splitCsv( header: true )
      .set { csv_ch }

   Channel.of( params.sample_sheet )
      .map { file( it, checkIfExists: true) }
      .set { sample_sheet }

   csv_ch
      .map { tuple( 
         it.experiment_id, 
         file( 
            "${params.inputs}/${it.lcms_filename}", 
            checkIfExists: true,
         ) 
      ) }
      .set { lcms_ch }  // expt_id, LCMS files
   csv_ch
      .map { 
         tuple( 
            it.experiment_id, 
            file( 
               "${params.inputs}/${it.compound_info}", 
               checkIfExists: true,
            ) 
      ) }
      .set { compound_ch }  // expt_id, compounds

   /*
   ========================================================================================
      Processing
   ========================================================================================
   */

   lcms_ch
      .transpose() 
      | parse_targetlynx_to_table 
      | targetlynx_table_summary
   
   compound_ch
      .map { it[1] }
      .unique()
      .map { tuple( it.name, it ) } 
      | annotate_compounds

   targetlynx_table_summary.out
      .combine( compound_ch, by: 0 )  // expt_id, LCMS file, compound file
      .map { [ it[2].name ] + it[0..1] }  // compound filename, expt_id, LCMS file
      .combine( annotate_compounds.out, by: 0 )  // compound filename, expt_id, LCMS file, compounds
      .map { it[1..-1] }  // expt_id, LCMS file, compounds
      .combine( sample_sheet )  // expt_id, LCMS file, compounds, sample_sheet
      | merge_compounds_and_sample_sheet

   // ANNOTATE_DATA.out | (QC, PLOTS, SUMMARIZE)

}

/*
 * Merge sample sheet with data files.
 */
process QC {

   tag "${expt_id}"

   publishDir( qc_o, 
               mode: 'copy' )

   input:
   tuple val( expt_id ), path( data )

   output:
   tuple val( expt_id ), path( "*.tsv" ), path( "*.png" )

   script:
   """
   hts qc "${data}" \
      --control ${params.control_column} \
      --positive ${params.positive} \
      --negative ${params.negative} \
      --grouping ${params.grouping} \
      --plot "${expt_id}" \
      --output "${expt_id}_qc.tsv"
   """

}

/*
 * Summarize replicates with statistics.
 */
process SUMMARIZE {

   tag "${expt_id}"

   publishDir( hits_o, 
               mode: 'copy' )

   input:
   tuple val( expt_id ), path( data )

   output:
   tuple val( expt_id ), path( "*.tsv" ), path( "*.png" )

   script:
   """
   hts summarize "${data}" \
      --control ${params.control_column} \
      --positive ${params.positive} \
      --negative ${params.negative} \
      --grouping ${params.hit_grouping} \
      --plot "${expt_id}_summary" \
      --output "${expt_id}_summary.tsv" 
   """

}

/*
 * Visualize plate data.
 */
process PLOTS {

   tag "${expt_id}"

   publishDir( plots_o, 
               mode: 'copy' )

   input:
   tuple val( expt_id ), path( data )

   output:
   tuple val( expt_id ), path( "*.png" )

   script:
   """
   hts plot-hm "${data}" \
      --grouping ${params.grouping} \
      --output "${expt_id}" 

   hts plot-rep "${data}" \
      --control ${params.control_column} \
      --positive ${params.positive} \
      --negative ${params.negative} \
      --grouping ${params.hit_grouping} \
      --output "${expt_id}" 

   hts plot-hist "${data}" \
      --control "${params.control_column}" \
      --positive ${params.positive} \
      --negative ${params.negative} \
      --output "${expt_id}" 

   """

}

/*
========================================================================================
   Workflow Event Handler
========================================================================================
*/

workflow.onComplete {

   println ( workflow.success ? """
      Pipeline execution summary
      ---------------------------
      Completed at: ${workflow.complete}
      Duration    : ${workflow.duration}
      Success     : ${workflow.success}
      workDir     : ${workflow.workDir}
      exit status : ${workflow.exitStatus}
      """ : """
      Failed: ${workflow.errorReport}
      exit status : ${workflow.exitStatus}
      """
   )
}

/*
========================================================================================
   THE END
========================================================================================
*/