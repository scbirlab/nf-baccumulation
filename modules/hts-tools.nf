process merge_compounds_and_sample_sheet {

   tag "${id}" 

   publishDir( 
      "${params.outputs}/lcms", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   tuple val( id ), path( data ), path( compound_info ), path( sample_sheet )

   output:
   tuple val( id ), path( "merged.tsv" )

   script:
   """
   hts join ${data} --right "${sample_sheet}" --format-right csv > first.tsv
   hts join first.tsv --right "${compound_info}" \
   > merged.tsv

   """

}

process normalize {

   tag "${id}:pos=${positive}:neg=${negative}:grouping=${grouping}" 

   publishDir( 
      "${params.outputs}/processed", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   tuple val( id ), path( data )
   val control_name
   val pos
   val neg
   val method

   output:
   tuple val( id ), path( "norm.tsv" )

   script:
   """
   hts normalize "${data}" \
      --control "${control_name}" \
      --positive "${pos}" \
      --negative "${neg}" \
      --grouping experiment_id \
      --method "${method}" \
   > norm.tsv

   """

}


process qc {

   tag "${id}:pos=${positive}:neg=${negative}:grouping=${grouping}" 

   publishDir( 
      "${params.outputs}/qc", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   tuple val( id ), path( data )
   val control_name
   val positive
   val negative
   val grouping

   output:
   tuple val( id ), path( "qc.tsv" ), emit: main
   tuple val( id ), path( "*.png" ), emit: plots

   script:
   """
   hts qc "${data}" \
      --control ${control_name} \
      --positive ${positive} \
      --negative ${negative} \
      --grouping experiment_id \
      --plot "qc" \
      --output "qc.tsv"

   """

}

/*
 * Summarize replicates with statistics.
 */
process summarize {

   ttag "${id}:pos=${positive}:neg=${negative}:grouping=${grouping}" 

   publishDir( 
      "${params.outputs}/summaries", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   tuple val( id ), path( data )
   val control_name
   val positive
   val negative
   val grouping

   output:
   tuple val( id ), path( "summary.tsv" ), emit: main
   tuple val( id ), path( "*.png" ), emit: plots


   script:
   """
   hts summarize "${data}" \
      --control ${control_name} \
      --positive ${positive} \
      --negative ${negative} \
      --grouping ${grouping} \
      --plot "summary" \
      --output "summary.tsv"
   
   """

}

/*
 * Visualize plate data.
 */
process plots {

   tag "${id}:pos=${positive}:neg=${negative}:grouping=${grouping}" 

   publishDir( 
      "${params.outputs}/plots", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   tuple val( id ), path( data )
   val control_name
   val positive
   val negative
   val grouping

   output:
   tuple val( id ), path( "*.png" )

   script:
   """
   #hts plot-hm "${data}" \
   #   --grouping ${grouping} \
   #   --output "${id}" 

   hts plot-rep "${data}" \
      --control ${control_name} \
      --positive ${positive} \
      --negative ${negative} \
      --grouping ${grouping} \
      --output "reps" 

   hts plot-hist "${data}" \
      --control ${control_name} \
      --positive ${positive} \
      --negative ${negative} \
      --output "hist" 

   """

}
