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
   tuple val( id ), path( "annotated.tsv" )

   script:
   """
   hts join ${data} --right "${sample_sheet}" \
   | hts join --right "${compound_info}" \
   > annotated.tsv

   """

}