// Annotate compounds systematically
process annotate_compounds {

   tag "${id}" 

   publishDir( 
      "${params.outputs}/compounds", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   tuple val( id ), path( compound_info )

   output:
   tuple val( id ), path( "annotated.tsv" )

   script:
   """
   NLINES=\$(cat "${compound_info}" | wc -l)
   NLINES=\$((\$NLINES-1))
   if [ "\$NLINES" -lt 10000 ]
   then
      pubchem_extras="pubchem_id pubchem_name"
   else
      pubchem_extras=
   fi

   schemist convert "${compound_info}" \
      --to smiles inchikey id \$pubchem_extras scaffold mwt clogp tpsa \
      --options prefix=SCB- \
      --output annotated.tsv

   """

}
