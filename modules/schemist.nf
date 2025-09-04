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
   clean_filename=clean.${compound_info.extension}
   sed '1s/^\\xEF\\xBB\\xBF//' "${compound_info}" | tr -d \$'\\r' > "\$clean_filename"
   NLINES=\$(cat "\$clean_filename" | wc -l)
   NLINES=\$((\$NLINES-1))
   pubchem_extras=
   #if [ "\$NLINES" -lt 10000 ]
   #then
   #   pubchem_extras="pubchem_id pubchem_name"
   #else
   #   pubchem_extras=
   #fi

   schemist convert "\$clean_filename" \
      --to smiles inchikey id \$pubchem_extras scaffold mwt clogp tpsa \
      --options prefix=SCB- \
   | awk -v OFS='\\t' -v fname="${compound_info}" '
      NR == 1 { print "compound_info", \$0 }
      NR > 0 { print fname, \$0 }
      ' \
   > annotated.tsv

   """

}
