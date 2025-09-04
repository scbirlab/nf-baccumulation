// Convert Xevo text file to table
process parse_targetlynx_to_table {

    tag "${id}"

    publishDir( 
        "${params.outputs}/lcms", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( data )

    output:
    tuple val( id ), path( "raw.tsv" )

    script:
    """
    set -euox pipefail

    sed '1s/^\\xEF\\xBB\\xBF//' "${data}"\
    | tr -d \$'\\r' \
    | grep -v 'MM-' \
    | sed 's/bb//g;s/MM//g;s/db//g;s/bd//g;s/dd//g' \
    | awk -v OFS='\\t' -v source="${data}" '
        BEGIN { write_header=0 } \
        /^Compound [0-9]/ { 
            split(\$0, s, ":  "); 
            lcms_cmpd_id=s[1]; 
            cmpd_id=s[2]; 
            if (!write_header) write_header=1 
        } 
        (write_header==1 && \$1=="#") { 
            \$1=\$1; 
            print "lcms_filename", "lcms_compound_id", "compound_id", \$0; 
            write_header=2 
        } 
        match(\$1, /[0-9].+/) { 
            \$1=source OFS lcms_cmpd_id OFS cmpd_id; 
            print \$0 
        }
        ' \
    | python -c '
    import sys, numpy as np, pandas as pd
    (
        pd.read_csv(sys.stdin, sep="\\t")
        .assign(
            lcms_name=lambda x: np.where(x["Type"] == "Analyte", x["Name"].str.split("_").str[:-1].str.join("_"), x["Name"]),
            lcms_technical_rep=lambda x: np.where(x["Type"] == "Analyte", x["Name"].str.split("_").str[-1], np.nan),
        )
        .to_csv(sys.stdout, index=False, sep="\\t")
    )
    ' \
    > raw.tsv

    """

}


process targetlynx_table_summary {

    tag "${id}"

    publishDir( 
        "${params.outputs}/lcms", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( data )

    output:
    tuple val( id ), path( "summary.tsv" )

    script:
    """
    #!/usr/bin/env python

    import sys, numpy as np, pandas as pd

    cols = ["nM", "RT"]
    grouping = []

    g = pd.read_csv("${data}", sep="\\t").groupby(grouping)
    print(g, file=sys.stderr)
    summary = g[cols].agg(["mean", "std"])

    summary.columns = ["_".join(a) for a in summary.columns.to_flat_index()]
    summary.reset_index().to_csv("summary.tsv", sep="\\t", index=False)

    """

}
