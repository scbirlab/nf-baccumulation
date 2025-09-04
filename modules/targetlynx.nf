// Convert Xevo text file to table
process parse_targetlynx_to_table {

    tag "${id}"

    publishDir( 
        "${params.outputs}/lcms-raw", 
        mode: 'copy',
        saveAs: { "${id}.${it}" },
    )

    input:
    tuple val( id ), path( data )

    output:
    tuple val( id ), path( "${data.baseName}.raw.tsv" )

    script:
    """
    set -euox pipefail

    sed '1s/^\\xEF\\xBB\\xBF//' "${data}" \
    | tr -d \$'\\r' \
    | grep -v 'MM-' \
    | sed 's/bb//g;s/MM//g;s/db//g;s/bd//g;s/dd//g' \
    | awk -v OFS='\\t' -v source="${data}" -v experiment_id="${id}" '
        BEGIN { write_header=0 } \
        /^Compound [0-9]/ { 
            split(\$0, s, ":  "); 
            lcms_cmpd_id=s[1]; 
            cmpd_id=s[2]; 
            if (!write_header) write_header=1 
        } 
        (write_header == 1 && \$1 == "#") { 
            \$1=\$1; 
            print "experiment_id", "lcms_source_file", "lcms_compound_id", "compound_id", \$0; 
            write_header=2 
        } 
        match(\$1, /[0-9].+/) { 
            \$1=experiment_id OFS source OFS lcms_cmpd_id OFS cmpd_id;
            print \$0 
        }
        ' \
    | python -c '
    import sys, numpy as np, pandas as pd
    (
        pd.read_csv(sys.stdin, sep="\\t")
        .assign(
            lcms_name=lambda x: np.where(x["Type"] == "Analyte", x["Name"].str.split("_").str[:-1].str.join("_"), x["Name"]),
            lcms_bio_rep=lambda x: np.where(x["Type"] == "Analyte", x["Name"].str.split("_").str[-2], np.nan),
            lcms_technical_rep=lambda x: np.where(x["Type"] == "Analyte", x["Name"].str.split("_").str[-1], np.nan),
            lcms_qc_conc=lambda x: np.where(x["Type"] == "QC", x["Name"].str.split("_").str[-1].str.replace("nM", "").astype(float), np.nan),
        )
        .rename(columns={
            "nM": "concentration",
            "RT": "retention_time",
        })
        .assign(
            concentration_unit="nM", 
            retention_time_unit="min",
            measured_concentration_ch1=lambda x: x["concentration"],
            measured_retention_ch1=lambda x: x["retention_time"],
            concentration_ch1_wavelength="0nm", 
            retention_ch1_wavelength="0nm", 
        )
        .to_csv(sys.stdout, index=False, sep="\\t")
    )
    ' \
    > "${data.baseName}.raw.tsv"

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
    tuple val( id ), path( data, stageAs: '*.raw.tsv' )

    output:
    tuple val( id ), path( "summary-analyte.tsv" ), emit: analyte
    tuple val( id ), path( "summary-qc.tsv" ), emit: qc
    tuple val( id ), path( "summary-analyte-tech.tsv" ), emit: tech_analyte
    tuple val( id ), path( "summary-qc-tech.tsv" ), emit: tech_qc

    script:
    """
    #!/usr/bin/env python

    import glob, sys, numpy as np, pandas as pd

    cols = ["concentration", "retention_time"]
    grouping = [
        "experiment_id",
        "lcms_source_file",
        "lcms_compound_id",
        "compound_id", 
        "lcms_qc_conc", 
        "Type",
        "concentration_unit",
        "retention_time_unit",
        "lcms_name", 
        "lcms_bio_rep",
    ]

    df = pd.concat([
        pd.read_csv(f, sep="\\t")
        for f in glob.glob("*.raw.tsv")
    ], axis=0)

    for groups, suffix in ((grouping, "-tech"), (grouping[:-2], "")):
        summary = (
            df
            .groupby(groups, dropna=False)
            [cols]
            .agg(["mean", "std", "count"])
        )
        summary.columns = [
            "_".join(a) for a in summary.columns.to_flat_index()
        ]
        if suffix == "-tech":
            summary = summary.assign(
                measured_concentration_ch1=lambda x: x["concentration_mean"],
                measured_retention_ch1=lambda x: x["retention_time_mean"],
                concentration_ch1_wavelength="0nm", 
                retention_ch1_wavelength="0nm", 
            )
        for query in ("Analyte", "QC"):
            (
                summary
                .query("Type == @query")
                .reset_index()
                .to_csv(f"summary-{query.casefold()}{suffix}.tsv", sep="\\t", index=False)
            )

    """

}
