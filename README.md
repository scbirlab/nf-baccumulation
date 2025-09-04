# Platereader screen analysis pipeline

![GitHub Workflow Status (with branch)](https://img.shields.io/github/actions/workflow/status/scbirlab/nf-baccumulation/nf-test.yml)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.10.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

**scbirlab/nf-baccumulation**  is a Nextflow pipeline process files from TargetLynx for intracellular accumulation.

**Table of contents**

- [Processing steps](#processing-steps)
- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [Issues, problems, suggestions](#issues-problems-suggestions)
- [Further help](#further-help)

## Processing steps

For each experiment in the `sample_sheet`:

### Other steps


## Requirements

### Software

You need to have Nextflow and either Anaconda, Singularity, or Docker installed on your system.

#### Crick users (and other HPC users)

If you're at the Crick or your shared cluster has it already installed, try:

```bash
module load Nextflow Singularity
```

#### Everyone else: installing Nextflow 

Otherwise, if it's your first time using Nextflow on your system and you have Conda installed, you can install it using `conda`:

```bash
conda install -c bioconda nextflow 
```

You may need to set the `NXF_HOME` environment variable. For example,

```bash
mkdir -p ~/.nextflow
export NXF_HOME=~/.nextflow
```

To make this a permanent change, you can do something like the following:

```bash
mkdir -p ~/.nextflow
echo "export NXF_HOME=~/.nextflow" >> ~/.bash_profile
source ~/.bash_profile
```

## Quick start

Make a [sample sheet (see below)](#sample-sheet) and, optionally, 
a [`nextflow.config` file](#inputs) in the directory where you want the 
pipeline to run. Then run Nextflow.

```bash 
nextflow run scbirlab/nf-baccumulation -latest
```

If you want to run a particular tagged version of the pipeline, such as `v0.0.1`, 
you can do so using

```bash 
nextflow run scbirlab/nf-baccumulation -r v0.0.1
```

For help, use `nextflow run scbirlab/nf-promotermap --help`.

The first time you run the pipeline, the software dependencies 
in `environment.yml` will be installed. This may take several minutes.

## Inputs

The following parameters are required:

- `sample_sheet`: path to a CSV with information about the samples and FASTQ files to be processed
- `fastq_dir`: path to where FASTQ files are stored
- `control_label`: the bin ID (from [sample sheet](#sample-sheet)) of background controls

The following parameters have default values which can be overridden if necessary.

- `inputs = "inputs"` : The folder containing your inputs.
- `outputs = "outputs"` : The folder to containing the pipeline outputs.

The parameters can be provided either in the `nextflow.config` file or on the `nextflow run` command.

Here is an example of the `nextflow.config` file:

```nextflow
params {
    sample_sheet = "/path/to/sample-sheet.csv"
    inputs = "/path/to/inputs"
}
```

Alternatively, you can provide the parameters on the command line:

```bash
nextflow run scbirlab/nf-promotermap \
    --sample_sheet /path/to/sample-sheet.csv \
    --inputs /path/to/inputs
``` 

### Sample sheet

The sample sheet is a CSV file providing information about which FASTQ files belong to which sample.

The file must have a header with the column names below (in any order), and one line per sample to be processed. 
You can have additional columns eith extra information if you like.

- `experiment_id`: Unique name of a peak-calling experiment. Peaks will be called across all samples with the same experiment ID.
- `lcms_filename`: Filename pattern of exported TargetLynx files.
- `compound_info`: Filename of compound library information.

Additional columns can be included. These will be included in output tables, so can be used for downstrwam analysis.

Here is an example of a couple of lines from a sample sheet:

| experiment_id | lcms_filename              | compound_info  |
| ------------- | -------------------------- | -------------- | 
| expt01        | TargetLynx/LCMS_*_Data.txt | compounds.xlsx |

### Example inputs

You cna find some examples in the `test` directory of this repository.

## Outputs

Outputs are saved in the directory specified by `--outputs` (`outputs` by default). 
They are organised into these directories:

## Issues, problems, suggestions

Add to the [issue tracker](https://www.github.com/scbirlab/nf-baccumulation/issues).

## Further help

Here are the help pages of the software used by this pipeline.

- [hts-tools](https://github.com/scbirlab/hts-tools)
- [Nextflow](https://www.nextflow.io/docs/latest/index.html)
- [schemist](https://github.com/scbirlab/schemist)