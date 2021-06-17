# Workflow

This folder contains a Nextflow script and configuration, intended to run on
the UPPMAX clusters.

```
workflow/
 | - bin/                            A folder for custom workflow scripts
 | - configs/                        Configuration files that govern workflow execution
 | - containers/                     Custom container definition files
 | - main.nf                         The primary analysis script
 | - modules/                        Nextflow process definitions (modules)
 | - nextflow.config                 General Nextflow configuration
 \ - params.config.TEMPLATE          A Template for parameter configuration
```

There can be more than one workflow script if one desires. The Nextflow DSL2 syntax also
supports modules / subworkflows, and it can be useful to make a folder for those.  

## Customisation for Uppmax

A custom profile named `uppmax` is available to run this workflow specifically
on UPPMAX clusters. The process `executor` is `slurm` so jobs are
submitted to the Slurm Queue Manager. All jobs submitted to slurm
must have a project allocation. This is automatically added to the `clusterOptions`
in the `uppmax` profile. All Uppmax clusters have node local disk space to do
computations, and prevent heavy input/output over the network (which
slows down the cluster for all).
The path to this disk space is provided by the `$SNIC_TMP` variable, used by
the `process.scratch` directive in the `uppmax` profile. Lastly
the profile enables the use of Singularity so that all processes must be
executed within Singularity containers. See [nextflow.config](nextflow.config)
for the profile specification.

The profile is enabled using the `-profile` parameter to nextflow:
```bash
nextflow run -profile uppmax <nextflow_script>
```
