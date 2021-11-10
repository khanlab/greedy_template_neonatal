# Snakemake workflow: `greedy_template_neonatal`

Workflow to build cohort templates based on multiple participants_tsv files

This workflow uses greedy instead of ANTS, for the sake of efficiency, in fact, the registrations as compared to `ants_build_template` seem to be more accurate (though that of course is likely just due to parameter selection, not an inherent limitation of ANTS).. There is ~20-30x speedup for a single pairwise registration compared to ANTS, making template-generation for hundreds of subjects a job that can be completed in around a day with modest resources (32cores). The 4-core 16gb memory greedy registration jobs take <30mins.


## Running with a cluster profile (e.g. cc-slurm)
This option will be the most # of jobs (`num_iters * num_subjects * num_cohorts`), but will maximize parallelization on a cluster
```
snakemake  --profile cc-slurm 
```

## Running with cluster profile, but grouping registration jobs together in chunks of 8
This will reduce the number of jobs by factor of 8 by grouping the 8 4core jobs to fill the 32-core nodes

HOWEVER -- if you are using the --group-components option to group the registrations in chunks of 8, you *MUST* run each iteration separately, using the `--config run_iter=#` option, where # is the iteration to run -- this seems to be a bug/limitation of the group-components option when a job has recursively-defined rules. (open issue as of Oct 3 2020: https://github.com/snakemake/snakemake/issues/656)
```
snakemake --config run_iter=2 --profile cc-slurm --group-components reg=8 composite=100
```
Note: composite=100 is for the composing warps to another ref space (e.g. mni) to combine 100 subjects in a single job (each is relatively quick and disk-bound)..


## Running as a single job
This is the most frugal for resources (32-cores max), but is simple in that only a single job is spawned.. This is the method to use if you are running this on a single machine. Four iterations of building a single 1mm template with 100 subjects takes < 24hrs on a 32core system.

regularSubmit -j Fat snakemake --use-singularity  -j32

## Running each cohort in parallel with single jobs and the --nolock flag
This runs each cohort separately with single-node jobs.. This is a happy medium for running on a cluster while not needing hundreds (or thousands) of short jobs.. 

WARNING: make sure your cohorts are mutually-exclusive when using this method, as it is running snakemake in parallel on the same directory with the --nolock option -- if cohorts are not mutually-exclusive, you can still use this method, but only after all the pre-processing (e.g. T2/T1 registration, is completed)
```
for cohort in young middle1 middle2 old1 old2; do  regularSubmit -j Fat snakemake --use-singularity -j32 --nolock --config run_cohort=$cohort; done
```

