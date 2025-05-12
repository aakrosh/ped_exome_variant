This workflow takes as input CRAM files from related individuals that includes at least one proband, and then performs the following tasks:
1. Realigns the reads to a new reference using bwa mem. 
2. Marks duplicates using SAMblaster. 
3. Calls SNV and indels using FreeBayes.
4. Runs Exomiser to prioritise the variants.

## Running the workflow
```{bash}
nextflow -c nextflow.config -profile uva -params-file params.yml run aakrosh/ped_exome_variant -resume 
```

## Notes
This workflow is not ready for public consumption.
