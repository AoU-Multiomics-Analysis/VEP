# VEP Pipeline with LOFTEE and dbNSFP

This WDL pipeline runs VEP (Variant Effect Predictor) with LOFTEE and dbNSFP plugins for comprehensive variant annotation.

## Overview

This pipeline combines:
- VEP annotation workflow structure from [evin-padhi/VEP](https://github.com/evin-padhi/VEP/blob/main/VEP.wdl)
- LOFTEE plugin integration from [sph17/loftee_plugin](https://github.com/sph17/loftee_plugin/blob/v0.2.1/wdl/vepAnnotate.wdl)
- dbNSFP plugin for functional annotations

## Features

- **LOFTEE Plugin**: Identifies loss-of-function variants with high confidence
- **dbNSFP Plugin**: Provides comprehensive functional annotations from dbNSFP database
- **Chromosome-level parallelization**: Splits VCF by chromosome for efficient processing
- **VCF normalization**: Normalizes variants before annotation
- **Genotype preservation**: Restores genotypes after annotation

## Workflow Steps

1. **Normalize VCF**: Normalizes variants using bcftools
2. **Split by Chromosome**: Splits VCF into per-chromosome shards
3. **VEP Annotation**: Runs VEP with LOFTEE and dbNSFP plugins on each chromosome
4. **Merge VCFs**: Concatenates chromosome-level annotated VCFs
5. **Add Genotypes**: Merges genotypes back into the annotated VCF

## Required Inputs

### Input VCF
- `vcf_file`: Input VCF file (File)
- `contigs`: Array of chromosome names to process (Array[String])

### Docker Images
- `vep_docker`: Docker image with VEP, LOFTEE, and dbNSFP plugins installed
- `sv_base_mini_docker`: Docker image with bcftools for VCF processing

### Reference Files
- `hg38_fasta`: Human reference genome (GRCh38) FASTA file
- `hg38_fasta_fai`: Index for the reference FASTA
- `top_level_fa`: Top-level FASTA file for VEP

### LOFTEE Plugin Files
- `human_ancestor_fa`: Human ancestor FASTA file
- `human_ancestor_fa_fai`: Index for human ancestor FASTA
- `gerp_conservation_scores`: GERP conservation scores (BigWig file)

### dbNSFP Plugin Files
- `dbnsfp_database`: dbNSFP database file (typically dbNSFP*.txt.gz)
- `dbnsfp_database_tbi`: Tabix index for dbNSFP database
- `dbnsfp_fields`: Array of dbNSFP fields to annotate (e.g., ["SIFT_score", "Polyphen2_HDIV_score"])

### Other Parameters
- `cohort_prefix`: Prefix for output files
- `vep_assembly`: Genome assembly version (default: "GRCh38")

## Example Input JSON

```json
{
  "VepWithLofteeAndDbNSFP.vcf_file": "gs://my-bucket/input.vcf.gz",
  "VepWithLofteeAndDbNSFP.contigs": ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"],
  "VepWithLofteeAndDbNSFP.vep_docker": "my-registry/vep-loftee-dbnsfp:latest",
  "VepWithLofteeAndDbNSFP.sv_base_mini_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:latest",
  "VepWithLofteeAndDbNSFP.hg38_fasta": "gs://my-bucket/references/Homo_sapiens_assembly38.fasta",
  "VepWithLofteeAndDbNSFP.hg38_fasta_fai": "gs://my-bucket/references/Homo_sapiens_assembly38.fasta.fai",
  "VepWithLofteeAndDbNSFP.top_level_fa": "gs://my-bucket/references/Homo_sapiens.GRCh38.dna.toplevel.fa.gz",
  "VepWithLofteeAndDbNSFP.human_ancestor_fa": "gs://my-bucket/loftee/human_ancestor.fa.gz",
  "VepWithLofteeAndDbNSFP.human_ancestor_fa_fai": "gs://my-bucket/loftee/human_ancestor.fa.gz.fai",
  "VepWithLofteeAndDbNSFP.gerp_conservation_scores": "gs://my-bucket/loftee/gerp_conservation_scores.homo_sapiens.GRCh38.bw",
  "VepWithLofteeAndDbNSFP.dbnsfp_database": "gs://my-bucket/dbnsfp/dbNSFP4.3a_grch38.gz",
  "VepWithLofteeAndDbNSFP.dbnsfp_database_tbi": "gs://my-bucket/dbnsfp/dbNSFP4.3a_grch38.gz.tbi",
  "VepWithLofteeAndDbNSFP.dbnsfp_fields": ["SIFT_score", "SIFT_pred", "Polyphen2_HDIV_score", "Polyphen2_HDIV_pred", "MutationTaster_score", "MutationTaster_pred"],
  "VepWithLofteeAndDbNSFP.cohort_prefix": "my_cohort",
  "VepWithLofteeAndDbNSFP.vep_assembly": "GRCh38"
}
```

## Outputs

- `vep_annotated_final_vcf`: Final annotated VCF with genotypes restored
- `vep_annotated_final_vcf_idx`: Tabix index for the final VCF

## Docker Image Requirements

The VEP docker image must have:
- VEP installed and configured
- LOFTEE plugin installed at `/opt/vep/.vep/Plugins/`
- dbNSFP plugin support (typically included with VEP)
- VEP cache at `/opt/vep/.vep/`

## Reference File Sources

### LOFTEE Files
- **human_ancestor_fa**: Download from [LOFTEE resources](https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/)
- **gerp_conservation_scores**: Download from [LOFTEE resources](https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/)

### dbNSFP Database
- Download from [dbNSFP database](https://sites.google.com/site/jpopgen/dbNSFP)
- Ensure proper indexing with tabix

### VEP Cache
- Download from [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html)

## Runtime Customization

You can override runtime attributes for each task using `RuntimeAttr` inputs:
- `runtime_attr_normalize`
- `runtime_attr_split_vcf`
- `runtime_attr_vep_annotate`
- `runtime_attr_merge_vcfs`
- `runtime_attr_add_genotypes`

Example:
```json
{
  "VepWithLofteeAndDbNSFP.runtime_attr_vep_annotate": {
    "mem_gb": 8,
    "cpu_cores": 2,
    "disk_gb": 100,
    "preemptible_tries": 2
  }
}
```

## Notes

- LOFTEE paths are pre-configured in the VEP command and assume standard plugin installation at `/opt/vep/.vep/Plugins/`
- The pipeline processes chromosomes in parallel for better performance
- If `dbnsfp_fields` is empty, the dbNSFP plugin will not be used
- The workflow preserves the original genotype information through the normalization and annotation process