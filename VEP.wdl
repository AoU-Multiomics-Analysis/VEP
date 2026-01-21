version 1.0

import "Structs.wdl"

workflow VepWithLofteeAndDbNSFP {

    input {
        File vcf_file
        Array[String] contigs
        
        # Docker images
        String vep_docker
        String sv_base_mini_docker
        
        # Reference files
        File hg38_fasta
        File hg38_fasta_fai
        
        # LOFTEE plugin files
        File human_ancestor_fa
        File human_ancestor_fa_fai
        File gerp_conservation_scores
        
        # dbNSFP plugin files
        File dbnsfp_database
        File dbnsfp_database_tbi
        
        # Additional files
        File top_level_fa
        
        # VEP cache
        File vep_cache_tar_gz
        
        # VEP options
        String vep_assembly = "GRCh38"
        Array[String] dbnsfp_fields = []
        
        String cohort_prefix
        
        # Runtime attributes
        RuntimeAttr? runtime_attr_normalize
        RuntimeAttr? runtime_attr_split_vcf
        RuntimeAttr? runtime_attr_vep_annotate
        RuntimeAttr? runtime_attr_merge_vcfs
        RuntimeAttr? runtime_attr_add_genotypes
    }

    # Normalize VCF file
    call normalizeVCF {
        input:
            vcf_file = vcf_file,
            sv_base_mini_docker = sv_base_mini_docker,
            hg38_fasta = hg38_fasta,
            hg38_fasta_fai = hg38_fasta_fai,
            runtime_attr_override = runtime_attr_normalize
    }

    # Split by chromosome, VEP annotate with LOFTEE and dbNSFP, and merge
    scatter (contig in contigs) {
        call splitVCF {
            input:
                vcf_file = normalizeVCF.vcf_no_genotype,
                vcf_idx_file = normalizeVCF.vcf_no_genotype_idx,
                chromosome = contig,
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_split_vcf
        }

        call vepAnnotate {
            input:
                vcf_file = splitVCF.vcf_output,
                top_level_fa = top_level_fa,
                human_ancestor_fa = human_ancestor_fa,
                human_ancestor_fa_fai = human_ancestor_fa_fai,
                gerp_conservation_scores = gerp_conservation_scores,
                dbnsfp_database = dbnsfp_database,
                dbnsfp_database_tbi = dbnsfp_database_tbi,
                dbnsfp_fields = dbnsfp_fields,
                vep_assembly = vep_assembly,
                vep_cache_tar_gz = vep_cache_tar_gz,
                vep_docker = vep_docker,
                runtime_attr_override = runtime_attr_vep_annotate
        }
    }

    call mergeVCFs {
        input:
            vcf_contigs = vepAnnotate.vep_vcf_file,
            sv_base_mini_docker = sv_base_mini_docker,
            cohort_prefix = cohort_prefix,
            runtime_attr_override = runtime_attr_merge_vcfs
    }

    call addGenotypes {
        input:
            vep_annotated_vcf = mergeVCFs.merged_vcf_file,
            normalized_vcf = normalizeVCF.vcf_normalized_file_with_genotype,
            normalized_vcf_idx = normalizeVCF.vcf_normalized_file_with_genotype_idx,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_add_genotypes
    }

    output {
        File vep_annotated_final_vcf = addGenotypes.combined_vcf_file
        File vep_annotated_final_vcf_idx = addGenotypes.combined_vcf_idx
    }
}

task normalizeVCF {
    input {
        File vcf_file
        String sv_base_mini_docker
        File hg38_fasta
        File hg38_fasta_fai
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, "GB")
    Float base_disk_gb = 10.0
    Float base_mem_gb = 2.0
    Float input_mem_scale = 3.0
    Float input_disk_scale = 5.0
    
    RuntimeAttr runtime_default = object {
        mem_gb: base_mem_gb + input_size * input_mem_scale,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String prefix = basename(vcf_file, ".vcf.gz")
    String vcf_normalized_file_name = "~{prefix}.normalized.vcf.gz"
    String vcf_normalized_nogeno_file_name = "~{prefix}.normalized.stripped.vcf.gz"

    output {
        File vcf_normalized_file_with_genotype = "~{prefix}.normalized.vcf.gz"
        File vcf_normalized_file_with_genotype_idx = "~{prefix}.normalized.vcf.gz.tbi"
        File vcf_no_genotype = "~{prefix}.normalized.stripped.vcf.gz"
        File vcf_no_genotype_idx = "~{prefix}.normalized.stripped.vcf.gz.tbi"
    }

    command <<<
        set -euo pipefail

        bcftools norm -m - -f ~{hg38_fasta} -Oz -o ~{vcf_normalized_file_name} ~{vcf_file}

        bcftools index -t ~{vcf_normalized_file_name}

        bcftools view -G ~{vcf_normalized_file_name} -Oz -o ~{vcf_normalized_nogeno_file_name}

        bcftools index -t ~{vcf_normalized_nogeno_file_name}
    >>>
}

task splitVCF {
    input {
        File vcf_file
        File vcf_idx_file
        String chromosome
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, "GB")
    Float base_disk_gb = 10.0
    Float base_mem_gb = 2.0
    Float input_mem_scale = 3.0
    Float input_disk_scale = 5.0
    
    RuntimeAttr runtime_default = object {
        mem_gb: base_mem_gb + input_size * input_mem_scale,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    output {
        File vcf_output = "~{chromosome}.vcf.gz"
    }

    command <<<
        set -euo pipefail  

        bcftools view ~{vcf_file} --regions ~{chromosome} -Oz -o ~{chromosome}.vcf.gz
    >>>
}

task vepAnnotate {
    input {
        File vcf_file
        File top_level_fa
        File human_ancestor_fa
        File human_ancestor_fa_fai
        File gerp_conservation_scores
        File dbnsfp_database
        File dbnsfp_database_tbi
        Array[String] dbnsfp_fields
        String vep_assembly
        File vep_cache_tar_gz
        String vep_docker
        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(vcf_file, ".vcf.gz")
    String vep_annotated_vcf_name = "~{prefix}.vep.loftee.dbnsfp.vcf.gz"

    Float input_size = size(vcf_file, "GB")
    Float ref_size = size([top_level_fa, human_ancestor_fa, gerp_conservation_scores, dbnsfp_database, vep_cache_tar_gz], "GB")
    Float base_disk_gb = 10.0
    Float base_mem_gb = 2.0
    Float input_mem_scale = 3.0
    Float input_disk_scale = 5.0
    
    RuntimeAttr runtime_default = object {
        mem_gb: base_mem_gb + input_size * input_mem_scale,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale + ref_size * 2.0),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: vep_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    output {
        File vep_vcf_file = vep_annotated_vcf_name
    }

    command <<<
        set -euo pipefail

        # Decompress VEP cache
        echo "Decompressing VEP cache..."
        mkdir -p vep_cache
        tar -xzf ~{vep_cache_tar_gz} -C vep_cache
        
        # Find the actual cache directory (it may be nested)
        HOMO_SAPIENS_DIR=$(find vep_cache -type d -name "homo_sapiens" | head -n 1)
        if [ -n "$HOMO_SAPIENS_DIR" ]; then
            VEP_CACHE_DIR=$(dirname "$HOMO_SAPIENS_DIR")
        else
            # If homo_sapiens directory not found, use the extracted directory
            VEP_CACHE_DIR="vep_cache"
        fi
        
        echo "Using VEP cache directory: $VEP_CACHE_DIR"
        
        # Look for synonyms file in the cache
        SYNONYMS_FILE=$(find "$VEP_CACHE_DIR" -type f \( -name "*synonym*" -o -name "*synonyms*" \) | head -n 1)
        
        # Move dbNSFP database files to current directory for easier access
        mv ~{dbnsfp_database} .
        mv ~{dbnsfp_database_tbi} .
        
        dbnsfp_basename=$(basename ~{dbnsfp_database})

        # Build dbNSFP plugin command if fields are provided
        dbnsfp_plugin=""
        if [ ~{length(dbnsfp_fields)} -gt 0 ]; then
            dbnsfp_plugin="--plugin dbNSFP,$dbnsfp_basename,~{sep=',' dbnsfp_fields}"
        fi
        
        # Build VEP command with optional synonyms
        vep_cmd="vep --vcf \
            --force_overwrite \
            --dir \"$VEP_CACHE_DIR\" \
            --format vcf \
            --everything \
            --allele_number \
            --no_stats \
            --cache \
            --offline \
            --minimal \
            --assembly ~{vep_assembly} \
            --fasta ~{top_level_fa} \
            --input_file ~{vcf_file} \
            --output_file ~{vep_annotated_vcf_name} \
            --compress_output bgzip \
            --plugin LoF,loftee_path:/opt/vep/.vep/Plugins/,human_ancestor_fa:~{human_ancestor_fa},gerp_score:~{gerp_conservation_scores} \
            --dir_plugins /opt/vep/.vep/Plugins/"
        
        # Add synonyms file if found
        if [ -n "$SYNONYMS_FILE" ]; then
            echo "Found synonyms file: $SYNONYMS_FILE"
            vep_cmd="$vep_cmd --synonyms \"$SYNONYMS_FILE\""
        fi
        
        # Add dbNSFP plugin if configured
        if [ -n "$dbnsfp_plugin" ]; then
            vep_cmd="$vep_cmd $dbnsfp_plugin"
        fi
        
        # Execute VEP command
        eval "$vep_cmd"
    >>>
}

task mergeVCFs {
    input {
        Array[File] vcf_contigs
        String sv_base_mini_docker
        String cohort_prefix
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_contigs, "GB")
    Float base_disk_gb = 10.0
    Float base_mem_gb = 2.0
    Float input_mem_scale = 3.0
    Float input_disk_scale = 5.0
    
    RuntimeAttr runtime_default = object {
        mem_gb: base_mem_gb + input_size * input_mem_scale,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String merged_vcf_name = "~{cohort_prefix}.vep.merged.vcf.gz"

    output {
        File merged_vcf_file = merged_vcf_name
    }

    command <<<
        set -euo pipefail
        
        VCFS="~{write_lines(vcf_contigs)}"
        cat $VCFS | awk -F '/' '{print $NF"\t"$0}' | sort -k1,1V | awk '{print $2}' > vcfs_sorted.list
        bcftools concat --no-version --naive -Oz --file-list vcfs_sorted.list --output ~{merged_vcf_name}
    >>>
}

task addGenotypes {
    input {
        File vep_annotated_vcf
        File normalized_vcf
        File normalized_vcf_idx
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
        Int? thread_num_override
    }

    Float vep_annotate_sizes = size(vep_annotated_vcf, "GB") 
    Float norm_vcf_sizes = size(normalized_vcf, "GB")
    Float base_disk_gb = 10.0

    RuntimeAttr runtime_default = object {
        mem_gb: 16,
        disk_gb: ceil(base_disk_gb + (vep_annotate_sizes + norm_vcf_sizes) * 5.0),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }
    
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String prefix = basename(vep_annotated_vcf, ".vcf.gz")
    String combined_vcf_name = "~{prefix}.vep.geno.vcf.gz"
    Int thread_num = select_first([thread_num_override, 1])

    output {
        File combined_vcf_file = combined_vcf_name
        File combined_vcf_idx = combined_vcf_name + ".tbi"
    }

    command <<<
        set -euo pipefail
        
        bcftools index -t ~{vep_annotated_vcf}

        bcftools merge \
            --no-version \
            --threads ~{thread_num} \
            -Oz \
            --output ~{combined_vcf_name} \
            ~{normalized_vcf} \
            ~{vep_annotated_vcf}

        bcftools index -t ~{combined_vcf_name}
    >>>
}
