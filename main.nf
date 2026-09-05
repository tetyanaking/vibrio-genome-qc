#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


include { RESOLVE_SRA }             from './modules/local/resolve_sra'
include { DOWNLOAD_READS }          from './modules/local/download_reads'
include { MERGE_RUNS }              from './modules/local/merge_runs'
include { FETCH_ORIGINAL_ASSEMBLY } from './modules/local/fetch_original_assembly'

include { FASTP }                   from './modules/local/fastp'
include { KRAKEN2 }                 from './modules/local/kraken2'
include { EXTRACT_VIBRIO_READS }    from './modules/local/extract_vibrio_reads'
include { UNICYCLER }               from './modules/local/unicycler'
include { COV_LENGTH_FILTER }       from './modules/local/cov_length_filter'

include { GUNC }                    from './modules/local/gunc'
include { APPLY_FCS_POLICY }        from './modules/local/apply_fcs_policy'

include { CHECKM2 as CHECKM2_ORIGINAL } from './modules/local/checkm2'
include { CHECKM2 as CHECKM2_PRE }      from './modules/local/checkm2'
include { CHECKM2 as CHECKM2_POST }     from './modules/local/checkm2'
include { CHECKM2_DELTA }               from './modules/local/checkm2_delta'

include { BUSCO as BUSCO_ORIGINAL } from './modules/local/busco'
include { BUSCO as BUSCO_POST }     from './modules/local/busco'

include { SEQKIT_STATS as STATS_ORIGINAL } from './modules/local/seqkit_stats'
include { SEQKIT_STATS as STATS_PRE }      from './modules/local/seqkit_stats'
include { SEQKIT_STATS as STATS_POST }     from './modules/local/seqkit_stats'

include { SAMPLE_SUMMARY }          from './modules/local/summary'
include { COHORT_SUMMARY }          from './modules/local/cohort_summary'

include { MULTIQC }                 from './modules/local/multiqc'


workflow {

    /* Input metadata manifest */
    def input_file = params.input ?: params.samplesheet
    ch_metadata = Channel
        .fromPath(input_file, checkIfExists: true)
        .splitCsv(header: true)
        .filter { row ->
            row.sample_id?.trim() &&
            row.assembly_accession?.trim() &&
            row.biosample?.trim()
        }
        .map { row ->
            tuple(
                row.sample_id.trim(),
                row.assembly_accession.trim(),
                row.biosample.trim()
            )
        }

    /* Separate channels for processes that consume the metadata*/
    ch_sra_metadata = ch_metadata.map {
        sample_id, assembly_accession, biosample ->

        tuple(sample_id, assembly_accession, biosample)
    }

    ch_assembly_metadata = ch_metadata.map {
        sample_id, assembly_accession, biosample ->

        tuple(sample_id, assembly_accession)
    }

    /*1. Resolve BioSample accessions to SRR runs*/
    RESOLVE_SRA(ch_sra_metadata)
    ch_srr = RESOLVE_SRA.out.runs
        .flatMap { sample_id, assembly_accession, run_table ->

            run_table
                .readLines()
                .drop(1)
                .findAll { it.trim() }
                .collect { line ->

                    def fields = line.split('\t')

                    tuple(sample_id, fields[1])
                }
        }

    /*2. Download reads*/
    DOWNLOAD_READS(ch_srr)
    ch_runs = DOWNLOAD_READS.out.reads
     .map { sample_id, srr, reads ->
        tuple(sample_id, reads)
     }
     .groupTuple()
     .map { sample_id, grouped_reads ->
        tuple(sample_id, grouped_reads.flatten())
     }

    /* 3. Merge multiple SRR runs belonging to one isolate*/
    MERGE_RUNS(ch_runs)

    /* 4. Download original GenBank assemblies*/
    FETCH_ORIGINAL_ASSEMBLY(ch_assembly_metadata)

    /* 5. Baseline statistics for original assemblies*/
    STATS_ORIGINAL(FETCH_ORIGINAL_ASSEMBLY.out.assembly, 'original')
    CHECKM2_ORIGINAL(FETCH_ORIGINAL_ASSEMBLY.out.assembly, 'original')
    BUSCO_ORIGINAL(
        FETCH_ORIGINAL_ASSEMBLY.out.assembly,
        params.busco_lineage,
        'original'
    )

    /* 6. Read QC*/
    FASTP(MERGE_RUNS.out.reads)

    /* 7. Kraken2 read classification*/
    KRAKEN2(
        FASTP.out.reads,
        file(params.kraken2_db, checkIfExists: true)
    )

    ch_extract = FASTP.out.reads
        .join(KRAKEN2.out.classified)
        .join(KRAKEN2.out.report)

    /* 8. Extract Vibrio-classified read pairs*/
    EXTRACT_VIBRIO_READS(
        ch_extract,
        params.vibrio_taxid
    )

    /* 9. Reassemble*/
    UNICYCLER(EXTRACT_VIBRIO_READS.out.reads)

    /* 10. Assess the reassembly before contig cleaning*/
    STATS_PRE(UNICYCLER.out.scaffolds, 'pre')
    CHECKM2_PRE(UNICYCLER.out.scaffolds, 'pre')

    /* 11. Map reads back and filter by coverage and length*/
    ch_cov = UNICYCLER.out.scaffolds
        .join(EXTRACT_VIBRIO_READS.out.reads)

    COV_LENGTH_FILTER(ch_cov)

    /* 12. GUNC — post-assembly chimera + cross-clade contamination on the
     * coverage/length-filtered assembly. Its per-contig contamination_portion
     * and genome-level maxCSS verdict drive the cleaning policy.
     */
    GUNC(
        COV_LENGTH_FILTER.out.filtered,
        file(params.gunc_db, checkIfExists: true)
    )

    /* 13. Combine filtered assembly, coverage, and GUNC calls*/
    ch_policy = COV_LENGTH_FILTER.out.filtered
        .join(COV_LENGTH_FILTER.out.coverage)
        .join(GUNC.out.summary)
        .join(GUNC.out.contigs)

    /* 14. Apply the GUNC-driven cleaning policy*/
    APPLY_FCS_POLICY(ch_policy)

    /*15. Final quality assessment*/
    STATS_POST(APPLY_FCS_POLICY.out.clean_assembly, 'post')
    CHECKM2_POST(APPLY_FCS_POLICY.out.clean_assembly, 'post')

    BUSCO_POST(
        APPLY_FCS_POLICY.out.clean_assembly,
        params.busco_lineage,
        'post'
    )

    /*16. Pre -> post CheckM2 delta*/
    ch_delta = CHECKM2_PRE.out.quality_report
        .join(CHECKM2_POST.out.quality_report)

    CHECKM2_DELTA(ch_delta)

    /*17. Per-sample summary row (outer join — failed samples still appear)*/
    def na_file = file("${projectDir}/assets/NA.placeholder")
    ch_sample_summary = STATS_ORIGINAL.out.tsv
        .join(CHECKM2_ORIGINAL.out.quality_report, remainder: true)
        .join(BUSCO_ORIGINAL.out.summary, remainder: true)
        .join(STATS_POST.out.tsv, remainder: true)
        .join(CHECKM2_DELTA.out.tsv, remainder: true)
        .join(BUSCO_POST.out.summary, remainder: true)
        .join(APPLY_FCS_POLICY.out.decisions, remainder: true)
        .map { row ->
            def sid = row[0]
            [sid] + row[1..-1].collect { it == null ? na_file : it }
        }

    SAMPLE_SUMMARY(ch_sample_summary)

    /* 18. Report */
    COHORT_SUMMARY(
        SAMPLE_SUMMARY.out.tsv
            .map { _sid, f -> f }
            .collect()
    )

    ch_multiqc_input = Channel.empty()
        .mix(FASTP.out.json.map { sid, f -> f })
        .mix(KRAKEN2.out.report.map { sid, f -> f })
        .mix(UNICYCLER.out.log.map { sid, f -> f })
        .mix(BUSCO_ORIGINAL.out.summary.map { sid, f -> f })
        .mix(BUSCO_POST.out.summary.map { sid, f -> f })
        .mix(CHECKM2_ORIGINAL.out.quality_report.map { sid, f -> f })
        .mix(CHECKM2_PRE.out.quality_report.map { sid, f -> f })
        .mix(CHECKM2_POST.out.quality_report.map { sid, f -> f })
        .mix(GUNC.out.summary.map { sid, f -> f })
        .collect()

    MULTIQC(ch_multiqc_input)
}
