#!/bin/bash
set -uex

checkPath() {
    if [ ! -d $1 ]; then mkdir -p $1; fi
}

PROJECT_DIR=/public/data_2/lxy/alto/project/brucella_wgs_Yunnan_2024

METADATA_DIR=${PROJECT_DIR}/metadata
FILE_LIST=${METADATA_DIR}/rawdata_filepath.tsv

SUB_PROJECT_DIR=${PROJECT_DIR}/pipelines/denovo_ssembly
DATA_DIR=${SUB_PROJECT_DIR}/data
SOURCE_DIR=${SUB_PROJECT_DIR}/source

FASTQC_OUT_DIR=${DATA_DIR}/fastqc
checkPath ${FASTQC_OUT_DIR}

MULTIQC_OUT_DIR=${DATA_DIR}/multiqc
checkPath ${MULTIQC_OUT_DIR}

FASTP_OUT_DIR=${DATA_DIR}/fastp
checkPath ${FASTP_OUT_DIR}

SPADES_OUT_DIR=${DATA_DIR}/spades_result
checkPath ${SPADES_OUT_DIR}

UNICYCLER_OUT_DIR=${DATA_DIR}/unicycler_result
checkPath ${UNICYCLER_OUT_DIR}

MLST_OUT_DIR=${DATA_DIR}/mlst_result
checkPath ${MLST_OUT_DIR}

PROKKA_OUT_DIR=${DATA_DIR}/prokka_result
checkPath ${PROKKA_OUT_DIR}

ROARY_OUT_DIR=${DATA_DIR}/roary_result
## roary 会自动创建文件夹
## checkPath ${ROARY_OUT_DIR}

TRIMAL_OUT_DIR=${DATA_DIR}/trimal_result
checkPath ${TRIMAL_OUT_DIR}

RAXML_OUT_DIR=${DATA_DIR}/raxml_result
checkPath ${RAXML_OUT_DIR}

TMP_DIR=${DATA_DIR}/seqkit
checkPath ${TMP_DIR}

REMOVED_GAP_DIR=${DATA_DIR}/core_genes_removed_gap
checkPath ${REMOVED_GAP_DIR}

UNICYCLER_ALL_FNA=${DATA_DIR}/unicycler_and_source_fna
checkPath ${UNICYCLER_ALL_FNA}

SNIPPY_CONTIG_OUT=${DATA_DIR}/snippy_contig_alignment
checkPath ${SNIPPY_CONTIG_OUT}

SNIPPY_CONTIG_GUBBINS_CLEAN=${DATA_DIR}/snippy_contig_gubbins
checkPath ${SNIPPY_CONTIG_GUBBINS_CLEAN}

RAXML_CONTIG_GUBBINS_DIR=${DATA_DIR}/raxml_gubbins_contig
checkPath ${RAXML_OUT_DIR}

BEAST_CONTIG_GUBBINS_DIR=${DATA_DIR}/beast_gubbins_contig
checkPath ${RAXML_OUT_DIR}

DIAMOND_DIR=${DATA_DIR}/diamond
checkPath ${DIAMOND_DIR}

COREGENE_MERGE_DIR=${DATA_DIR}/core_genes_fasta
checkPath ${COREGENE_MERGE_DIR}

SNIPPY_COREGENE_OUT_DIR=${DATA_DIR}/coregene_merge_snippy_alignment
checkPath ${SNIPPY_COREGENE_OUT_DIR}

PANAROO_OUT_DIR=${DATA_DIR}/panaroo_result
checkPath ${PANAROO_OUT_DIR}

# PANAROO_OUTGROUP_OUT_DIR=${DATA_DIR}/panaroo_result_add_outgroup
# checkPath ${PANAROO_OUTGROUP_OUT_DIR}

PANAROO_TREE_DIR=${DATA_DIR}/raxml_panaroo
checkPath ${PANAROO_TREE_DIR}

FASTANI_OUT=${DATA_DIR}/fastANI_result
checkPath ${FASTANI_OUT}

PANAROO_EMAPPER_DIAMOND_OUT_DIR=${DATA_DIR}/panaroo_emapper_diamond
checkPath ${PANAROO_EMAPPER_DIAMOND_OUT_DIR}

PANAROO_SENSITIVE_OUT_DIR=${DATA_DIR}/panaroo_result_sensitive
checkPath ${PANAROO_SENSITIVE_OUT_DIR}

PANAROO_SENSITIVE_75_OUT_DIR=${DATA_DIR}/panaroo_result_sensitive_75
checkPath ${PANAROO_SENSITIVE_75_OUT_DIR}

PANAROO_SENSITIVE_80_OUT_DIR=${DATA_DIR}/panaroo_result_sensitive_80
checkPath ${PANAROO_SENSITIVE_80_OUT_DIR}

PANAROO_SENSITIVE_85_OUT_DIR=${DATA_DIR}/panaroo_result_sensitive_85
checkPath ${PANAROO_SENSITIVE_85_OUT_DIR}

PANAROO_SENSITIVE_90_OUT_DIR=${DATA_DIR}/panaroo_result_sensitive_90
checkPath ${PANAROO_SENSITIVE_90_OUT_DIR}

PANAROO_SENSITIVE_95_OUT_DIR=${DATA_DIR}/panaroo_result_sensitive_95
checkPath ${PANAROO_SENSITIVE_95_OUT_DIR}

PANAROO_75_OUT_DIR=${DATA_DIR}/panaroo_result_75
checkPath ${PANAROO_75_OUT_DIR}

PANAROO_80_OUT_DIR=${DATA_DIR}/panaroo_result_80
checkPath ${PANAROO_80_OUT_DIR}

PANAROO_85_OUT_DIR=${DATA_DIR}/panaroo_result_85
checkPath ${PANAROO_85_OUT_DIR}

PANAROO_90_OUT_DIR=${DATA_DIR}/panaroo_result_90
checkPath ${PANAROO_90_OUT_DIR}

PANAROO_95_OUT_DIR=${DATA_DIR}/panaroo_result_95
checkPath ${PANAROO_95_OUT_DIR}

PANAROO_PAN_OUT_DIR=${DATA_DIR}/panaroo_result_strict_pan
checkPath ${PANAROO_PAN_OUT_DIR}

PANAROO_EMAPPER_DIAMOND_PAN_OUT_DIR=${DATA_DIR}/panaroo_result_strict_pan
checkPath ${PANAROO_EMAPPER_DIAMOND_PAN_OUT_DIR}
