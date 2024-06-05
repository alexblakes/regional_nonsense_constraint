.ONESHELL:
.PHONY: all test

SHELL=/bin/bash

gnomad_coverage = data/interim/gnomad_v4.1_median_coverage.tsv.gz \
                  data/interim/gnomad_v4.1_median_coverage.tsv.gz.tbi \
                  data/manual/header_lines_coverage.txt \

vep_vcf = data/interim/cds_all_possible_snvs_vep.vcf.gz \
          data/interim/cds_all_possible_snvs_vep.vcf.gz.tbi \
		  data/manual/header_lines_vep_consequence.txt \

gnomad_tsv = data/manual/header_lines_gnomad.txt \
             data/interim/gnomad_v4.1_pass_snvs.tsv.gz \

snvs_annotated = data/interim/cds_all_possible_snvs_annotated.tsv.gz \
                 data/interim/cds_all_possible_snvs_annotated.tsv.gz.tbi \
				 data/logs/annotate_snvs.log \

test : $(gnomad_coverage)

all : $(gnomad_coverage) \
      $(vep_vcf) \
      $(gnomad_tsv) \
	  $(snvs_annotated) \

# Reformat gnomAD v4.1 exome coverage data for annotation
$(gnomad_coverage) &: src/snv_annotation/get_coverage.sh
	source activate vep
	bash src/snv_annotation/get_coverage.sh
	# Takes about 5 minutes to run

# Convert combined VEP output to VCF
$(vep_vcf) &: data/interim/cds_all_possible_snvs_vep.tsv.gz \
              src/snv_annotation/vep_tsv_to_vcf.sh
	source activate vep
	bash src/snv_annotation/vep_tsv_to_vcf.sh
	# Takes about 5 minutes to run

# Convert the gnomAD PASS SNVs to TSV
$(gnomad_tsv) &: data/interim/gnomad_v4.1_pass_snvs.vcf.gz \
                 src/snv_annotation/gnomad_pass_snvs_to_tsv.sh
	source activate vep
	bash src/snv_annotation/gnomad_pass_snvs_to_tsv.sh

# Merge annotations for all SNVs
$(snvs_annotated) &: data/interim/cds_all_possible_snvs_vep.vcf.gz \
                    data/interim/cds_all_possible_snvs_vep.tsv.gz \
                    data/final/nmd_annotations_simple.tsv.gz \
                    data/interim/gnomad_v4.1_pass_snvs.tsv.gz \
                    data/interim/gnomad_v4.1_median_coverage.tsv.gz \
                    data/interim/roulette/raw_cds_combined.tsv.gz \
                    data/interim/roulette/scaled_cds_combined_sorted.tsv.gz \
                    data/manual/header_lines_vep_consequence.txt \
                    data/manual/header_lines_nmd_annotation.txt \
                    data/manual/header_lines_gnomad.txt \
                    data/manual/header_lines_coverage.txt \
                    data/manual/header_lines_roulette_raw_tsv.txt \
                    data/manual/header_lines_roulette_scaled.txt \
					src/snv_annotation/annotate_snvs.sh
	source activate vep
	bash src/snv_annotation/annotate_snvs.sh
	# Takes about 10 minutes to run