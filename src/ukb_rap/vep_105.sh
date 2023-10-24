#! usr/bin/env bash

docker pull ensemblorg/ensembl-vep:release_105.0

dx download -r "$DX_PROJECT_CONTEXT_ID:/.vep"
tar xzf .vep/homo_sapiens_vep_105_GRCh38.tar.gz --directory=.vep
rm -rf .vep/homo_sapiens_vep_105_GRCh38.tar.gz

