vep \
    --assembly GRCh38 \
    --offline \
    --cache \
    --dir_cache /mnt/bmh01-rds/Ellingford_gene/.vep \
    --no_stats \
    --canonical \
    --id 'chr14 102001345 102001345 A G' \
    --vcf \
    -o STDOUT \
| filter_vep \
    --filter "CANONICAL is YES" \
    --filter "Consequence regex synonymous|missense|stop_gained" \
| less -S