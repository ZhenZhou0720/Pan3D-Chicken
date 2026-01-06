plink --vcf merge.vcf.snp.sv/merge.snp.sv.noZW.sorted.vcf.gz \
  --double-id \
  --r2 gz\
  --ld-snp-list sv_ids.txt \
  --ld-window-kb 50 \
  --ld-window 100000 \
  --ld-window-r2 0 \
  --chr-set 40 \
  --out sv.snp
