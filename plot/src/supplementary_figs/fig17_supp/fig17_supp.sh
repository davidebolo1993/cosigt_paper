mkdir -p fig17_supp
Rscript \
  fig17_supp.r \
  fig17_supp/fig17_supp \
  ../../../../resources/benchmarking/leave_zero_out/svs_predictions/cosigt_locityper_svs_hgsvcv3_30X.tsv.gz \
  ../../../../resources/benchmarking/svs/svs_refined.annotation.bed
