docker run \
-v /data/docker_local/psp/psp_on_clue.cfg:/usr/src/psp/proteomics-signature-pipeline/broadinstitute_psp/clue/psp_on_clue.cfg \
-v /data/docker_local/psp-data/qcnorm:/proteomics/qcnorm_dir \
-v /data/docker_local/psp-data/sim:/proteomics/sim_dir \
-v /data/docker_local/psp-data/examples/test_external_query_many_single_sample.gct:/proteomics/test_external_query_many_single_sample.gct \
-v /data/docker_local/psp-data/out:/out \
cmap/psp GCP True /proteomics/test_external_query_many_single_sample.gct /out /usr/src/psp/proteomics-signature-pipeline/broadinstitute_psp/clue/psp_on_clue.cfg pert_id cell_id

