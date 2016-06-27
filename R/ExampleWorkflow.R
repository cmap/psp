source('~/code/broadinstitute.psp/R/PCCSE_Connectivity.R')
source('~/code/broadinstitute.psp/R/p100_processing.R')
source('~/code/broadinstitute.psp/R/PCCSE_scratch_development.R')
# library(reshape2)

A<-ProteomicConnectivityMapFromGCTFiles('~/psp_data/all_gcp_QCNORM_n1330x59.gct',filter_level=-1,connection_type='GCP',zscore=TRUE)
D<-doAllRankAnalysis(A$tabulatedConnections,mode='datadump')
K<-ksTestForAll(D,assay='GCP')

d.PC3<-D[D$cell_id.1=='PC3',]

ct.PC3.KSstat<-dcast(k.PC3, queryDrug ~ otherDrug, value.var='KSstat')
ct.PC3.pval<-dcast(k.PC3, queryDrug ~ otherDrug, value.var='pval')
ct.PC3.directedKSstat<-dcast(k.PC3, queryDrug ~ otherDrug, value.var='directedKSstat')

write.csv(ct.PC3.KSstat,"PC3.connectivityMatrix.KSstat.csv")
write.csv(ct.PC3.directedKSstat,"PC3.connectivityMatrix.directedKSstat.csv")
write.csv(ct.PC3.pval,"PC3.connectivityMatrix.pval.csv")

q()
