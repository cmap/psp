P100_ProvideAnalysisObjectFromPanoramaWeb <- function (repAnnot, probeAnnot, dataTable) {
	dimsRA<-dim(repAnnot);
	dimsPA<-dim(probeAnnot);
	dimsDT<-dim(dataTable);
	local_surviving_headers<-repAnnot;
	local_dt<-dataTable;
	local_surviving_rowAnnots<-data.frame(id=rownames(probeAnnot),probeAnnot);
	local_static_headers<-as.data.frame(matrix(data='NA',ncol=dim(local_surviving_rowAnnots)[2],nrow=dim(local_surviving_headers)[1]));
	rownames(local_static_headers)<-rownames(local_surviving_headers);
	colnames(local_static_headers)<-colnames(local_surviving_rowAnnots);
	local_static_headers[,1]<-rownames(local_static_headers);
	return(list(surviving_headers=local_surviving_headers,static_headers=local_static_headers,surviving_rowAnnots=local_surviving_rowAnnots,dt=local_dt,colsAnnot=dim(local_static_headers)[2],rowsAnnot=dim(local_static_headers)[1]));
}