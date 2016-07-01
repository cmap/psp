#MAJOR VERSION 1.0
#MINOR VERSION 1.0
#8-JUL-2015

ProteomicConnectivityMapFromGCTFiles <- function (...,filter_level=0.7,filter_position='pre',connection_type='cnxn',zscore=FALSE,isochronic=FALSE) {
	##MAIN PROCESSING TO LAY FOUNDATION FOR NETWORK GRAPHS OR LAY FOUNDATION FOR KS-BASED QUERY
	##LIKELY EXTREMELY INEFFICIENTLY CODED
	##PROBABLY CAN BE SEPARATED INTO DIFFERENT FUNCTIONS
	##isochronic flag will only output similarities for perturbations at the same time point...not sure if we want to keep this here or move later in pipeline
	##I think several steps can be skipped if we only care about a query later.


	#going to aggregate multiple GCTs...either named as individual files or pulled from a directory
	#parse parameters
	inpars<-list(...);
	delisted<-unlist(inpars)
	gctind<-grep('gct',delisted)
	straightgcts<-delisted[gctind] #named GCT files
	evaluatefurther<-delisted[-gctind]
	if (length(gctind) == 0) {
		evaluatefurther<-delisted;
	}
	FI<-file.info(evaluatefurther);
	dirstocheck<-rownames(FI[which(FI$isdir==TRUE),])
	#checks any named directories for GCTs and adds to list
	if (length(dirstocheck>0)) {
		for (j in 1:length(dirstocheck)) {
		  morefiles<-paste(dirstocheck[j],list.files(dirstocheck[j])[grep('gct',list.files(dirstocheck[j]))],sep='/');
		  straightgcts<-c(straightgcts,morefiles)
		}
	}
	inpars<-straightgcts

	#get first file
	print("Reading first file...\n");
	MO<-P100provideGCTlistObjectFromFile(as.character(inpars[1]));

	#NOTE: look how simple z-scoring is here and why I'm not sure if it's a module in its own right
	if (zscore) {
		MO$dt<-t(apply(MO$dt,1,.zscore));
	}

	#merge GCT data tables and annotations
	if (length(inpars) > 1) {
		for (j in 2:length(inpars)) {
			t<-P100provideGCTlistObjectFromFile(as.character(inpars[j]))
			if (zscore) {
				t$dt<-t(apply(t$dt,1,.zscore));
			}
			MO$dt<-merge(MO$dt,t$dt,by=0,all=TRUE,sort=FALSE);
			MO$dt<-data.frame(MO$dt[-1],row.names=MO$dt$Row.names);
			MO$surviving_headers<-merge(MO$surviving_headers,t$surviving_headers,by=0,all=TRUE,sort=FALSE);
			MO$surviving_headers<-data.frame(MO$surviving_headers[-1],row.names=MO$surviving_headers$Row.names);

			#NOTE: surviving_rowAnnots is a hard case...we don't really want a merge, we just want extension of the data.frame if a new set is added
			#This is probably not critical for computing the connectivity map, but maybe can be better addressed in re-write.
			#MO$surviving_rowAnnots<-merge(MO$surviving_rowAnnots,t$surviving_rowAnnots,by=0,all=TRUE,sort=FALSE);
		}
	}

	if (zscore) {
		MO$surviving_headers<-.updateProvenanceCode(MO$static_headers,MO$surviving_headers,"ZSC");
	}

	#compute SIMILARITY metric - right now just correlation, but maybe replace with function later
	#NOTE - RYAN USES SPEARMAN
	C<-cor(MO$dt,use='pairwise.complete.obs')
	C[upper.tri(C,diag=TRUE)]<-NA

	#filter connections - maybe add option to do this prior to collapse or after collapse?
	#filtering is only useful for limiting stuff that goes into graph networks
	#for downstream queries the filter should always be set to -1 or whatever the lowest value of the similarity score metric might be
	fC<-C;
	if (filter_position=='pre') {
	  fC[fC<filter_level]<-NA;
	}
	#NOTE: I HAVE USED BAD NOMENCLATURE!  I REALLY MEAN "similarityScore" below.  Please forgive me!  And fix in your rewrite.
	fpC<-as.data.frame.table(fC,responseName='connectivityScore');
	fpC<-fpC[!(is.na(fpC$connectivityScore)),];

	#encode desired connectivity groups
	annotTable<-t(MO$surviving_headers);
	connGroupAnnots<-c('pert_iname','cell_id','pert_time'); #this defines the unique entity level for collapse of replicates...consider including dose?
	grp1<-paste0(connGroupAnnots,".1");
	grp2<-paste0(connGroupAnnots,".2");
	M1<-merge(fpC,annotTable[,connGroupAnnots],by.x='Var1',by.y=0,all.x=TRUE,sort=FALSE);
	M2<-merge(M1,annotTable[,connGroupAnnots],by.x='Var2',by.y=0,all.x=TRUE,sort=FALSE,suffixes=c('.1','.2'));
	M2$encoding.1<-paste(M2[,grp1[1]],M2[,grp1[2]],sep="_");
	M2$encoding.2<-paste(M2[,grp2[1]],M2[,grp2[2]],sep="_");
    M2$interaction<-apply(M2,1,regularizeNamingConvention) #THIS SEEMS TO BE TIME CONSUMING
    if (isochronic) {
    	M2<-M2[M2$pert_time.1==M2$pert_time.2,];
    }
	time.encoding<-paste(M2$pert_time.1,M2$pert_time.2,sep="-")
	M2<-data.frame(M2,time.encoding=time.encoding)

	#collapse based on connectivity groups - again, only useful for network graph visualizations at the similarity level
	cI<-data.frame(intxnName=unique(M2$interaction));
	cI$collapsedScore<-apply(cI,1,collapseInteractionScores,fulldata=M2);
	tM<-matrix(unlist(strsplit(as.character(cI$intxnName),split=":")),ncol=2,byrow=TRUE);
	cI$source<-tM[,1];
	cI$dest<-tM[,2];
	cI$type<-connection_type;

	if (filter_position=='post') {
		cI<-cI[which(cI$collapsedScore >= filter_level),];
	}

	#output interaction map - writes some files that are mostly just setup for Cytoscape
	dT<-cI[order(cI$collapsedScore,decreasing=TRUE),]
	stamp<-format(Sys.time(), "%Y-%m-%d_%H-%M-%S");
	fN<-paste0("network_",stamp,".txt");
	afN<-paste0("annot_",stamp,".txt");
	cfN<-paste0("funccall_",stamp,".txt");
	write.table(dT,fN,row.names=FALSE,sep="\t",quote=FALSE);
	fJ<-paste(c(...),collapse=",");
	cat('ProteomicConnectivityMapFromGCTFiles(',fJ,',filter_level=',filter_level,',filter_position=\'',filter_position,'\')',file=cfN,sep=''); #I think this somehow writes the command that was issued, but I'm not sure where it goes!
	#output annotation file
	#todo, this could be tricky to make sure all nodes are annotated once and only once.
	uN<-unique(c(cI$source,cI$dest));
	SuN<-matrix(unlist(strsplit(uN,split="_")),ncol=2,byrow=TRUE);
	AT<-data.frame(name=uN,drug=SuN[,1],cell=SuN[,2]);
	write.table(AT,afN,row.names=FALSE,sep="\t",quote=FALSE);


	#return a graph structure?  plot?  return a data frame?
	#SEVERAL THINGS GET RETURNED HERE:
	#mergedObject - useful for writing a merged GCT if you want to, but maybe there's something better in gctoo
	#corrMatrix - self-explanitory, simlarity matrix
	#filteredConnections - useful for network graph
	#tabulatedConnections=M2 - useful for KS-query
	#collapsedConnectivities - useful for network graph

	return(list(mergedObject=MO,corrMatrix=C,filteredConnections=fpC,tabulatedConnections=M2,collapsedConnectivities=cI));


}

regularizeNamingConvention <- function (tableObject) {
	eN<-c(tableObject['encoding.1'],tableObject['encoding.2']);
	sN<-sort(eN);
	rS<-paste(sN[1],sN[2],sep=":");
	return(rS);
}

collapseInteractionScores <- function (comboName,fulldata) {
	return(mean(fulldata$connectivityScore[which(fulldata$interaction==comboName)]));
}

PCCSEwriteMergedGCT <- function (output.name,mergedObject) {

	PCCSEwriteGCT (output.name, mergedObject$surviving_headers, rownames(mergedObject$dt), mergedObject$dt, 1, dim(mergedObject$surviving_headers)[1])
}

PCCSEwriteGCT <- function (output.name,surviving_headers,surviving_rowAnnots,normalizedData,colsAnnot=1,rowsAnnot) {

  all_data<-c(as.data.frame(surviving_rowAnnots),as.data.frame(normalizedData));
  line1<-'#1.3';
  line2<-c(dim(normalizedData),(colsAnnot-1),(rowsAnnot-1));
  line2<-t(line2);

  write.table(line1,output.name,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
  write.table(line2,output.name,sep="\t",row.names=FALSE,col.names=FALSE,append=TRUE,quote=FALSE)
  write.table(surviving_headers,output.name,sep="\t",row.names=TRUE,col.names=FALSE,append=TRUE,quote=FALSE)
  write.table(all_data,output.name,sep="\t",row.names=FALSE,col.names=FALSE,append=TRUE,quote=FALSE)
}

PCCSEmergeGCTFiles <- function (...,zscore=FALSE) {

	#parse parameters
	inpars<-list(...);
	omitnames<-c('pr_probe_normalization_group','pr_probe_suitability_manual');

	#get first file
	print("Reading first file...\n");
	MO<-P100provideGCTlistObjectFromFile(as.character(inpars[1]));
	if (zscore) {
		MO$dt<-t(apply(MO$dt,1,.zscore));
	}
	omitcols<-which(colnames(MO$surviving_rowAnnots) %in% omitnames);
	MO$surviving_rowAnnots<-MO$surviving_rowAnnots[,-omitcols];
	MO$static_headers<-MO$static_headers[,-omitcols];
	MO$colsAnnot<-MO$colsAnnot-length(omitcols);

	#merge GCT data tables and annotations
	if (length(inpars) > 1) {
		for (j in 2:length(inpars)) {
			t<-P100provideGCTlistObjectFromFile(as.character(inpars[j]))
			if (zscore) {
				t$dt<-t(apply(t$dt,1,.zscore));
			}
			MO$dt<-merge(MO$dt,t$dt,by=0,all=TRUE,sort=FALSE);
			MO$dt<-data.frame(MO$dt[-1],row.names=MO$dt$Row.names);
			MO$surviving_headers<-merge(MO$surviving_headers,t$surviving_headers,by=0,all=TRUE,sort=FALSE);
			MO$surviving_headers<-data.frame(MO$surviving_headers[-1],row.names=MO$surviving_headers$Row.names);
			
			omitcols<-which(colnames(t$surviving_rowAnnots) %in% omitnames);
			t$surviving_rowAnnots<-t$surviving_rowAnnots[,-omitcols];
			MO$surviving_rowAnnots<-merge(MO$surviving_rowAnnots,t$surviving_rowAnnots,all=TRUE,sort=FALSE);
			rownames(MO$surviving_rowAnnots)<-MO$surviving_rowAnnots$id;
			#NOTE: surviving_rowAnnots is a hard case...we don't really want a merge, we just want extension of the data.frame if a new set is added
			#This is probably not critical for computing the connectivity map, though, so lets ignore for now
			#MO$surviving_rowAnnots<-merge(MO$surviving_rowAnnots,t$surviving_rowAnnots,by=0,all=TRUE,sort=FALSE);
		}
	}

	stamp<-format(Sys.time(), "%Y-%m-%d_%H-%M-%S");
	fN<-paste0("mergedGCT_",stamp,".gct");

	#maybe need to doublecheck that all rowsAnnot ids match dt ids?
	MO$outputData<-MO$dt
    
	P100writeGCTForProcessedObject(fN,MO);
	return(MO);
}