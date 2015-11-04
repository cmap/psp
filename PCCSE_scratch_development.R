mmSpecialPlotModelComponents <- function(mclustResult, inputData, histo=FALSE, nbreaks=40) {
	modelParams <- mmModelSummary(mclustResult);
	nComp <- mclustResult$G;
	upperRange = ceiling(as.numeric(quantile(inputData,0.99)))+0.3;
	lowerRange = floor(as.numeric(quantile(inputData,0.01)))-0.3;
	rangeSpan=c(lowerRange,upperRange);
	xcoords <- seq(lowerRange,upperRange,by=0.01);
	ydens = matrix(0,length(xcoords),nComp+1);
	colorCycle = c("red","green","blue","purple","orange","cyan","magenta","yellow");
	for (n in 1:nComp) {
		ydens[,n] <- dnorm(xcoords,mean=modelParams[2,n],sd=modelParams[3,n])*modelParams[1,n];
		ydens[,nComp+1] = ydens[,nComp+1]+ydens[,n];
	}
	ymax = max(ydens[,nComp+1]);
	for (k in 1:nComp) {
		plot(xcoords,ydens[,k],type="l",lty=2,lwd=2,col=colorCycle[k],xlim=rangeSpan,ylim=c(0,ymax),ylab="",xlab="");
		par(new=T);
	}
	plot(xcoords,ydens[,nComp+1],type="l",lwd=3,col="black",xlim=rangeSpan,ylim=c(0,ymax),main="Mixture Model Components",xlab="Pearson correlation coefficient",ylab="density");
	legendText <-c(colnames(modelParams),"Combined");
	legend("topright",legend=legendText,col=c(colorCycle[1:nComp],"black"),lwd=2);
	if (histo) {
		par(new=T);
		hist(inputData,breaks=nbreaks,xlim=rangeSpan,main="",xlab="",ylab="",axes=FALSE);
		axis(4);
	}
}


doClassComparisons <- function(d) {
	cmpClasses<-union(unique(d$pert_iname.1.class),unique(d$pert_iname.2.class))
	nCC<-length(cmpClasses);
	m<-matrix(data=0,nrow=nCC,ncol=nCC,dimnames=list(cmpClasses,cmpClasses))
	s<-m;
	for (j in 1:nCC) {
		for (k in 1:nCC) {
			c<-d$connectivityScore[(d$pert_iname.1.class==cmpClasses[j] & d$pert_iname.2.class==cmpClasses[k]) | (d$pert_iname.2.class==cmpClasses[j] & d$pert_iname.1.class==cmpClasses[k])];
			m[j,k]<-mean(c);
			s[j,k]<-sd(c);
		}
	}
	return(list(means=m,sds=s));
}


doSingleRankAnalysis <- function(dataset,drug,cell) {
	encoding<-paste(drug,cell,sep="_")
	cellDataInd<-dataset$cell_id.1==cell & dataset$cell_id.2==cell
	drugDataInd<-dataset$pert_iname.1==drug | dataset$pert_iname.2==drug
	cellAndDrugInd<-cellDataInd & drugDataInd
	if (sum(cellAndDrugInd)>0) {
		workingData<-data.frame(dataset[cellAndDrugInd,],workingDrug=drug,otherDrug=drug,stringsAsFactors=FALSE);
		wI<-which(workingData$pert_iname.1==drug)
		workingData$otherDrug[-wI]<-as.character(unlist(workingData$pert_iname.1[-wI]));
		workingData$otherDrug[wI]<-as.character(unlist(workingData$pert_iname.2[wI]));
		selfCorrelInd<-workingData$interaction==paste(encoding,encoding,sep=":")
		percentileRanks<-1-rank(workingData$connectivityScore)[selfCorrelInd]/max(rank(workingData$connectivityScore))
		return(list(drug=drug,cell=cell,replicateData=workingData[selfCorrelInd,],replicatePercentileRanks=percentileRanks,meanPercentileRank=mean(percentileRanks),workingData=workingData,workingIndex=wI));
	} else {
		return(NULL)
	}
}

doAllRankAnalysis <- function(dataset,mode='summary') {
	alldrugs<-unique(union(dataset$pert_iname.1,dataset$pert_iname.2));
	allcells<-unique(unique(dataset$cell_id.1,dataset$cell_id.2));
	total_combos<-length(alldrugs)*length(allcells);
	counter=1;
	R<-data.frame(drug=character(total_combos),cell=character(total_combos),percentile=numeric(total_combos),stringsAsFactors = FALSE);
	J<-data.frame();
	for (j in 1:length(alldrugs)) {
		for (k in 1:length(allcells)) {
			r<-doSingleRankAnalysis(dataset=dataset,drug=alldrugs[j],cell=allcells[k]);
			if (!(is.null(r))) {
				R$drug[counter]<-r$drug;
				R$cell[counter]<-as.character(r$cell);
				R$percentile[counter]<-r$meanPercentileRank;
				J<-rbind(J,r$workingData)
			} else {
				R$drug[counter]<-alldrugs[j];
				R$cell[counter]<-as.character(allcells[k]);
				R$percentile[counter]<-NA;
			}
			counter = counter + 1;
		}
	}
	if (mode=='datadump') {
		return(J)
	} else {
		return(R);
	}
}

queryDrugAndPlotPCCSEData <- function(raDataDump, drug, assay="Assay") {
  qd<-raDataDump[which(raDataDump$workingDrug==drug),];
  factorOrder<-sort(levels(factor(union(raDataDump$pert_iname.1,raDataDump$pert_iname.2))));
  print(paste0('FACTOR ORDER SIZE ',length(factorOrder)))
  #wp<-ggplot(raDataDump, aes(factor(workingDrug,levels=factorOrder),connectivityScore));  ####NOT RECOGNIZING factorOrder for some reason?
  wp<-ggplot(raDataDump, aes(factor(workingDrug),connectivityScore));
  wp + geom_violin() + 
    geom_boxplot(outlier.size=0,width=0.2) + 
    #geom_point( aes(factor(otherDrug,levels=factorOrder),connectivityScore,color=factor(cell_id.1)),data=qd,size=4, shape=124 ) + ####NOT RECOGNIZING factorOrder for some reason?
    geom_point( aes(factor(otherDrug),connectivityScore,color=factor(cell_id.1)),data=qd,size=4, shape=124 ) + 
    theme(axis.text.x = element_text(angle=60, vjust=1)) + xlab("compound") +
    ggtitle(paste(assay," Query of ",drug," Mapped on cScore Distribution",sep="")) + ylab("cScore Metric = Pearson of Z-score") + coord_flip()

}

queryDrugAndKStest <- function(raDataDump, drug, assay="Assay") {
  qd<-raDataDump[which(raDataDump$workingDrug==drug),];
  od<-unique(qd$otherDrug);
  nod<-length(od)
  rdf<-data.frame(queryDrug=character(nod),otherDrug=character(nod),KSstat=numeric(nod),pval=numeric(nod),stringsAsFactors=FALSE)
  for (j in 1:length(od)) {
  	rdf$queryDrug[j]<-drug;
  	rdf$otherDrug[j]<-od[j];
  	ksr<-ks.test(qd$connectivityScore[which(qd$workingDrug==drug & qd$otherDrug==od[j])],qd$connectivityScore[which(qd$workingDrug==drug)]);
  	rdf$KSstat[j]<-ksr$statistic;
  	rdf$pval[j]<-ksr$p.value;
  }
  return(rdf);
}

ksTestForAll <- function(raDataDump,assay="Assay") {
	drugs<-unique(raDataDump$workingDrug);
	ksTable<-queryDrugAndKStest(raDataDump,drug=drugs[1],assay=assay);
	for (j in 2:length(drugs)) {
		ks_temp<-queryDrugAndKStest(raDataDump,drug=drugs[1],assay=assay);
		ksTable<-rbind(ksTable,ks_temp)
	}
	return(ksTable);
}

replicateRecallPlotPCCSEData <- function(raDataDump, assay="Assay", metric="Pearson of Z-score",colorFactor='cell_id.1') {
	reps<-raDataDump[which(raDataDump$encoding.1==raDataDump$encoding.2),];
	ggp<-ggplot(raDataDump, aes(factor(workingDrug),connectivityScore));
	ggp +
	geom_violin() + 
	geom_boxplot(outlier.size=0,width=0.2) + 
	geom_point( aes(factor(workingDrug),connectivityScore,color=factor(cell_id.1)),data=reps,size=4,shape=124 ) + 
	theme(axis.text.x = element_text(angle=60, vjust=1)) + 
	ggtitle(paste(assay," Replicate Recall Mapped on cScore Distribution",sep="")) + 
	xlab("compound") + 
	ylab(paste0("cScore Metric = ",metric)) + 
	coord_flip()
}

timeSeriesReplicateRecallPlotPCCSEData <- function(raDataDump, assay="Assay", metric="Pearson of Z-score",colorFactor='cell_id.1') {
	reps<-raDataDump[which(raDataDump$encoding.1==raDataDump$encoding.2),];
	ggp<-ggplot(raDataDump, aes(factor(workingDrug),connectivityScore));
	ggp +
	geom_violin() + 
	geom_boxplot(outlier.size=0,width=0.2) + 
	geom_point( aes(factor(workingDrug),connectivityScore,color=factor(time.encoding)),data=reps,size=4,shape=124 ) + 
	theme(axis.text.x = element_text(angle=60, vjust=1)) + 
	ggtitle(paste(assay," Replicate Recall Mapped on cScore Distribution",sep="")) + 
	xlab("compound") + 
	ylab(paste0("cScore Metric = ",metric)) + 
	coord_flip()
}

timeSeriesQueryDrugAndPlotPCCSEData <- function(raDataDump, drug, assay="Assay") {
  qd<-raDataDump[which(raDataDump$workingDrug==drug),];
  factorOrder<-sort(levels(factor(union(raDataDump$pert_iname.1,raDataDump$pert_iname.2))));
  print(paste0('FACTOR ORDER SIZE ',length(factorOrder)))
  #wp<-ggplot(raDataDump, aes(factor(workingDrug,levels=factorOrder),connectivityScore));  ####NOT RECOGNIZING factorOrder for some reason?
  wp<-ggplot(raDataDump, aes(factor(workingDrug),connectivityScore));
  wp + geom_violin() + 
    geom_boxplot(outlier.size=0,width=0.2) + 
    #geom_point( aes(factor(otherDrug,levels=factorOrder),connectivityScore,color=factor(cell_id.1)),data=qd,size=4, shape=124 ) + ####NOT RECOGNIZING factorOrder for some reason?
    geom_point( aes(factor(otherDrug),connectivityScore,color=factor(time.encoding)),data=qd,size=4, shape=124 ) + 
    theme(axis.text.x = element_text(angle=60, vjust=1)) + xlab("compound") +
    ggtitle(paste(assay," Query of ",drug," Mapped on cScore Distribution",sep="")) + ylab("cScore Metric = Pearson of Z-score") + coord_flip()

}