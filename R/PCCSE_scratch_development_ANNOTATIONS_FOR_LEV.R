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
		workingData<-data.frame(dataset[cellAndDrugInd,],workingDrug=drug,otherDrug=drug,stringsAsFactors=FALSE); #this is the only thing we need if doAllRankAnalysis is in 'datadump' mode
		wI<-which(workingData$pert_iname.1==drug)
		workingData$otherDrug[-wI]<-as.character(unlist(workingData$pert_iname.1[-wI]));
		workingData$otherDrug[wI]<-as.character(unlist(workingData$pert_iname.2[wI]));
		selfCorrelInd<-workingData$interaction==paste(encoding,encoding,sep=":")
		percentileRanks<-1-rank(workingData$connectivityScore)[selfCorrelInd]/max(rank(workingData$connectivityScore))
		return(list(drug=drug,cell=cell,replicateData=workingData[selfCorrelInd,],replicatePercentileRanks=percentileRanks,meanPercentileRank=mean(percentileRanks),workingData=workingData,workingIndex=wI));
	} else {
		return(NULL) # if can't find any potential matches
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
			r<-doSingleRankAnalysis(dataset=dataset,drug=alldrugs[j],cell=allcells[k]); #if mode is 'datadump' we don't need most of what goes on in this functinon and it can be made vastly simpler
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
  #wp<-ggplot(raDataDump, aes(factor(workingDrug),connectivityScore));
  wp<-ggplot(raDataDump, aes(factor(otherDrug),connectivityScore));
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
  	ksr<-ks.test(qd$connectivityScore[which(qd$otherDrug==od[j])],raDataDump$connectivityScore[which(raDataDump$otherDrug==od[j])]);
  	
  	# Print some output
  	print(paste0("query: ", drug))
  	print(paste0("target: ", od[j]))
  	print(paste0("test_distrib: ", qd$connectivityScore[which(qd$otherDrug==od[j])]))
  	print(paste0("test_distrib_interactions: ", qd$interaction[which(qd$otherDrug==od[j])]))
  	print(paste0("num null_distrib elements: ", length(raDataDump$connectivityScore[which(raDataDump$otherDrug==od[j])])))
  	print(paste0("null_distib: ", raDataDump$connectivityScore[which(raDataDump$otherDrug==od[j])]))
  	print(paste0("null_distib_interactions: ", raDataDump$interaction[which(raDataDump$otherDrug==od[j])]))
  	print(paste0("KS-test stat: ", ksr$statistic))
  	print(paste0("KS-test pval: ", ksr$p.value))
  	
  	med.test<-median(qd$connectivityScore[which(qd$otherDrug==od[j])]);
  	med.bkg<-median(raDataDump$connectivityScore[which(raDataDump$otherDrug==od[j])])
  	rdf$KSstat[j]<-ksr$statistic;
  	rdf$pval[j]<-ksr$p.value;
  	rdf$med.test[j]<-med.test;
  	rdf$med.bkg[j]<-med.bkg;
  	rdf$direction[j]<-(1);
  	if (med.test < med.bkg) {
  		rdf$direction[j]<-(-1);
  	}
  	rdf$directedKSstat[j]<-rdf$direction[j]*rdf$KSstat[j]
  }
  rdf<-rdf[order(-rdf$directedKSstat),];
  return(rdf);
}

ksTestForAll <- function(raDataDump,assay="Assay") {
	drugs<-unique(raDataDump$workingDrug);
	ksTable<-queryDrugAndKStest(raDataDump,drug=drugs[1],assay=assay);
	for (j in 2:length(drugs)) {
		ks_temp<-queryDrugAndKStest(raDataDump,drug=drugs[j],assay=assay);
		ksTable<-rbind(ksTable,ks_temp)
	}
	kso<-order(ksTable$queryDrug,-ksTable$directedKSstat);
	return(ksTable[kso,]);
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

KSenabledTimeSeriesQueryDrugAndPlotPCCSEData <- function(raDataDump, drug, presortedKS, assay="Assay") {
  qd<-raDataDump[which(raDataDump$workingDrug==drug),];
  FO<-rev(presortedKS$otherDrug[which(presortedKS$queryDrug==drug)]);
  raDataDump<-raDataDump[order(factor(raDataDump$otherDrug,levels=FO)),]
  print(paste0('FACTOR ORDER SIZE ',length(FO)))
  wp<-ggplot(raDataDump, aes(factor(otherDrug,levels=factor(otherDrug)),connectivityScore));  ####NOT RECOGNIZING factorOrder for some reason?
  #wp<-ggplot(raDataDump, aes(factor(workingDrug),connectivityScore));
  wp + geom_violin() + 
    geom_boxplot(outlier.size=0,width=0.2) + 
    #geom_point( aes(factor(otherDrug,levels=factorOrder),connectivityScore,color=factor(cell_id.1)),data=qd,size=4, shape=124 ) + ####NOT RECOGNIZING factorOrder for some reason?
    geom_point( aes(factor(otherDrug),connectivityScore,color=factor(time.encoding)),data=qd,size=4, shape=124 ) + 
    theme(axis.text.x = element_text(angle=60, vjust=1)) + xlab("compound") +
    ggtitle(paste(assay," Query of ",drug,sep="")) + ylab("cScore Metric = Pearson of Z-score") + coord_flip()

}

KSenabledQueryDrugAndPlotPCCSEData <- function(raDataDump, drug, presortedKS, assay="Assay") {
  qd<-raDataDump[which(raDataDump$workingDrug==drug),];
  FO<-rev(presortedKS$otherDrug[which(presortedKS$queryDrug==drug)]);
  raDataDump<-raDataDump[order(factor(raDataDump$otherDrug,levels=FO)),]
  print(paste0('FACTOR ORDER SIZE ',length(FO)))
  wp<-ggplot(raDataDump, aes(factor(otherDrug,levels=factor(otherDrug)),connectivityScore));  ####NOT RECOGNIZING factorOrder for some reason?
  #wp<-ggplot(raDataDump, aes(factor(workingDrug),connectivityScore));
  wp + geom_violin() + 
    geom_boxplot(outlier.size=0,width=0.2) + 
    #geom_point( aes(factor(otherDrug,levels=factorOrder),connectivityScore,color=factor(cell_id.1)),data=qd,size=4, shape=124 ) + ####NOT RECOGNIZING factorOrder for some reason?
    geom_point( aes(factor(otherDrug),connectivityScore,color=factor(cell_id.1)),data=qd,size=4, shape=124 ) + 
    theme(axis.text.x = element_text(angle=60, vjust=1)) + xlab("compound") +
    ggtitle(paste(assay," Query of ",drug,sep="")) + ylab("cScore Metric = Pearson of Z-score") + coord_flip()

}



queryDrugResamplingKStest <- function(raDataDump, drug, assay="Assay",iterations=100) {
  ##RESAMPLES TO SMALLER BACKGROUND DISTRIBUTION
  qd<-raDataDump[which(raDataDump$workingDrug==drug),];
  od<-unique(qd$otherDrug);
  nod<-length(od)
  rdf<-data.frame(queryDrug=character(nod),otherDrug=character(nod),KSstat=numeric(nod),pval=numeric(nod),stringsAsFactors=FALSE)
  statList<-list();
  for (j in 1:length(od)) {
  	rdf$queryDrug[j]<-drug;
  	rdf$otherDrug[j]<-od[j];
  	rsKandP<-data.frame(statistic=numeric(iterations),p.value=numeric(iterations));
  	sampleSize<-length(qd$connectivityScore[which(qd$workingDrug==drug & qd$otherDrug==od[j])])
  	for (k in 1:iterations) {
	  	ksr<-ks.test(qd$connectivityScore[which(qd$otherDrug==od[j])],sample(raDataDump$connectivityScore[which(raDataDump$otherDrug==od[j])],sampleSize));
	  	rsKandP$statistic[k]<-ksr$statistic;
	  	rsKandP$p.value[k]<-ksr$p.value;
  	}
  	statList[[j]]<-rsKandP;
  	rdf$KSstat[j]<-mean(rsKandP$statistic);
  	rdf$pval[j]<-mean(rsKandP$p.value);
  }
  return(list(rdf=rdf,statList=statList));
}

queryDrugResamplingKStest2 <- function(raDataDump, drug, assay="Assay",iterations=100) {
  #RESAMPLES BACKGROUND DISTRIBUTION AT RANDOM TO EQUAL SIZE OF QUERY MATCH
  qd<-raDataDump[which(raDataDump$workingDrug==drug),];
  od<-unique(qd$otherDrug);
  nod<-length(od)
  canonical_rdf<-data.frame(queryDrug=character(nod),otherDrug=character(nod),KSstat=numeric(nod),pval=numeric(nod),ecdf_percentile=numeric(nod),sample_size=numeric(nod),stringsAsFactors=FALSE)
  rdf<-data.frame(queryDrug=character(nod),otherDrug=character(nod),KSstat=numeric(nod),pval=numeric(nod),stringsAsFactors=FALSE)
  statList<-list();
  for (j in 1:length(od)) {
  	canonical_rdf$queryDrug[j]<-drug;
  	canonical_rdf$otherDrug[j]<-od[j];
  	rdf$queryDrug[j]<-drug;
  	rdf$otherDrug[j]<-od[j];
  	rsKandP<-data.frame(statistic=numeric(iterations),p.value=numeric(iterations));
  	sampleSize<-length(qd$connectivityScore[which(qd$workingDrug==drug & qd$otherDrug==od[j])])

  	#cannonical KS test
  	ksr<-ks.test(qd$connectivityScore[which(qd$otherDrug==od[j])],raDataDump$connectivityScore[which(raDataDump$otherDrug==od[j])]);
  	med.test<-median(qd$connectivityScore[which(qd$otherDrug==od[j])]);
    med.bkg<-median(raDataDump$connectivityScore[which(raDataDump$otherDrug==od[j])])
    canonical_rdf$KSstat[j]<-ksr$statistic;
  	canonical_rdf$pval[j]<-ksr$p.value;
    canonical_rdf$med.test[j]<-med.test;
    canonical_rdf$med.bkg[j]<-med.bkg;
    canonical_rdf$direction[j]<-(1);
    if (med.test < med.bkg) {
      canonical_rdf$direction[j]<-(-1);
    }
    canonical_rdf$directedKSstat[j]<-canonical_rdf$direction[j]*canonical_rdf$KSstat[j]

    # now make a background distribution of KS scores sampling from the same size as the query
  	for (k in 1:iterations) {
	  	temp_ksr<-ks.test(sample(raDataDump$connectivityScore[which(raDataDump$otherDrug==od[j])],sampleSize),raDataDump$connectivityScore[which(raDataDump$otherDrug==od[j])]);
	  	rsKandP$statistic[k]<-temp_ksr$statistic;
	  	rsKandP$p.value[k]<-temp_ksr$p.value;
  	}

  	statList[[j]]<-rsKandP;
  	tempECDF<-ecdf(rsKandP$statistic);
  	canonical_rdf$ecdf_percentile[j]<-tempECDF(canonical_rdf$KSstat[j]);
  	canonical_rdf$sample_size[j]<-sampleSize;

  	rdf$KSstat[j]<-mean(rsKandP$statistic);
  	rdf$pval[j]<-mean(rsKandP$p.value);
  }
  return(list(canonical_rdf=canonical_rdf,rdf=rdf,statList=statList));
}

queryDrugResamplingKStest3 <- function(raDataDump, drug, assay="Assay",iterations=100) {
  #RESAMPLES BACKGROUND DISTRIBUTION AT RANDOM TO EQUAL SIZE OF QUERY MATCH
  qd<-raDataDump[which(raDataDump$workingDrug==drug),];
  od<-unique(qd$otherDrug);
  nod<-length(od)
  canonical_rdf<-data.frame(queryDrug=character(nod),otherDrug=character(nod),KSstat=numeric(nod),pval=numeric(nod),ecdf_percentile=numeric(nod),sample_size=numeric(nod),stringsAsFactors=FALSE)
  rdf<-data.frame(queryDrug=character(nod),otherDrug=character(nod),KSstat=numeric(nod),pval=numeric(nod),stringsAsFactors=FALSE)
  statList<-list();
  for (j in 1:length(od)) {
  	canonical_rdf$queryDrug[j]<-drug;
  	canonical_rdf$otherDrug[j]<-od[j];
  	rdf$queryDrug[j]<-drug;
  	rdf$otherDrug[j]<-od[j];
  	rsKandP<-data.frame(statistic=numeric(iterations),p.value=numeric(iterations));
  	sampleSize<-length(qd$connectivityScore[which(qd$workingDrug==drug & qd$otherDrug==od[j])])

  	#cannonical KS test
  	ksr<-ks.test(qd$connectivityScore[which(qd$otherDrug==od[j])],raDataDump$connectivityScore[which(raDataDump$otherDrug==od[j])]);
  	med.test<-median(qd$connectivityScore[which(qd$otherDrug==od[j])]);
    med.bkg<-median(raDataDump$connectivityScore[which(raDataDump$otherDrug==od[j])])
    canonical_rdf$KSstat[j]<-ksr$statistic;
    canonical_rdf$pval[j]<-ksr$p.value;
    canonical_rdf$med.test[j]<-med.test;
    canonical_rdf$med.bkg[j]<-med.bkg;
    canonical_rdf$direction[j]<-(1);
    if (med.test < med.bkg) {
      canonical_rdf$direction[j]<-(-1);
    }
    canonical_rdf$directedKSstat[j]<-canonical_rdf$direction[j]*canonical_rdf$KSstat[j]

  	bkgDistOtherDrug<-raDataDump$connectivityScore[which(raDataDump$otherDrug==od[j])];
  	ecdf_bkg<-ecdf(bkgDistOtherDrug);
  	queryDist<-qd$connectivityScore[which(qd$otherDrug==od[j])];
  	sdQueryDist<-sd(queryDist);

  	#set.seed(1234567);

    # now make a background distribution of KS scores sampling from the same size as the query
  	for (k in 1:iterations) {
  		currentCentroid<-sample(bkgDistOtherDrug,1);
  		fakePoints<-rnorm(sampleSize,mean=currentCentroid,sd=sdQueryDist)
  		temp_ksr<-ks.test(fakePoints,bkgDistOtherDrug);
      rsKandP$statistic[k]<-temp_ksr$statistic;
	  	rsKandP$p.value[k]<-temp_ksr$p.value;
      if (median(fakePoints) > med.bkg) {
        rsKandP$directedKSstat[k]<-temp_ksr$statistic;
      } else {
        rsKandP$directedKSstat[k]<-(-1)*temp_ksr$statistic;
      }
  	}

  	statList[[j]]<-rsKandP;
  	tempECDF<-ecdf(rsKandP$statistic);
  	canonical_rdf$ecdf_percentile[j]<-(1-tempECDF(canonical_rdf$KSstat[j]));
  	canonical_rdf$sample_size[j]<-sampleSize;

  	rdf$KSstat[j]<-mean(rsKandP$statistic);
  	rdf$pval[j]<-mean(rsKandP$p.value);
  }
  return(list(canonical_rdf=canonical_rdf,rdf=rdf,statList=statList));
}

RS3wrapper<-function (raDataDump, iterations=1000) {
	drugs<-unique(raDataDump$otherDrug);
	results_df<-queryDrugResamplingKStest3(raDataDump,drugs[1],iterations=iterations);
	canonical_rdf_object<-results_df$canonical_rdf;
  	for (k in 2: length(drugs)) {
		result<-queryDrugResamplingKStest3(raDataDump,drugs[k],iterations=iterations);
		canonical_rdf_object<-rbind(canonical_rdf_object,result$canonical_rdf);
	}
	return(canonical_rdf_object);
}


queryDrugResamplingKStest4 <- function(raDataDump, drug, assay="Assay",iterations=100) {
  #RESAMPLES BACKGROUND DISTRIBUTION AT RANDOM TO EQUAL SIZE OF QUERY MATCH
  qd<-raDataDump[which(raDataDump$workingDrug==drug),];
  od<-unique(qd$otherDrug);
  nod<-length(od)
  canonical_rdf<-data.frame(queryDrug=character(nod),otherDrug=character(nod),KSstat=numeric(nod),pval=numeric(nod),ecdf_percentile=numeric(nod),sample_size=numeric(nod),stringsAsFactors=FALSE)
  rdf<-data.frame(queryDrug=character(nod),otherDrug=character(nod),KSstat=numeric(nod),pval=numeric(nod),stringsAsFactors=FALSE)
  statList<-list();
  for (j in 1:length(od)) {
    canonical_rdf$queryDrug[j]<-drug;
    canonical_rdf$otherDrug[j]<-od[j];
    rdf$queryDrug[j]<-drug;
    rdf$otherDrug[j]<-od[j];
    rsKandP<-data.frame(statistic=numeric(iterations),p.value=numeric(iterations));
    sampleSize<-length(qd$connectivityScore[which(qd$workingDrug==drug & qd$otherDrug==od[j])])

    #cannonical KS test
    ksr<-ks.test(qd$connectivityScore[which(qd$otherDrug==od[j])],raDataDump$connectivityScore[which(raDataDump$otherDrug==od[j])]);
    med.test<-median(qd$connectivityScore[which(qd$otherDrug==od[j])]);
    med.bkg<-median(raDataDump$connectivityScore[which(raDataDump$otherDrug==od[j])])
    canonical_rdf$KSstat[j]<-ksr$statistic;
    canonical_rdf$pval[j]<-ksr$p.value;
    canonical_rdf$med.test[j]<-med.test;
    canonical_rdf$med.bkg[j]<-med.bkg;
    canonical_rdf$direction[j]<-(1);
    if (med.test < med.bkg) {
      canonical_rdf$direction[j]<-(-1);
    }
    canonical_rdf$directedKSstat[j]<-canonical_rdf$direction[j]*canonical_rdf$KSstat[j]

    bkgDistOtherDrug<-raDataDump$connectivityScore[which(raDataDump$otherDrug==od[j])];
    ecdf_bkg<-ecdf(bkgDistOtherDrug);
    queryDist<-qd$connectivityScore[which(qd$otherDrug==od[j])];
    sdQueryDist<-sd(queryDist);

    #set.seed(1234567);

    # sample label permutation
    n<-length(raDataDump$workingDrug);
    shuffledraDataDump<-raDataDump;
    for (k in 1:iterations) {
      shuffledraDataDump$connectivityScore<-sample(raDataDump$connectivityScore,n);
      sqd<-shuffledraDataDump[which(shuffledraDataDump$workingDrug==drug),];

      temp_ksr<-ks.test(sqd$connectivityScore[which(sqd$otherDrug==od[j])],shuffledraDataDump$connectivityScore[which(shuffledraDataDump$otherDrug==od[j])]);
      rsKandP$statistic[k]<-temp_ksr$statistic;
      rsKandP$p.value[k]<-temp_ksr$p.value;
      if (median(sqd$connectivityScore[which(sqd$otherDrug==od[j])]) > med.bkg) {
        rsKandP$directedKSstat[k]<-temp_ksr$statistic;
      } else {
        rsKandP$directedKSstat[k]<-(-1)*temp_ksr$statistic;
      }
    }

    statList[[j]]<-rsKandP;
    tempECDF<-ecdf(rsKandP$statistic);
    canonical_rdf$ecdf_percentile[j]<-(1-tempECDF(canonical_rdf$KSstat[j]));
    canonical_rdf$sample_size[j]<-sampleSize;

    rdf$KSstat[j]<-mean(rsKandP$statistic);
    rdf$pval[j]<-mean(rsKandP$p.value);
  }
  return(list(canonical_rdf=canonical_rdf,rdf=rdf,statList=statList));
}
