#MAJOR VERSION 1.2
#MINOR VERSION 1
#17-FEB-2015

P100processGCTPanorama <- function (repAnnot, probeAnnot, dataTable, fileStub, optim=TRUE,fileOutput=TRUE,log2=TRUE,samplePctCutoff=0.8, probePctCutoff=0.9, probeSDCutoff=3, distSDcutoff=3) {
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
  a<-list(surviving_headers=local_surviving_headers,static_headers=local_static_headers,surviving_rowAnnots=local_surviving_rowAnnots,dt=local_dt,colsAnnot=dim(local_static_headers)[2],rowsAnnot=dim(local_static_headers)[1],gctFileName=fileStub);
  P100processGCT(a,optim=optim,fileOutput=fileOutput,log2=log2,samplePctCutoff=samplePctCutoff, probePctCutoff=probePctCutoff, probeSDCutoff=probeSDCutoff, distSDcutoff=distSDcutoff);
}

P100processGCTFromFile <- function (gctFileName,optim=TRUE,fileOutput=TRUE,log2=TRUE,samplePctCutoff=0.8, probePctCutoff=0.9, probeSDCutoff=3, distSDcutoff=3) {
  g<-P100_ReadGCT(gctFileName);
  local_static_headers<-g$headers[,1:g$colsAnnot];
  local_surviving_headers<-g$headers[,(g$colsAnnot+1):(g$colsAnnot+g$colsData)];
  local_surviving_rowAnnots<-g$rowAnnot;
  local_dt<-g$data;
  a<-list(surviving_headers=local_surviving_headers,static_headers=local_static_headers,surviving_rowAnnots=local_surviving_rowAnnots,dt=local_dt,colsAnnot=g$colsAnnot,rowsAnnot=g$rowsAnnot,gctFileName=gctFileName);
  P100processGCT(a,optim=optim,fileOutput=fileOutput,log2=log2,samplePctCutoff=samplePctCutoff, probePctCutoff=probePctCutoff, probeSDCutoff=probeSDCutoff, distSDcutoff=distSDcutoff);
}

P100processGCT <- function (g,optim=TRUE,fileOutput=TRUE,log2=TRUE,samplePctCutoff=0.8, probePctCutoff=0.9, probeSDCutoff=3, distSDcutoff=3) {

  static_headers<-g$static_headers;
  surviving_headers<-g$surviving_headers;
  surviving_rowAnnots<-g$surviving_rowAnnots;
  dt<-g$dt;
  gctFileName=g$gctFileName;

  #log2 transform
  if (log2) {
    dt[dt==0]=NA;
    dt<-log(dt)/log(2);
    surviving_headers<-.updateProvenanceCode(static_headers,surviving_headers,"L2X");
  }

  s<-P100filterSamplesForPoorCoverage(dt, pctFilterCutoff=samplePctCutoff)
  surviving_headers<-surviving_headers[,s$colsPassing];
  surviving_headers<-.updateProvenanceCode(static_headers,surviving_headers,"SF8");

  f<-P100filterProbes(s$filteredData, pctFilterCutoff=probePctCutoff, sdFilterCutoff=probeSDCutoff);
  surviving_rowAnnots<-surviving_rowAnnots[f$rowsPassing,];
  surviving_headers<-.updateProvenanceCode(static_headers,surviving_headers,"PF9");

  o<-P100optimizeSampleBalance(f$filteredData);

  b<-P100filterOutlierSamplesAndApplyCorrections(dt=f$filteredData,optData=o,sdFilterCutoff=distSDcutoff,optim=optim);
  surviving_headers<-surviving_headers[,b$colsPassing];
  if (optim) {
    surviving_headers<-.updateProvenanceCode(static_headers,surviving_headers,"LLB");
  }
  surviving_headers<-.updateProvenanceCode(static_headers,surviving_headers,"OSF");

  n<-P100rowMedianNormalize(b$filteredData);
  surviving_headers<-.updateProvenanceCode(static_headers,surviving_headers,"RMN");

  if (fileOutput) {
    ofn<-paste(gctFileName,'.processed.gct',sep='');
    P100writeGCT(output.name=ofn,static_headers=static_headers,surviving_headers=surviving_headers,surviving_rowAnnots=surviving_rowAnnots,normalizedData=n$normalizedData,colsAnnot=g$colsAnnot,rowsAnnot=g$rowsAnnot);
  }

  q<-list(inputData=dt,initialSampleFiltering=s,probeFiltering=f,optimizationParams=o,secondarySampleCorrectionAndFiltering=b,normalizedData=n,outputData=n$normalizedData,static_headers=static_headers,surviving_headers=surviving_headers,surviving_rowAnnots=surviving_rowAnnots);
  return(q);
}

quickP100processGCT <- function (gctFileName,optim=TRUE,fileOutput=TRUE,log2=TRUE) {

  g<-P100_ReadGCT(gctFileName);
  static_headers<-g$headers[,1:g$colsAnnot];
  surviving_headers<-g$headers[,(g$colsAnnot+1):(g$colsAnnot+g$colsData)];
  surviving_rowAnnots<-g$rowAnnot;
  dt<-g$data;

  #log2 transform
  if (log2) {
    dt[dt==0]=NA;
    dt<-log(dt)/log(2);
  }
  n<-P100rowMedianNormalize(dt);
  surviving_headers<-.updateProvenanceCode(static_headers,surviving_headers,"RMN");

  if (fileOutput) {
    ofn<-paste(gctFileName,'.quick.processed.gct',sep='');
    P100writeGCT(output.name=ofn,static_headers=static_headers,surviving_headers=surviving_headers,surviving_rowAnnots=surviving_rowAnnots,normalizedData=n$normalizedData,colsAnnot=g$colsAnnot,rowsAnnot=g$rowsAnnot);
  }

}

P100evaluateGCTasPlate <- function (gctFileName,optim=TRUE,fileOutput=TRUE,log2=TRUE,samplePctCutoff=0.8, probePctCutoff=0.9, probeSDCutoff=3, distSDcutoff=3) {

  require(ggplot2);
  require(grid);

  g<-P100_ReadGCT(gctFileName);
  static_headers<-g$headers[,1:g$colsAnnot];
  surviving_headers<-g$headers[,(g$colsAnnot+1):(g$colsAnnot+g$colsData)];
  surviving_rowAnnots<-g$rowAnnot;
  dt<-g$data;

  qq1<-plotSurvivingRowsLogical(surviving_headers,stage='Initial Data Read')
  #readline(prompt='next');

  #log2 transform
  if (log2) {
    dt[dt==0]=NA;
    dt<-log(dt)/log(2);
    surviving_headers<-.updateProvenanceCode(static_headers,surviving_headers,"L2X");
  }

  s<-P100filterSamplesForPoorCoverage(dt, pctFilterCutoff=samplePctCutoff)
  surviving_headers<-surviving_headers[,s$colsPassing];
  surviving_headers<-.updateProvenanceCode(static_headers,surviving_headers,"SF8");

  qq2<-plotSurvivingRowsLogical(surviving_headers,stage='Filter Poor Coverage')
  #readline(prompt='next');


  f<-P100filterProbes(s$filteredData, pctFilterCutoff=probePctCutoff, sdFilterCutoff=probeSDCutoff);
  surviving_rowAnnots<-surviving_rowAnnots[f$rowsPassing,];
  surviving_headers<-.updateProvenanceCode(static_headers,surviving_headers,"PF9");

  o<-P100optimizeSampleBalance(f$filteredData);
  qq3<-plotSurvivingRowsNumeric(surviving_headers,o$constants,stage='Correction Factors')
  b<-P100filterOutlierSamplesAndApplyCorrections(dt=f$filteredData,optData=o,sdFilterCutoff=distSDcutoff,optim=optim);
  surviving_headers<-surviving_headers[,b$colsPassing];
  if (optim) {
    surviving_headers<-.updateProvenanceCode(static_headers,surviving_headers,"LLB");
  }
  surviving_headers<-.updateProvenanceCode(static_headers,surviving_headers,"OSF");

  qq4<-plotSurvivingRowsLogical(surviving_headers,stage='Outliers Removed')
  #readline(prompt='next');

  multiplot(qq1,qq3,qq2,qq4,cols=2);


  n<-P100rowMedianNormalize(b$filteredData);
  surviving_headers<-.updateProvenanceCode(static_headers,surviving_headers,"RMN");

  if (fileOutput) {
    ofn<-paste(gctFileName,'.processed.gct',sep='');
    P100writeGCT(output.name=ofn,static_headers=static_headers,surviving_headers=surviving_headers,surviving_rowAnnots=surviving_rowAnnots,normalizedData=n$normalizedData,colsAnnot=g$colsAnnot,rowsAnnot=g$rowsAnnot);
  }

  q<-list(inputData=dt,initialSampleFiltering=s,probeFiltering=f,optimizationParams=o,secondarySampleCorrectionAndFiltering=b,normalizedData=n,outputData=n$normalizedData,static_headers=static_headers,surviving_headers=surviving_headers,surviving_rowAnnots=surviving_rowAnnots,corrFactors=o);
  return(q);
}

P100readAndDivideGCT <- function (gctFileName,log2=TRUE) {
  g<-P100_ReadGCT(gctFileName);
  static_headers<-g$headers[,1:g$colsAnnot];
  surviving_headers<-g$headers[,(g$colsAnnot+1):(g$colsAnnot+g$colsData)];
  surviving_rowAnnots<-g$rowAnnot;
  dt<-g$data;

  #log2 transform
  if (log2) {
    dt[dt==0]=NA;
    dt<-log(dt)/log(2);
  }

  q<-list(wholeGCT=g,dataTable=dt,static_headers=static_headers,surviving_headers=surviving_headers,surviving_rowAnnots=surviving_rowAnnots);
  return(q);
}

P100prepareForImputationModel <- function (p100object) {
  p100Data<-p100object$outputData;
  p100SampleNames<-p100object$surviving_headers[2,];
  p100Probes<-p100object$surviving_rowAnnots[,8];
  colnames(p100Data)<-p100SampleNames;
  rownames(p100Data)<-p100Probes;

  set.seed(8736697);

  p100DataImp<-p100Data;

  for (i in 1:dim(p100DataImp)[1]) {
    tmp_mean<-mean((p100DataImp[i, ]), na.rm = TRUE);
    tmp_std<-sd((p100DataImp[i, ]), na.rm = TRUE);
    p100DataImp[i, is.na(p100DataImp[i, ])] <- rnorm(length(p100DataImp[i, is.na(p100DataImp[i, ])]),mean=tmp_mean,sd=tmp_std)
  }

  q<-list(simpleFormat=p100Data,simpleFormatNAsImputed=p100DataImp);
  return(q);
}

P100dataNormalization <- function (dt,optim=TRUE) {
  #0. Filter out bad samples based on high number of missing values
  #1. Filter out bad rows based on #of NAs or high stdevs?
  #2. Optimize an additive constant for each column - only using non-NA cells (load balancing).
  #2.5 Filter out more bad samples based on far distance from optimum (might need to look at distro)
  #3. Reprocess filtered matrix using these constants.
  #4. Return an object with all of this stuff.

  s<-P100filterSamplesForPoorCoverage(dt)
  f<-P100filterProbes(s$filteredData);
  o<-P100optimizeSampleBalance(f$filteredData);
  b<-P100filterOutlierSamplesAndApplyCorrections(dt=f$filteredData,optData=o,optim=optim);
  n<-P100rowMedianNormalize(b$filteredData);

  q<-list(inputData=dt,initialSampleFiltering=s,probeFiltering=f,optimizationParams=o,secondarySampleCorrectionAndFiltering=b,normalizedData=n,outputData=n$normalizedData);

  return(q);
}

P100rowMedianNormalize <- function (dt) {
  retList<-list(rowMedians=NULL,normalizedData=dt,originalData=dt);
  dt<-as.matrix(dt);
  ddt<-dim(dt);
  nRows<-ddt[1];
  nCols<-ddt[2];
  retList$rowMedians<-.getRowMedians(dt);
  for (n in 1:nCols) {
    dt[,n]<-dt[,n]-retList$rowMedians;
  }
  retList$normalizedData<-dt;
  return(retList);
}

P100filterSamplesForPoorCoverage <- function (dt, pctFilterCutoff=0.8) {
  #Filter out columns missing x% of data;
  retList<-list(colsPassing=NULL,filteredData=dt,originalData=dt);
  dt<-as.matrix(dt);
  ddt<-dim(dt);
  nRows<-ddt[1];
  nCols<-ddt[2];

  numNAs<-.countColNAs(dt);
  pctNotNAs<-1-(numNAs/nCols);
  retList$colsPassing<-(pctNotNAs>=pctFilterCutoff);

  retList$filteredData<-dt[,retList$colsPassing];
  return(retList);

}

P100filterOutlierSamplesAndApplyCorrections <- function (dt, optData, sdFilterCutoff=3, optim=TRUE) {
  #Filter out columns missing x% of data;
  retList<-list(colsPassing=NULL,filteredData=dt,originalData=dt);
  dt<-as.matrix(dt);
  ddt<-dim(dt);
  nRows<-ddt[1];
  nCols<-ddt[2];

  if (optim) {
    print("Optimization of sample load balance performed.");
    for (n in 1:nCols) {
      dt[,n]<-dt[,n]+optData$constants[n];
    }
    cutoff<-mean(optData$optDistances)+sd(optData$optDistances)*sdFilterCutoff;
    retList$colsPassing<-optData$optDistances <= cutoff;
    retList$filteredData<-dt[,retList$colsPassing];
  } else {
    print("No sample load balancing was performed with these settings.");
    cutoff<-mean(optData$origDistances)+sd(optData$origDistances)*sdFilterCutoff;
    retList$colsPassing<-optData$origDistances <= cutoff;
    retList$filteredData<-dt[,retList$colsPassing];
  }

  return(retList);

}


P100filterProbes <-function (dt, pctFilterCutoff=0.9, sdFilterCutoff=3) {
  retList<-list(rowsPassing=NULL,filteredData=dt,originalData=dt);
  dt<-as.matrix(dt);
  ddt<-dim(dt);
  nRows<-ddt[1];
  nCols<-ddt[2];

  #Filter out rows missing x% of data;
  numNAs<-.countRowNAs(dt);
  pctNotNAs<-1-(numNAs/nRows);
  retList$rowsPassing<-(pctNotNAs>=pctFilterCutoff);

  #Filter out rows with high SD - remember each unit is 2 fold, multiplicative! 2,4,8
  rowSDs<-.getRowSDs(dt);
  retList$rowsPassing<-(retList$rowsPassing & rowSDs <= sdFilterCutoff);

  retList$filteredData<-dt[retList$rowsPassing,];
  return(retList);
}

P100optimizeSampleBalance <- function (dt) {
  dt<-as.matrix(dt);
  ddt<-dim(dt);
  nCols<-ddt[2];
  optim<-list(constants=array(dim=nCols),optDistances=array(dim=nCols),origDistances=array(dim=nCols))
  medianProfile<-.getRowMedians(dt);
  for (n in 1:nCols) {
    #determine non NA rows
    .good<-!(is.na(dt[,n]));
    originalDist<-.objectiveSSD(x=0,v1=dt[.good,n],v2=medianProfile[.good]);
    temp<-optimize(.objectiveSSD,c(-7,7),tol=0.0001,v1=dt[.good,n],v2=medianProfile[.good]);
    optim$constants[n]<-temp$minimum;
    optim$optDistances[n]<-temp$objective;
    optim$origDistances[n]<-originalDist;
  }
  return(optim);
}

P100mineAnnotation <- function (idlist,term="chromatin") {
	require(UniProt.ws)
	require(GO.db)
	getFromGO<-c('TERM','ONTOLOGY')
	uniprotResults<-select(UniProt.ws,idlist,'GO-ID','UNIPROTKB')
	uniAccessions<-uniprotResults$'UNIPROTKB'
	foundTerm<-vector(mode='logical',length=length(uniAccessions))
	uniGOcats<-strsplit(uniprotResults$'GO-ID',"; ")
	for (i in 1:length(uniGOcats)) {
		if (is.na(uniGOcats[[i]])) {
			foundTerm[i] = FALSE
		} else {
			GOresult<-select(GO.db,uniGOcats[[i]],getFromGO,"GOID")
			grepResult<-grep(term,GOresult$TERM)
			if (length(grepResult) > 0 ) {
				foundTerm[i] = TRUE
			} else {
				foundTerm[i] = FALSE
			}
		}
	}
	return(data.frame(UNIPROTKB=uniAccessions,FOUNDTERM=foundTerm))
}

.getRowMedians <- function(x) {
  x<-as.matrix(x);
  dx<-dim(x);
  nRows<-dx[1];
  y<-array(dim=nRows);
  for (n in 1:nRows) {
    y[n] <- median(x[n,],na.rm=TRUE)
  }
  return(y);
}

.getRowSDs <- function(x) {
  x<-as.matrix(x);
  dx<-dim(x);
  nRows<-dx[1];
  y<-array(dim=nRows);
  for (n in 1:nRows) {
    y[n] <- (sd(x[n,],na.rm=TRUE));
  }
  return(y);
}

.countRowNAs <- function(x) {
  x<-as.matrix(x);
  dx<-dim(x);
  nRows<-dx[1];
  y<-array(dim=nRows);
  for (n in 1:nRows) {
    y[n]<-sum(is.na(x[n,]));
  }
  return(y);
}

.countColNAs <- function(x) {
  x<-as.matrix(x);
  dx<-dim(x);
  nCols<-dx[2];
  y<-array(dim=nCols);
  for (n in 1:nCols) {
    y[n]<-sum(is.na(x[,n]));
  }
  return(y);
}


.objectiveSSD <- function(x, v1, v2) {
  return(sum(((v1+x)-v2)**2));
}

P100_ReadGCT <- function(gctFile) {
  firstPass<-read.delim(gctFile,header=FALSE,as.is=TRUE);
  rowsData<-as.numeric(firstPass[2,1]);
  colsData<-as.numeric(firstPass[2,2]);
  colsAnnot<-as.numeric(firstPass[2,3])+1;
  rowsAnnot<-as.numeric(firstPass[2,4])+1;
  colCls<-character(colsData+colsAnnot);
  colCls[1:(colsData+colsAnnot)]<-'character';
  headerRows<-read.delim(gctFile,header=FALSE,skip=2,colClasses=colCls,as.is=TRUE,nrows=rowsAnnot);
  #finish this later
  dataColumnClasses<-character(colsData+colsAnnot);
  dataColumnClasses[1:colsAnnot]<-'character';
  dataColumnClasses[(colsAnnot+1):(colsData+colsAnnot)]<-'numeric';
  dataRows<-read.delim(gctFile,header=FALSE,skip=2+rowsAnnot,colClasses=dataColumnClasses,as.is=TRUE,na.strings=c('NaN','#N/A!',"#N/A"));
  #dataSet<-as.numeric(unlist(dataRows[,(colsAnnot+1):(colsAnnot+colsData)]));
  #dim(dataSet)<-c(rowsData,colsData);
  dataSet<-dataRows[,(colsAnnot+1):(colsAnnot+colsData)];
  rowAnnotation<-dataRows[,1:colsAnnot];
  rownames(headerRows)<-headerRows[,1];
  parsed<-list(headers=headerRows,rowAnnot=rowAnnotation,data=dataSet,rowsData=rowsData,colsData=colsData,colsAnnot=colsAnnot,rowsAnnot=rowsAnnot);
}


P100writeGCT <- function (output.name,static_headers,surviving_headers,surviving_rowAnnots,normalizedData,colsAnnot,rowsAnnot) {

  all_headers<-c(static_headers,surviving_headers);
  all_data<-c(as.data.frame(surviving_rowAnnots),as.data.frame(normalizedData));
  line1<-'#1.3';
  line2<-c(dim(normalizedData),(colsAnnot-1),(rowsAnnot-1));
  line2<-t(line2);

  write.table(line1,output.name,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
  write.table(line2,output.name,sep="\t",row.names=FALSE,col.names=FALSE,append=TRUE,quote=FALSE)
  write.table(all_headers,output.name,sep="\t",row.names=FALSE,col.names=FALSE,append=TRUE,quote=FALSE)
  write.table(all_data,output.name,sep="\t",row.names=FALSE,col.names=FALSE,append=TRUE,quote=FALSE)
}

P100overlapNames<-function(somenames) {
  #for doing correlation matrices and overlapping by name
  a<-length(somenames);
  b<-matrix(data=FALSE,nrow=a,ncol=a);
  for (x in 1:a) {
    for (y in 1:a) {
      if (somenames[x] == somenames[y]) {
        b[x,y] = TRUE;
      }
    }
  }
  rownames(b)<-somenames;
  colnames(b)<-somenames;
  return(b);
}

genericRowNormalizeByColumnAnnot<-function(dt,columnAnnot,groupingField,subset=NULL) {
  normGroups<- unique(unlist(columnAnnot[groupingField,]));
  for (x in 1:length(normGroups)) {
    dtTemp<-dt[,columnAnnot[groupingField,]==normGroups[x]];
    groupNorm<-apply(dtTemp,1,median,na.rm=TRUE);
    dtTemp<-dtTemp-groupNorm;
    dt[,columnAnnot[groupingField,]==normGroups[x]]<-dtTemp;
  }
  return(dt);
}

.updateProvenanceCode <- function (static_annot,annot,code) {
	indarr<-static_annot[,1]=="provenance_code";
	annot[indarr,]<-paste(annot[indarr,],code,sep="+");
	return(annot);
}

theme_bdc_microtiter <- function (base_size=14, base_family="")
{
    theme(complete=TRUE,
          line = element_line(colour="grey70", size=0.5, linetype=1,
                              lineend="square"),
          rect = element_rect(fill="white", colour="grey70", size=0.5,
                              linetype=1),
          text = element_text(family=base_family, face="plain", colour="black",
                              size=base_size, hjust=0.5, vjust=0.5, angle=0,
                              lineheight=0.9),
          title = element_text(family=base_family, face="bold", colour="black",
                               vjust=0.0, hjust=0.5, angle=0),

          plot.background = element_rect(fill="transparent", colour=NA),
          plot.title = element_text(size=rel(1.2), vjust=0.8),
          plot.margin = unit(c(1, 1, 1, 1), "lines"),

          panel.background = element_rect(color='grey20', fill='white'),
          panel.border = element_rect(fill="transparent", color='grey20'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.margin = unit(0, "lines"),

          strip.background = element_rect(fill="grey70"),
          strip.text = element_text(size=rel(0.8)),
          strip.text.x = element_text(),
          strip.text.y = element_text(angle=-90),

          axis.text = element_text(face="bold"),
          axis.line = element_blank(),
          axis.text.x = element_text(),
          axis.text.y = element_text(hjust=1),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.ticks.length = unit(0, "cm"),
          axis.ticks.margin = unit(0.2, "cm"),

          legend.background = element_rect(fill="transparent", colour=NA),
          legend.margin = unit(0, "cm"),
          legend.key = element_rect(fill="transparent", color=NA),
          legend.key.size = unit(0.3, "lines"),
          legend.key.height = unit(0.5, "lines"),
          legend.key.width = unit(0.7, "lines"),
          legend.text = element_text(size=rel(0.6), colour="grey40"),
          legend.text.align = 0.5,
          legend.title = element_text(size=rel(0.7)),
          legend.title.align = 0,
          legend.position = "bottom",
          legend.direction = "vertical",
          legend.justification = "center",
          legend.box = "horizontal"
          )
}

plotSurvivingRowsLogical <- function(d,stage='platemap') {
  wells<-as.data.frame(t(d['det_well',]));
  platemap<-transform(wells,Row=as.numeric(match(toupper(substr(wells$det_well, 1, 1)), LETTERS)),Column=as.numeric(substr(wells$det_well, 2, 5)));
  ggplot(data=platemap, aes(x=Column, y=Row)) +
    geom_point(data=expand.grid(seq(1, 12), seq(1, 8)), aes(x=Var1, y=Var2),
               color="grey90", fill="white", shape=21, size=6) +
    geom_point(size=10) +
    coord_fixed(ratio=(13/12)/(9/8), xlim=c(0.5, 12.5), ylim=c(0.5, 8.5)) +
    scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8]) +
    scale_x_continuous(breaks=seq(1, 12)) +
    labs(title=stage) +
    theme_bdc_microtiter()
}

plotSurvivingRowsNumeric <- function(d,vals,stage='platemap') {
  wells<-as.data.frame(t(d['det_well',]));
  platemap<-transform(wells,Data=vals,Row=as.numeric(match(toupper(substr(wells$det_well, 1, 1)), LETTERS)),Column=as.numeric(substr(wells$det_well, 2, 5)));
  ggplot(data=platemap, aes(x=Column, y=Row)) +
    geom_point(data=expand.grid(seq(1, 12), seq(1, 8)), aes(x=Var1, y=Var2),
               color="grey90", fill="white", shape=21, size=6) +
    geom_point(size=10,aes(colour=Data)) +
    geom_text(aes(label=as.character(round(Data,2)),size=1)) +
    coord_fixed(ratio=(13/12)/(9/8), xlim=c(0.5, 12.5), ylim=c(0.5, 8.5)) +
    scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8]) +
    scale_x_continuous(breaks=seq(1, 12)) +
    labs(title=stage) +
    scale_color_gradientn(colours = c('blue','white','red')) +
    theme_bdc_microtiter()
}

plotSurvivingRowsNumericWithErrors <- function(d,vals,errs,stage='platemap') {
  wells<-as.data.frame(t(d['det_well',]));
  platemap<-transform(wells,Data=vals,MADs=errs,Row=as.numeric(match(toupper(substr(wells$det_well, 1, 1)), LETTERS)),Column=as.numeric(substr(wells$det_well, 2, 5)));
  ggplot(data=platemap, aes(x=Column, y=Row)) +
    geom_point(data=expand.grid(seq(1, 12), seq(1, 8)), aes(x=Var1, y=Var2),
               color="grey90", fill="white", shape=21, size=6) +
    geom_point(size=10,aes(colour=Data)) +
    geom_errorbarh(aes(xmax=Column + MADs,xmin=Column - MADs, height=0.2)) +
    geom_text(aes(label=as.character(round(Data,2)),size=1,vjust=2.5)) +
    coord_fixed(ratio=(13/12)/(9/8), xlim=c(0.5, 12.5), ylim=c(0.5, 8.5)) +
    scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8]) +
    scale_x_continuous(breaks=seq(1, 12)) +
    labs(title=stage) +
    scale_color_gradientn(limits=c(0,1),colours = c('blue','white','red')) +
    theme_bdc_microtiter()
}

#scale_color_gradientn(limits=c(-5,5),colours = c('black','blue','white','red','lightgreen')) +


simple_wrapper <- function (q) {
  plotSurvivingRowsLogical(q,stage='FromSimpleWrapper');
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
