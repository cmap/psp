#MAJOR VERSION 2.0
#MINOR VERSION 1
#04-MAR-2015


######################################
#                                    #
#    MASTER PROCESSING FUNCTIONS     #
#                                    #
######################################

P100processGCTMaster <- function (gctFileName=NULL,repAnnot=NULL, probeAnnot=NULL, dataTable=NULL,
                                  fileOutput=TRUE,outputFileName=NULL, processMode='full',
                                  optim=TRUE,log2=TRUE,samplePctCutoff=0.8, probePctCutoff=0.9, probeSDCutoff=3, distSDcutoff=3)
{
  ######################################################################################
  #  This function is the major entry point for automated data processing of P100.     #
  #  It supports two modes of data input and output.                                   #
  #  If the parameter 'gctFileName' is not NULL, it will look for a local file of that #
  #  name.                                                                             #
  #  If the paramters 'repAnnot', 'probeAnnot', and 'dataTable' are not NULL, it will  #
  #  assume that these are data objects from Panorama report views and proceed         #
  #  accordingly.                                                                      #
  #  Local users andor PANAORAMA should pass these parameters as named arguments!!!    #
  #  Also, now supporting processMode ('full','quick') switch for less redundancy      #
  ######################################################################################

  o<-list();
  
  #CHECK THAT AT LEAST ONE MODE OF INPUT IS FULLY EMPLOYED
  if (is.null(gctFileName) && (is.null(repAnnot) || is.null(repAnnot) || is.null(probeAnnot))) {
    stop('Either provide a gctFileName or data objects from Panorama');
  }
  
  #GET THE DATA (LIST) OBJECT AND SET METHOD SPECIFIC PARAMETERS
  if (!(is.null(gctFileName))) {
    #IF FROM LOCAL FILE
    o<-P100provideGCTlistObjectFromFile(gctFileName);
    if (is.null(outputFileName)) {
      outputFileName<-paste(gctFileName,'.processed.gct',sep='');
    }
  } else {
    #IF FROM PANORAMA DATA OBJECTS
    o<-P100provideGCTlistObjectFromPanorama(repAnnot, probeAnnot, dataTable);
    if(is.null(outputFileName))
    {
      outputFileName <- "${fileout:p100.processed.gct}";
    }
  }

  po<-list();
  #DO THE PROCESSING AND RETURN AN OBJECT READY TO BE WRITTEN TO GCT FILE FORMAT
  if (processMode == 'full') {
    po<-P100processGCT(o,optim=optim,log2=log2,samplePctCutoff=samplePctCutoff, probePctCutoff=probePctCutoff, probeSDCutoff=probeSDCutoff, distSDcutoff=distSDcutoff)
  } else if (processMode == 'quick') {
    po<-P100processGCTquick(o,log2=log2)
  } else {
    stop(paste('That processing mode is not supported: ',as.character(processMode),sep=''));
  } 

  #WRITE THE GCT FILE (IF FILE WRITING ENABLED)
  if (fileOutput) {
    P100writeGCTForProcessedObject(output.name=outputFileName,processedObject=po);
  }
}

GCPprocessGCTMaster <- function (gctFileName=NULL,repAnnot=NULL, probeAnnot=NULL, dataTable=NULL,
                                  fileOutput=TRUE,outputFileName=NULL, processMode='full', normalization_peptide_id = 'BI10052',
                                  log2=TRUE,samplePctCutoff=0.5, probePctCutoff=0.5, probeSDCutoff=4, probeGroupNormalization=FALSE)
{
  ######################################################################################
  #  This function is the major entry point for automated data processing of GCP.      #
  #  It supports two modes of data input and output.                                   #
  #  If the parameter 'gctFileName' is not NULL, it will look for a local file of that #
  #  name.                                                                             #
  #  If the paramters 'repAnnot', 'probeAnnot', and 'dataTable' are not NULL, it will  #
  #  assume that these are data objects from Panorama report views and proceed         #
  #  accordingly.                                                                      #
  #  Local users andor PANAORAMA should pass these parameters as named arguments!!!    #
  #  Also, now supporting processMode ('full','quick') switch for less redundancy      #
  ######################################################################################

  o<-list();
  
  #CHECK THAT AT LEAST ONE MODE OF INPUT IS FULLY EMPLOYED
  if (is.null(gctFileName) && (is.null(repAnnot) || is.null(repAnnot) || is.null(probeAnnot))) {
    stop('Either provide a gctFileName or data objects from Panorama');
  }
  
  #GET THE DATA (LIST) OBJECT AND SET METHOD SPECIFIC PARAMETERS
  if (!(is.null(gctFileName))) {
    #IF FROM LOCAL FILE
    o<-P100provideGCTlistObjectFromFile(gctFileName);
    if (is.null(outputFileName)) {
      outputFileName<-paste(gctFileName,'.processed.gct',sep='');
    }
  } else {
    #IF FROM PANORAMA DATA OBJECTS
    o<-P100provideGCTlistObjectFromPanorama(repAnnot, probeAnnot, dataTable);
    if(is.null(outputFileName))
    {
      outputFileName <- "${fileout:p100.processed.gct}";
    }
  }

  po<-list();
  #DO THE PROCESSING AND RETURN AN OBJECT READY TO BE WRITTEN TO GCT FILE FORMAT
  if (processMode == 'full') {
    po<-GCPprocessGCT(o,log2=log2,samplePctCutoff=samplePctCutoff, probePctCutoff=probePctCutoff, probeSDCutoff=probeSDCutoff, normalization_peptide_id=normalization_peptide_id, probeGroupNormalization=probeGroupNormalization)
  } else {
    stop(paste('That processing mode is not supported: ',as.character(processMode),sep=''));
  } 

  #WRITE THE GCT FILE (IF FILE WRITING ENABLED)
  if (fileOutput) {
    P100writeGCTForProcessedObject(output.name=outputFileName,processedObject=po);
  }
}

P100processGCT <- function (g,optim=TRUE,log2=TRUE,samplePctCutoff=0.8, probePctCutoff=0.9, probeSDCutoff=3, distSDcutoff=3) {

  static_headers<-g$static_headers;
  surviving_headers<-g$surviving_headers;
  surviving_rowAnnots<-g$surviving_rowAnnots;
  colnames(surviving_rowAnnots)<-static_headers[1,];
  dt<-g$dt;
  #gctFileName=g$gctFileName;

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

  q<-list(inputData=dt,initialSampleFiltering=s,probeFiltering=f,optimizationParams=o,secondarySampleCorrectionAndFiltering=b,
          normalizedData=n,outputData=n$normalizedData,static_headers=static_headers,surviving_headers=surviving_headers,
          surviving_rowAnnots=surviving_rowAnnots,colsAnnot=g$colsAnnot,rowsAnnot=g$rowsAnnot);
  return(q);
}

P100processGCTquick <- function (g,log2=TRUE) {

  static_headers<-g$static_headers;
  surviving_headers<-g$surviving_headers;
  surviving_rowAnnots<-g$surviving_rowAnnots;
  colnames(surviving_rowAnnots)<-static_headers[1,];
  dt<-g$dt;

  #log2 transform
  if (log2) {
    dt[dt==0]=NA;
    dt<-log(dt)/log(2);
    surviving_headers<-.updateProvenanceCode(static_headers,surviving_headers,"L2X");
  }

  n<-P100rowMedianNormalize(dt);
  surviving_headers<-.updateProvenanceCode(static_headers,surviving_headers,"RMN");

  q<-list(inputData=dt,
          normalizedData=n,outputData=n$normalizedData,static_headers=static_headers,surviving_headers=surviving_headers,
          surviving_rowAnnots=surviving_rowAnnots,colsAnnot=g$colsAnnot,rowsAnnot=g$rowsAnnot);
  return(q);
}

GCPprocessGCT <- function (g,log2=TRUE,samplePctCutoff=0.8, probePctCutoff=0.9, probeSDCutoff=3, normalization_peptide_id = 'BI10052', probeGroupNormalization=FALSE) {

  static_headers<-g$static_headers;
  surviving_headers<-g$surviving_headers;
  surviving_rowAnnots<-g$surviving_rowAnnots;
  colnames(surviving_rowAnnots)<-static_headers[1,];
  dt<-g$dt;
  #gctFileName=g$gctFileName;

  #log2 transform
  if (log2) {
    dt[dt==0]=NA;
    dt<-log(dt)/log(2);
    surviving_headers<-.updateProvenanceCode(static_headers,surviving_headers,"L2X");
  }

  h<-GCPhistoneNormalize(dt,surviving_rowAnnots,normalization_peptide_id)
  surviving_headers<-.updateProvenanceCode(static_headers,surviving_headers,"H3N");
  surviving_rowAnnots<-h$survivingRows;

  s<-P100filterSamplesForPoorCoverage(h$filteredData, pctFilterCutoff=samplePctCutoff)
  surviving_headers<-surviving_headers[,s$colsPassing];
  surviving_headers<-.updateProvenanceCode(static_headers,surviving_headers,"SF5");

  f<-P100filterProbes(s$filteredData, pctFilterCutoff=probePctCutoff, sdFilterCutoff=probeSDCutoff);
  surviving_rowAnnots<-surviving_rowAnnots[f$rowsPassing,];
  surviving_headers<-.updateProvenanceCode(static_headers,surviving_headers,"PF5");

  n<-P100rowMedianNormalize(f$filteredData);
  if (length(unique(surviving_rowAnnots$pr_probe_normalization_group)) > 1) {
    probeGroupNormalization<-TRUE;
  }

  if (probeGroupNormalization) {
    n<-GCPprobeGroupSpecificRowMedianNormalize(data=f$filteredData,ra=surviving_rowAnnots,sth=static_headers,sh=surviving_headers)
    surviving_headers<-.updateProvenanceCode(static_headers,surviving_headers,"GMN");
  } else {
    surviving_headers<-.updateProvenanceCode(static_headers,surviving_headers,"RMN");
  }

  q<-list(inputData=dt,histoneNormedData=h, initialSampleFiltering=s,probeFiltering=f,
          normalizedData=n,outputData=n$normalizedData,static_headers=static_headers,surviving_headers=surviving_headers,
          surviving_rowAnnots=surviving_rowAnnots,colsAnnot=g$colsAnnot,rowsAnnot=g$rowsAnnot);
  return(q);
}


######################################
#                                    #
#          FILE I/O SUPPORT          #
#                                    #
######################################

# Example GCT header: 
  # 
  # #1.3
  # 96  36  9 15
  # 
  # 96 = number of probes (data rows)
  # 36 = number of replicates (in columns)
  # 9 = number of probe annotations (+1 for the 'id' column -- first column)
  # 15 = number of replicate annotations (+1 for the 'id' row -- first row)
  # ---------------------------------------------------
  # |id|  probe_annots - 9   | replicate names - 36   |
  # ---------------------------------------------------
  # |r |                     |                        |
  # |e |                     |                        |
  # |p |                     |                        |
  # |- |                     |                        |
  # |a |    blank            |     replicate          |
  # |n |  (static_headers)   |  annotation values     |
  # |n |                     |  (surviving_headers)   |
  # |o |                     |                        |
  # |t |                     |                        |
  # |- |                     |                        |
  # |15|                     |                        |
  #----------------------------------------------------
  # |p |                     |                        |
  # |r |                     |                        |
  # |o |                     |                        |
  # |b |  probe annotation   |    data - 96 rows      |
  # |e |  values             |      (dt)              |
  # |s |(surviving_rowannots)|                        |
  # |  |                     |                        |
  #----------------------------------------------------

P100provideGCTlistObjectFromPanorama <- function (repAnnot, probeAnnot, dataTable) {
  local_surviving_headers<-repAnnot;
  local_dt<-dataTable;
  local_surviving_rowAnnots<-data.frame(id=rownames(probeAnnot),probeAnnot);
  local_static_headers<-as.data.frame(matrix(data='NA',ncol=dim(local_surviving_rowAnnots)[2],nrow=dim(local_surviving_headers)[1]),stringsAsFactors=FALSE);
  rownames(local_static_headers)<-rownames(local_surviving_headers);
  colnames(local_static_headers)<-colnames(local_surviving_rowAnnots);
  local_static_headers[,1]<-rownames(local_static_headers);
  local_static_headers[1,]<-colnames(local_static_headers);
  a<-list(surviving_headers=local_surviving_headers,static_headers=local_static_headers,surviving_rowAnnots=local_surviving_rowAnnots,dt=local_dt,colsAnnot=dim(local_static_headers)[2],rowsAnnot=dim(local_static_headers)[1],gctFileName='panorama');
  return(a)
}

P100provideGCTlistObjectFromFile <- function (gctFileName) {
  g<-P100_ReadGCT(gctFileName);
  local_static_headers<-g$headers[,1:g$colsAnnot];
  local_surviving_headers<-g$headers[,(g$colsAnnot+1):(g$colsAnnot+g$colsData)];
  local_surviving_rowAnnots<-g$rowAnnot;
  local_dt<-g$data;
  a<-list(surviving_headers=local_surviving_headers,static_headers=local_static_headers,surviving_rowAnnots=local_surviving_rowAnnots,dt=local_dt,colsAnnot=g$colsAnnot,rowsAnnot=g$rowsAnnot,gctFileName=gctFileName);
  return(a)
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

P100writeGCTForProcessedObject <- function (output.name,processedObject) {
  
  P100writeGCT(output.name, 
               processedObject$static_headers, 
               processedObject$surviving_headers, 
               processedObject$surviving_rowAnnots, 
               processedObject$outputData,
               processedObject$colsAnnot,
               processedObject$rowsAnnot
               );
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

######################################
#                                    #
#         DATA Q/C FUNCTIONS         #
#                                    #
######################################

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

GCPprobeGroupSpecificRowMedianNormalize <- function (data,ra,sth,sh) {
  pr_id_equiv<-'pr_id';
  if (!('pr_id' %in% rownames(sth))) {
    pr_id_equiv<-'id';
  }
  colnames(ra)<-sth[pr_id_equiv,];
  probe_normalization_assignments<-as.numeric(ra$pr_probe_normalization_group);
  probe_normalization_group<-unique(probe_normalization_assignments);
  num_probe_groups<-length(probe_normalization_group);
  sample_group_vectors<-sh['det_normalization_group_vector',];
  sample_group_matrix<-matrix(as.numeric(unlist(strsplit(unlist(sample_group_vectors),split=','))),nrow=num_probe_groups);
  sample_group_maxes<-apply(sample_group_matrix,1,max);
  copy_of_data<-data;
  G<-list();

  for (j in 1:num_probe_groups) {
    for (k in 1:sample_group_maxes[j]) {
      if (sample_group_maxes[j] > 1) {
        G[[j]]<-list();
      }
      working_data<-data[probe_normalization_assignments==j,sample_group_matrix[j,]==k];
      working_data_medians<-apply(working_data,1,median,na.rm=TRUE);
      working_data<-working_data-working_data_medians;
      copy_of_data[probe_normalization_assignments==j,sample_group_matrix[j,]==k]<-working_data;
      if (sample_group_maxes[j] > 1) {
        #print(c(j,k))
        #print(working_data_medians);
        G[[j]][[k]]<-working_data_medians;
        #print(G);
      } else {
        G[[j]]<-working_data_medians;
      }
    }
  }
  q<-list(originalData=data,normalizedData=copy_of_data,pna=probe_normalization_assignments,sgm=sample_group_matrix,rowMedians=G);
  return(q);
}

P100filterSamplesForPoorCoverage <- function (dt, pctFilterCutoff=0.8) {
  #Filter out columns missing x% of data;
  retList<-list(colsPassing=NULL,filteredData=dt,originalData=dt);
  dt<-as.matrix(dt);
  ddt<-dim(dt);
  nRows<-ddt[1];
  nCols<-ddt[2];

  numNAs<-.countColNAs(dt);
  pctNotNAs<-1-(numNAs/nRows);
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
  pctNotNAs<-1-(numNAs/nCols);
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

GCPhistoneNormalize <- function (dt, surv_rows, norm_peptide_id) {
  #Filter out columns missing x% of data;
  retList<-list(filteredData=dt,originalData=dt,normConstants=dt[1,],survivingRows=surv_rows);
  dt<-as.matrix(dt);
  ddt<-dim(dt);
  nRows<-ddt[1];
  nCols<-ddt[2];

  #Find the norm peptide row:
  indarr<-surv_rows[,1]==norm_peptide_id;
  norm_consts<-dt[indarr,];
  retList$normConstants<-norm_consts;
  norm_mat<-matrix(data=retList$normConstants,nrow=nRows,ncol=nCols,byrow=TRUE)
  retList$filteredData<-dt-norm_mat;
  retList$filteredData<-retList$filteredData[(!indarr),];
  retList$survivingRows<-retList$survivingRows[(!indarr),];

  return(retList);

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

genericRowNormalizeByColumnAnnot<-function(dt,columnAnnot,groupingField,subset=NULL) {
  retList<-list(groupRowMedians=NULL,normalizedData=dt,originalData=dt);
  normGroups<- unique(unlist(columnAnnot[groupingField,]));
  grm<-matrix(data=FALSE,nrow=length(dt[,1]),ncol=length(normGroups));
  colnames(grm)<-normGroups;
  for (x in 1:length(normGroups)) {
    dtTemp<-dt[,columnAnnot[groupingField,]==normGroups[x]];
    groupNorm<-apply(dtTemp,1,median,na.rm=TRUE);
    grm[,x]<-groupNorm;
    dtTemp<-dtTemp-groupNorm;
    dt[,columnAnnot[groupingField,]==normGroups[x]]<-dtTemp;
  }
  retList$groupRowMedians<-grm;
  retList$normalizedData<-dt;
  return(retList);
}

.updateProvenanceCode <- function (static_annot,annot,code) {
  indarr<-static_annot[,1]=="provenance_code";
  annot[indarr,]<-paste(annot[indarr,],code,sep="+");
  return(annot);
}


######################################
#                                    #
#   PLATE MAP EVALUATION FUNCTIONS   #
#                                    #
######################################


P100evaluateGCTasPlate <- function (gctFileName,optim=TRUE,fileOutput=TRUE,log2=TRUE,samplePctCutoff=0.8, probePctCutoff=0.9, probeSDCutoff=3, distSDcutoff=3) {
  ###THIS FUNCTION PROBABLY NEEDS REWRITING
  ###TRY TO USE THE UNIVERSAL P100P

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

P100evaluateEnrichment <- function (gctFileName=NULL,repAnnot=NULL, probeAnnot=NULL, dataTable=NULL, all=FALSE) {
  o<-list();
  
  if (is.null(gctFileName) && (is.null(repAnnot) || is.null(repAnnot) || is.null(probeAnnot))) {
    stop('Either provide a gctFileName or data objects from Panorama');
  }
  
  if (!(is.null(gctFileName))) {
    o<-P100provideGCTlistObjectFromFile(gctFileName);
  } else {
    o<-P100provideGCTlistObjectFromPanorama(repAnnot, probeAnnot, dataTable);
  }

  so<-P100scaleGCTdataTable(o);
  cso<-P100collapseProbesToSingleValues(so);
  if (!(all)) {
    plotSurvivingRowsNumericWithErrorsFromCollapsedPlateListObject(cso,stage='Enrichment Map');
  } else {
    pl<-list();
    pl[[1]]<-plotSurvivingRowsNumericWithErrorsFromCollapsedPlateListObject(cso,stage='Enrichment Map');
    for (j in 1:dim(so$surviving_rowAnnots)[1]) {
      pl[[j+1]]<-plotAnIndexedRowNumericFromPlateListObject(so,j,stage=so$surviving_rowAnnots[j,1]);
    }
    multiplot(plotlist=pl,cols=4);
  }
}

P100evaluateProbeIntensities <- function (gctFileName=NULL,repAnnot=NULL, probeAnnot=NULL, dataTable=NULL, all=FALSE, labelType='heavy') {
  o<-list();
  
  if (is.null(gctFileName) && (is.null(repAnnot) || is.null(repAnnot) || is.null(probeAnnot))) {
    stop('Either provide a gctFileName or data objects from Panorama');
  }
  
  if (!(is.null(gctFileName))) {
    o<-P100provideGCTlistObjectFromFile(gctFileName);
  } else {
    o<-P100provideGCTlistObjectFromPanorama(repAnnot, probeAnnot, dataTable);
  }

  colnames(o$surviving_rowAnnots)<-o$static_headers['pr_id',];
  fo<-o;
  fo$surviving_rowAnnots<-o$surviving_rowAnnots[o$surviving_rowAnnots$IsotopeLabelType==labelType,]
  fo$dt<-o$dt[o$surviving_rowAnnots$IsotopeLabelType==labelType,]

  so<-P100scaleGCTdataTable(fo);
  cso<-P100collapseProbesToSingleValues(so);
  if (!(all)) {
    plotSurvivingRowsNumericWithErrorsFromCollapsedPlateListObject(cso,stage='Heavy Intensity Map');
  } else {
    pl<-list();
    pl[[1]]<-plotSurvivingRowsNumericWithErrorsFromCollapsedPlateListObject(cso,stage='Heavy Intensity Map');
    for (j in 1:dim(so$surviving_rowAnnots)[1]) {
      pl[[j+1]]<-plotAnIndexedRowNumericFromPlateListObject(so,j,stage=so$surviving_rowAnnots[j,1]);
    }
    multiplot(plotlist=pl,cols=4);
  }

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

P100scaleGCTdataTable <- function (o,method='zeroToMax',by='row') {
  #o is a GCT list object
  MARGIN<-1;
  if (by == 'column') {
    MARGIN<-2;
  }
  o$dt[is.na(o$dt)]<-0;
  mx<-apply(o$dt,MARGIN,max,na.rm=TRUE);
  if (method=='zeroToMax') {
    o$dt<-o$dt/mx;
  } else {
    #NOTE: IMPLEMENT OTHER METHODS LATER AS ELSIFs
    #do nothing
  }
  return(o);
}

P100collapseProbesToSingleValues <- function (o, method='mean',by='column') {
  MARGIN<-1;
  if (by == 'column') {
    MARGIN<-2;
  }
  cfunc<-mean;
  efunc<-sd;
  if (method == 'median') {
    cfunc<-median;
    efunc<-mad;
  }
  o$sv<-apply(o$dt,MARGIN,cfunc,na.rm=TRUE);
  o$err<-apply(o$dt,MARGIN,efunc,na.rm=TRUE);
  return(o);
}

plotSurvivingRowsNumericWithErrorsFromCollapsedPlateListObject <- function(o,stage='platemap') {
  wells<-as.data.frame(t(o$surviving_headers['det_well',]));
  platemap<-transform(wells,Data=o$sv,MADs=o$err,Row=as.numeric(match(toupper(substr(wells$det_well, 1, 1)), LETTERS)),Column=as.numeric(substr(wells$det_well, 2, 5)));
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

plotAnIndexedRowNumericFromPlateListObject <- function(o,ind,stage='platemap') {
  wells<-as.data.frame(t(o$surviving_headers['det_well',]));
  platemap<-transform(wells,Data=as.numeric(o$dt[ind,]),Row=as.numeric(match(toupper(substr(wells$det_well, 1, 1)), LETTERS)),Column=as.numeric(substr(wells$det_well, 2, 5)));
  ggplot(data=platemap, aes(x=Column, y=Row)) +
    geom_point(data=expand.grid(seq(1, 12), seq(1, 8)), aes(x=Var1, y=Var2),
               color="grey90", fill="white", shape=21, size=6) +
    geom_point(size=10,aes(colour=Data)) +
    geom_text(aes(label=as.character(round(Data,2)),size=1,vjust=2.5)) +
    coord_fixed(ratio=(13/12)/(9/8), xlim=c(0.5, 12.5), ylim=c(0.5, 8.5)) +
    scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8]) +
    scale_x_continuous(breaks=seq(1, 12)) +
    labs(title=stage) +
    scale_color_gradientn(limits=c(0,1),colours = c('blue','white','red')) +
    theme_bdc_microtiter()
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
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


#######################################
#                                     #
# IMPUTATION/INFERENCE PREP FUNCTIONS #
#                                     #
#######################################

#VERY BETA - CONSIDER REMOVING FROM PRODUCTION CODE
#UNSURE OF USAGE
#WOULD FEED INTO PYTHON CODE

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

#######################################
#                                     #
#    DATA MINING / MISC FUNCTIONS     #
#                                     #
#######################################

#THESE FUNCTIONS ARE OF UNKNOWN PROVENANCE
#LOOKS LIKE THEY REQUIRE SOME EXTERNAL LIBRARIES

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


