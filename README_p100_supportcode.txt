P100writeGCT('good_by_eye_working_corrected_targeted_data.gct',r$static_headers,gbe_survhead,r$surviving_rowAnnots,gbe_working_targeted_data,10,17)

P100processGCT('gbe_working_corrected_targeted_data.gct',log2=FALSE,samplePctCutoff=0,probePctCutoff=0.8)


perts<-c(
"PD0325901",
"PD0325901",
"PD0325901",
"PD0325901",
"WYE125132",
"WYE125132",
"WYE125132",
"PD0325901",
"PD0325901",
"PD0325901",
"PD0325901",
"PD0325901",
"PD0325901",
"WYE125132",
"WYE125132"
)


diffmat <- function (x) {
  l<-length(x);
  o<-matrix(nrow=l,ncol=l);
  for (i in 1:l) {
    for (j in 1:l) {
      o[i,j]=x[i]-x[j];
    }
  }
  rownames(o)<-rownames(x);
  colnames(o)<-rownames(x);
  return(o);
}

sna <- function (x) {
  sum(is.na(x));
}

compute_differences <- function (x) {
  p<-dim(x);
  nam<-rownames(x);
  resarr<-array(dim=c(p[1],p[1],p[2]));
  for (k in 1:p[2]) {
    resarr[,,k]<-diffmat(x[,k]);
  }
  return(resarr);
}

#for each probe, pick top 5 (or n) highest correlated non-self probes that EXIST IN DATASET, triangulate difference based on times of known probes?
# straight average? weighted average (by corr, by distance, by distance variation)?

#input is a probes list, and model components (as an object)?

#get all non-na-data
#foreach na in scheduling run

testprobe<-'8982_EMK1_T597_ST[+80]FHAGQLR'
allgood.meanrts[,testprobe]
allgood.meanrts[testprobe]
allgood.cormatrix[,testprobe]
sort(allgood.cormatrix[,testprobe])
help(sort)
sort(allgood.cormatrix[,testprobe],decreasing=TRUE)
sort(allgood.cormatrix[,testprobe],decreasing=TRUE)[1:10]
testtop10<-sort(allgood.cormatrix[,testprobe],decreasing=TRUE)[2:10]
testtop10<-sort(allgood.cormatrix[,testprobe],decreasing=TRUE)[1:10]
testset<-sort(allgood.cormatrix[,testprobe],decreasing=TRUE)[2:10]
allgood.meanrts[testtop10]
testtop10<-names(sort(allgood.cormatrix[,testprobe],decreasing=TRUE)[1:10])
testset<-names(sort(allgood.cormatrix[,testprobe],decreasing=TRUE)[2:10])
allgood.meanrts[testtop10]
testprobe
allgood.meanrts[testset]
allgood.meandiffs[testset,testprobe]
allgood.diffvars[testset,testprobe]
dim(allgood)
allgood[,1]
rownames(allgood)
testcase.1<-allgood[,1]
names(testcase.1)<-rownames(allgood)
testcase.1[testset]
testcase.1[testset]-allgood.meandiffs[testset,testprobe]
mean(testcase.1[testset]-allgood.meandiffs[testset,testprobe])
testcase.1[testprobe]



testcase.2<-t_rt[,3]
names(tescase.2)<-rownames(t_rt)
names(testcase.2)<-rownames(t_rt)
is.na(testcase.2)
is.na(testcase.2)==TRUE
testcase[is.na(testcase.2)==TRUE]
testcase.2[is.na(testcase.2)==TRUE]
names(testcase.2[is.na(testcase.2)==TRUE])
probes_to_get<-names(testcase.2[is.na(testcase.2)==TRUE])
probes_to_use<-names(testcase.2[is.na(testcase.2)==FALSE])
history(500)
probes_to_get[1]
allgood.cormatrix[probes_to_use,probes_to_get]
relevant_cors<-allgood.cormatrix[probes_to_use,probes_to_get]
relevant_dists<-allgood.meandiffs[probes_to_use,probes_to_get]
relevant_dists
sort(relevant_cors[,1],descending=TRUE)
sort(relevant_cors[,1],decreasing=TRUE)
topcors.1<-names(sort(relevant_cors[,1],decreasing=TRUE))[2:9]
relevant_dists[topcors.1,1]
testcase.2
testcase.2[topcors]-relevant_dists[topcors.1,1]
testcase.2[topcors.1]-relevant_dists[topcors.1,1]
mean(testcase.2[topcors.1]-relevant_dists[topcors.1,1])
var(testcase.2[topcors.1]-relevant_dists[topcors.1,1])
var(testcase.2[topcors.1]-relevant_dists[topcors.1,1])**2

    #consider enforcing minimum correlation?

P100reschedule <- function (measured,model,topN=9) {
  #measured is a vector containing some NAs with names as P100 pr_id
  probes_to_get<-names(measured[is.na(measured)==TRUE]);
  probes_to_use<-names(measured[is.na(measured)==FALSE]);
  scheduled<-measured;
  for (n in 1:length(probes_to_get)) {
    relevant_cors<-model$cormat[probes_to_use,probes_to_get[n]];
    relevant_dists<-model$diffmat[probes_to_use,probes_to_get[n]];
    sorted_cors<-sort(relevant_cors[,1],decreasing=TRUE);
    best_cors<-names(sorted_cors[2:topN]);
    calcRT.mean<-mean(measured[best_cors]-relevant_dists[best_cors,probes_to_get[n]);
    calcRT.var<-var(measured[best_cors]-relevant_dists[best_cors,probes_to_get[n]);
    scheduled[probes_to_get[n]]=calcRT.mean;
  }
  return(scheduled);
}


schedModel<-list(cormat=allgood.cormatrix,diffmat=allgood.meandiffs,varmat=allgood.diffvars,sourcedata=allgood,alldiffs=allgood.diffs,allvars=allgood.var,meanrts=allgood.meanrts)


for (j in 1:17) {
	tmat<-blankMatrix[,j]
	names(tmat)<-rownames(blankMatrix)
	plot(schedModel$sourcedata[,j],P100reschedule(tmat,schedModel,topN=3)-schedModel$sourcedata[,j])
}


leave_out = c(8,16,24,32,40,48,56,64,72,80);
top_n = c(1,2,3,4,5,6,7,8,9);
a<-b<-c<-array(dim=length(leave_out)*length(top_n));
counter = 0;
for (i in 1:length(leave_out)) {
  for (j in 1:length(top_n)) {
    counter=counter+1;
    residuals<-vector(mode='numeric',length=96);
    for (nreps in 1:20) {
      temp<-rs02.times2;
      temp[sample(1:96,leave_out[i])]<-NA;
      residuals=residuals+P100reschedule(temp,schedModel,topN=j)-rs02.times2;
    }
    residuals=residuals/100;
    a[counter]<-leave_out[i];
    b[counter]<-top_n[j];
    c[counter]<-sqrt(mean(residuals**2));
  }
}



#######################
Genetic ALGORITHM STUFF
#######################

library('GA')

response<-t(bs.targ[1,c(-1,-2)])
modvars<-t(bs.dia[bs.dia$pr_id=='7849_NANS_S275_ALGS[+80]PTKQLLPC[+57]EMAC[+57]NEK',c(-1,-2)])
colnames(modvars)<-bs.dia[bs.dia$pr_id=='7849_NANS_S275_ALGS[+80]PTKQLLPC[+57]EMAC[+57]NEK',2]
modvars
response
colnames(response)<-'targ'
response
testdf<-data.frame(response,1/modvars)
testdf<-log(testdf)/log(2)

dft<-function(respdat,preddat) {
  predmean<-apply(preddat,1,mean)
  return(sqrt(sum((respdat-predmean)**2)))
}

apply(testdf[,-1],1,mean)
apply(testdf[,-1],1,mean,na.rm=TRUE)

dft<-function(respdat,preddat) {
  predmean<-apply(preddat,1,mean,na.rm=TRUE)
  return(sqrt(sum((respdat-predmean)**2)))
}

pddm<-apply(testdf[,-1],1,mean,na.rm=TRUE)
pddm
testdf[,1]-pddm
(testdf[,1]-pddm)**2
sum((testdf[,1]-pddm)**2)
sqrt(sum((testdf[,1]-pddm)**2))
dft(testdf[,1],testdf[,-1])

fitness <- function(string) {
  inc<-which(string==1)
  X<-cbind(1,QQ[,inc])
  resdist<-dft(PP,X)
  return(1/resdist)
}

QQ<-testdf[,-1]
PP<-testdf[,1]
QQ
PP
GA<-ga("binary",fitness=fitness,nBits=ncol(QQ),names=colnames(QQ),monitor=plot)
dev.new()
GA<-ga("binary",fitness=fitness,nBits=ncol(QQ),names=colnames(QQ),monitor=plot)
dev.new()
plot(GA)
summary(GA)

colnames(bs.dia)

bs.areadia<-read.csv('BestSamples_DIA_transition_simple_area.csv',na.strings='#N/A')
head(bs.areadia)
therats<-bs.areadia[,3:14]/bs.areadia[,15:26]
head(therats)
bs.areadia[,3:13]
bs.areadia[1:6,3:14]
sum(bs.areadia[1:6,3:14])
apply(bs.areadia[1:6,3:14],2,sum,na.rm=TRUE)
apply(bs.areadia[1:6,3:14],2,sum,na.rm=TRUE)/apply(bs.areadia[1:6,15:26],2,sum,na.rm=TRUE)
1/(apply(bs.areadia[1:6,3:14],2,sum,na.rm=TRUE)/apply(bs.areadia[1:6,15:26],2,sum,na.rm=TRUE))
1/(apply(bs.areadia[c(1,4,5),3:14],2,sum,na.rm=TRUE)/apply(bs.areadia[c(1,4,5),15:26],2,sum,na.rm=TRUE))
bs.areadia[1:6,3:14]
bs.areadia[1:6,3:14]==0
bs.areadia[1:6,3:14]==NA
bs.areadia[1:6,3:14]==0
is.na(bs.areadia[1:6,3:14])
sur.1<-bs.areadia[1:6,3:14]
sur.2<-bs.areadia[1:6,15:26]
sur.1
sur.2
sur.1==0
sur.1[sur.1==0]
sur.1[sur.1==0]=NA
is.na(sur.1)
sur.2[is.na(sur.1)]=NA
sur.2
apply(sur.2,2,sum,na.rm=TRUE)/apply(sur.1,2,sum,na.rm=TRUE)

----
HOW TO DO GENETIC ALGORITHM (requires GA package)
-----

ahnak.csv<-read.csv('AHNAK_CL01_DIA_simple.csv')
dim(ahnak.csv)
ahnak.areas<-t(ahnak.csv[,3:26])
colnames(ahnak.areas)<-ahnak.csv[,2]
predictor.order<-rownames(potential.predictors)
ahnak.ordered<-ahnak.areas[predictor.order,]

measured.values<-bs.targ[2,3:14]
measured.values<-log(measured.values)/log(2)
measured.values<-unlist(measured.values)

trivial.transitions<-apply(ahnak.ordered,2,sum,na.rm=TRUE)==0
ahnak.nontrivial<-ahnak.ordered[,!(trivial.transitions)]
potential.predictors<-ahnak.nontrivial
ahnak.GA<-ga("binary",fitness=fitness,nBits=ncol(potential.predictors),names=colnames(potential.predictors),monitor=plot, maxiter=200)
summary(ahnak.GA)
ahnak.solutions<-ahnak.GA@solution
apply(ahnak.solutions,1,sum)
#mininum number of tx occurs in row 6 (20 tx)
ahnak.minimal.transitions<-names(which(ahnak.solutions[6,]==1))
ahnak.dia_calculated<-log(P100replicateSkylineRatio(ahnak.nontrivial[,which(ahnak.solutions[6,]==1)]))/log(2)
plot(ahnak.measured.values ~ ahnak.dia_calculated)
ahnak.lm<-lm(ahnak.measured.values ~ ahnak.dia_calculated)
summary(ahnak.lm)

-with jenn-

ourdata<-read.csv('PFKF_p100_transition_area_report.csv')
head(ourdata)
pfkf.csv<-fix_skyline_csv(ourdata)
head(pfkf.csv)
dim(pfkf.csv)
pfkf.areas<-t(pfkf.csv[,3:26])
colnames(pfkf.areas)<-pfkf.csv[,2]
colnames(pfkf.csv)
colnames(pfkf.areas)<-pfkf.csv[,'Fragment.Ion']
head(pfkf.areas)
predictor.order
pfkf.ordered<-pfkf.areas[predictor.order,]
rownames(pfkf.ordered)
head(bs.targ)
measured.values<-bs.targ[5,3:14]
colnames(measured.values)
rownames(pfkf.ordered)
measured.values<-log(measured.values)/log(2)
measured.values<-unlist(measured.values)
measured.values
trivial.transitions<-apply(pfkf.ordered,2,sum,na.rm=TRUE)==0
trivial.transitions
pfkf.nontrivial<-pfkf.ordered[,!(trivial.transitions)]
pfkf.nontrivial
potential.predictors<-pfkf.nontrivial
pfkf.GA<-ga("binary",fitness=fitness,nBits=ncol(potential.predictors),names=colnames(potential.predictors),monitor=plot, maxiter=200)
summary(pfkf.GA)
pfkf.solutions<-pfkf.GA@solution



----
support functions:
----

P100replicateSkylineRatio <- function (intensities) {
  #assume top half is light and bottom half is heavies
  dat.size<-dim(intensities);
  light.ind<-1:(dat.size[1]/2);
  heavy.ind<-(dat.size[1]/2+1):(dat.size[1]);
  dat.light<-intensities[light.ind,];
  dat.heavy<-intensities[heavy.ind,];

  #exlude all zero values;
  dat.light.zero.ind<-dat.light==0;
  dat.heavy.zero.ind<-dat.heavy==0;
  dat.light[dat.light.zero.ind]=NA;
  dat.heavy[dat.light.zero.ind]=NA;
  dat.light[dat.heavy.zero.ind]=NA;
  dat.heavy[dat.heavy.zero.ind]=NA;

  #compute ratios
  ret.ratios<-apply(dat.light,1,sum,na.rm=TRUE)/apply(dat.heavy,1,sum,na.rm=TRUE);
  names(ret.ratios)<-rownames(dat.light);
  return(ret.ratios);

}

fitness <- function(string) {
  inc<-which(string==1)
  subset.of.predictors<-cbind(1,potential.predictors[,inc])
  resdist<-dft_use_int(measured.values,subset.of.predictors);
  return(1/resdist);
}

dft_use_int <- function(respdat,preddat) {
  computed.ratio<-P100replicateSkylineRatio(preddat);
  computed.ratio<-log(computed.ratio)/log(2);
  return(sqrt(sum((respdat-computed.ratio)**2)))
}

fix_skyline_csv <- function (csvfile) {
  newcolumn<-paste(csvfile[,3],csvfile[,4],sep="_");
  csvfile[,3]<-newcolumn;
  csv.dims<-dim(csvfile);
  newcsv<-csvfile[,c(1,3,5:csv.dims[2])];
  return(newcsv);
}