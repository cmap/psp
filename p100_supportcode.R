#CONTAINS SUPPORT FUNCTIONS FOR RESCHEDULING/TRANSITION PICKING
#MAJOR VERSION 1.0
#MINOR VERSION 1
#21-JAN-2015


P100reschedule <- function (measured,model,topN=9) {
  #measured is a vector containing some NAs with names as P100 pr_id
  #the schedule model object can be found in:
  #P:\Projects\PhosphoLINCS\03_ASSAY_PIPELINE\DATASETS\08_LINCS_Transcenter\Outputs\P100schedulingData.RData
  #where P: is \\argon-cifs\prot_proteomics
  #OR in linux /prot/proteomics
  probes_to_get<-names(measured[is.na(measured)==TRUE]);
  probes_to_use<-names(measured[is.na(measured)==FALSE]);
  relevant_cors<-model$cormat[probes_to_use,probes_to_get];
  relevant_dists<-model$diffmat[probes_to_use,probes_to_get];
  scheduled<-measured;
  for (n in 1:length(probes_to_get)) {
    sorted_cors<-sort(relevant_cors[,probes_to_get[n]],decreasing=TRUE);
    best_cors<-names(sorted_cors[2:topN]);
    calcRT.mean<-mean(measured[best_cors]-relevant_dists[best_cors,probes_to_get[n]]);
    calcRT.var<-var(measured[best_cors]-relevant_dists[best_cors,probes_to_get[n]]);
    scheduled[probes_to_get[n]]=calcRT.mean;
  }
  return(scheduled);
}

P100correlateDIATxToTargeted <- function (targ.data,dia.data) {
  #targ.data and dia.data are data frames
  #Col 1 is probe id (pr_id)
  #Col 2 is transition type (called 'Fragment.Ion'; yn or bn for DIA, or 'proberatio' for targeted)
  #remaining columns are ratios in each sample
  #targ.data has one line per probe
  #dia.data has as many lines per probe as there are transitions
  #IMPORTANT NOTE: Targeted data is L/H; DIA data is H/L so must flip! (done during log2 step below)
  correlation.list<-list();
  probes.all<-as.character(unlist(targ.data$pr_id));
  for (j in 1:length(probes.all)) {
    probe.name<-probes.all[j];
    targ.local<-as.numeric(targ.data[targ.data$pr_id==probe.name,3:dim(targ.data)[2]]);
    dia.local<-as.data.frame(dia.data[dia.data$pr_id==probe.name,3:dim(dia.data)[2]],row.names=as.character(dia.data[dia.data$pr_id==probe.name,2]));
    #NOT SURE ABOUT LOG2 here...leave in for now
    targ.local<-log(targ.local)/log(2);
    dia.local<-log(1/dia.local)/log(2);
    cormat.local<-cor(targ.local,t(dia.local),use='pair');
    corvec.local<-as.vector(cormat.local);
    names(corvec.local)<-colnames(cormat.local);
    print(probe.name);
    print(dia.local);
    num.dia.tx<-dim(dia.local)[1];
    #par(mfrow=c(3,2));
    m<-matrix(c(1,1,2,3,4,5,6,7),ncol=2,byrow=TRUE);
    layout(m,widths = c(0.5,0.5), heights = c(0.04,0.32,0.32,0.32));
    par(mar = c(0,0,0,0));
    plot(1,1,type = "n",frame.plot = FALSE,axes = FALSE)
    u <- par("usr")
    text(1,u[4],labels = probe.name,col = "red",pos = 1,cex=2)
    par(mar = c(2,4,2,2) + 0.1)
    n.graphs<-min(6,num.dia.tx);
    tx.order<-names(sort(corvec.local,decreasing=TRUE));
    #print(tx.order);
    for (k in 1:n.graphs) {
      #print(tx.order[k]);
      if (sum(is.na(dia.local[tx.order[k],])) == length(dia.local[tx.order[k],])) {
         next;
      }
      plot(targ.local,dia.local[tx.order[k],],xlab="Targeted Ratio",ylab="DIA Transition Ratio",main=paste(tx.order[k],as.character(format(corvec.local[tx.order[k]],digits=3))));
    }
    #print(cormat.local);
    #print(sort(corvec.local,decreasing=TRUE));
    readline(prompt="Press enter.");
    correlation.list[[probe.name]]<-sort(corvec.local,decreasing=TRUE);
  }
  return(correlation.list);
}

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
  ret.ratios<-apply(dat.light,1,function(x){sum(as.numeric(x),na.rm=TRUE)})/apply(dat.heavy,1,function(x){sum(as.numeric(x),na.rm=TRUE)});
  names(ret.ratios)<-rownames(dat.light);
  return(ret.ratios);

}

dft<-function(respdat,preddat) {
  predmean<-apply(preddat,1,mean,na.rm=TRUE)
  return(sqrt(sum((respdat-predmean)**2,na.rm=TRUE)))
}

dft_use_int<-function(respdat,preddat) {
  computed.ratio<-P100replicateSkylineRatio(preddat);
  computed.ratio<-log(computed.ratio)/log(2);
  return(sqrt(sum((respdat-computed.ratio)**2,na.rm=TRUE)))
}

fitness <- function(string) {
  inc<-which(string==1)
  subset.of.predictors<-cbind(1,potential.predictors[,inc])
  resdist<-dft_use_int(measured.values,subset.of.predictors);
  return(1/resdist);
}

fitness_with_params <- function(string,pot_pred,meas_val) {
  inc<-which(string==1)
  subset.of.predictors<-cbind(1,pot_pred[,inc])
  resdist<-dft_use_int(log(meas_val)/log(2),subset.of.predictors);
  return(1/resdist);
}

P100selectBestDIATx <- function(diatxfilename,targvals) {
  ###NOT FINISHED###
  robj = list();
  robj$csv<-read.csv(diatxfilename);
  robj$csv.areacols<-grep('Area',colnames(robj$csv));
  robj$ionlabels<-paste(robj$csv$Fragment.Ion,robj$csv$Product.Charge,sep="_");
  ###NOT FINISHED###

  return(robj);
}

P100pickDIAtransitions <- function(probe_name_code,dia.datafile,targeted.datafile,replicatename.translationtablefile) {
  #read files
  dia.csv<-fix_skyline_csv(read.csv(dia.datafile,na="#N/A"));
  targeted.csv<-read.csv(targeted.datafile,na="#N/A");
  replicatenames.csv<-read.csv(replicatename.translationtablefile);

  #extract data and set order
  dia.probedata<-dia.csv[dia.csv[,1]==probe_name_code,]
  targeted.probedata<-targeted.csv[targeted.csv[,1]==probe_name_code,]
  dia.dim<-dim(dia.probedata)
  targeted.dim<-dim(targeted.probedata)
  dia.data<-t(dia.probedata[,3:dia.dim[2]]);
  colnames(dia.data)<-dia.probedata[,2];
  dia.data.lines<-dim(dia.data)[1];
  dia.light.data<-dia.data[as.character(unlist(replicatenames.csv['dia.light.names'])),];
  dia.heavy.data<-dia.data[as.character(unlist(replicatenames.csv['dia.heavy.names'])),];
  dia.data[1:(dia.data.lines/2),]<-dia.light.data;
  dia.data[(dia.data.lines/2+1):dia.data.lines,]<-dia.heavy.data;
  targeted.data<-t(targeted.probedata[,3:targeted.dim[2]]);
  targeted.data<-targeted.data[as.character(unlist(replicatenames.csv['targeted.names'])),];

  #
  targeted.data.log2<-log(targeted.data)/log(2);
  dia.data.prefit.log2<-log(P100replicateSkylineRatio(dia.data))/log(2);
  linearmodel.prefit <- lm(dia.data.prefit.log2 ~ targeted.data.log2);

  #make model on subset of data, test on remainder
  model_fraction<-2/3;
  total_samples<-length(targeted.data);
  model_samples<-round(model_fraction * total_samples);
  sampling_indices<-sample(1:total_samples,model_samples);
  dia_sampling_indices<-c(sampling_indices,sampling_indices+total_samples);
  targeted.data.model.basis<-targeted.data[sampling_indices];
  dia.data.model.basis<-dia.data[dia_sampling_indices,];
  targeted.data.holdout<-targeted.data[-sampling_indices];
  dia.data.holdout<-dia.data[-dia_sampling_indices,];

  #
  trivial.transitions<-apply(dia.data.model.basis,2,function(x){sum(as.numeric(x),na.rm=TRUE)})==0;
  dia.data.model.basis.nontrivial<-dia.data.model.basis[,!(trivial.transitions)];
  potential.predictors<-dia.data.model.basis.nontrivial;
  probe.GA<-ga("binary",fitness=fitness_with_params,pot_pred=potential.predictors,meas_val=targeted.data.model.basis,nBits=ncol(potential.predictors),names=colnames(potential.predictors), maxiter=500,run=50, popSize=500, pcrossover=0.9, pmutation=0.5, elitism=0.1)
  best_fitness_value<-probe.GA@fitnessValue;
  #for (j in 2:2) {
  #  print(paste('ITERATION: ',j,sep=''));
  #  temp.GA<-ga("binary",fitness=fitness_with_params,pot_pred=potential.predictors,meas_val=targeted.data.model.basis,nBits=ncol(potential.predictors),names=colnames(potential.predictors), maxiter=500,run=50, popSize=500, pcrossover=0.9, pmutation=0.5, elitism=0.1)
  #  if (temp.GA@fitnessValue > best_fitness_value) {
  #     best_fitness_value<-temp.GA@fitnessValue;
  #     probe.GA<-temp.GA;
  #  }
  #}
  probe.solutions<-probe.GA@solution
  probe.solutions.txcounts<-apply(probe.solutions,1,sum)
  first_best_solution<-which(probe.solutions.txcounts==min(probe.solutions.txcounts))[1]
  probe.minimal.transition.indices<-which(probe.solutions[first_best_solution,]==1);
  probe.minimal.transitions<-names(probe.minimal.transition.indices);

  #linear model on all / modeled / holdout,plots
  dia.data.modeled.recalculated.log2<-log(P100replicateSkylineRatio(dia.data[,probe.minimal.transitions]))/log(2);
  dia.training.data.recalculated.log2<-log(P100replicateSkylineRatio(dia.data.model.basis[,probe.minimal.transitions]))/log(2);
  dia.holdout.data.recalculated.log2<-log(P100replicateSkylineRatio(dia.data.holdout[,probe.minimal.transitions]))/log(2);
  targeted.training.data.log2<-log(targeted.data.model.basis)/log(2);
  targeted.holdout.data.log2<-log(targeted.data.holdout)/log(2);

  linearmodel.postfit.all <- lm(dia.data.modeled.recalculated.log2 ~ targeted.data.log2);
  linearmodel.postfit.training <- lm(dia.training.data.recalculated.log2 ~ targeted.training.data.log2);
  linearmodel.postfit.holdout <- lm(dia.holdout.data.recalculated.log2 ~ targeted.holdout.data.log2);

  modelsummaries<-c(summary(linearmodel.prefit),summary(linearmodel.postfit.all),summary(linearmodel.postfit.training),summary(linearmodel.postfit.holdout));

  #dev.new();
  dev.current<-dev.cur();
  par(mfrow=c(3,2));
  ad.cols<-vector('character',total_samples);
  ad.cols[1:total_samples]<-'red';
  ad.cols[sampling_indices]<-'black';
  plot(dia.data.prefit.log2 ~ targeted.data.log2,main='Pre-pick Raw Fit');
  plot(probe.GA,main='Genetic Algorithm Results')
  plot(dia.data.modeled.recalculated.log2 ~ targeted.data.log2,main='Post-pick All Data Fit',col=ad.cols)
  plot(dia.training.data.recalculated.log2 ~ targeted.training.data.log2,main='Post-pick Training Data Fit')
  plot(dia.holdout.data.recalculated.log2 ~ targeted.holdout.data.log2,main='Post-pick Holdout Data Fit')

  return.models<-list(prefit=linearmodel.prefit,postfit.all=linearmodel.postfit.all,postfit.training=linearmodel.postfit.training,postfit.holdout=linearmodel.postfit.holdout)
  return.object<-list(probe=probe_name_code,transitions=probe.minimal.transitions,models=return.models,GA=probe.GA);
}

fix_skyline_csv <- function (csvfile) {
  newcolumn<-paste(csvfile[,3],csvfile[,4],sep="_");
  csvfile[,3]<-newcolumn;
  csv.dims<-dim(csvfile);
  newcsv<-csvfile[,c(1,3,5:csv.dims[2])];
  return(newcsv);
}
