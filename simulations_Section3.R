library(TreeSim)
library(mvSLOUCH)
library(RMaCzek)
library(grDevices)
library(randtests)
library(seriation)


f_cluster_assess<-function(ordered_obs,obs_clusters,obs_labels,b_onlydiag=TRUE){
## The function checks how well a seriation corresponds to a provided true clustering.
## The observations are "arranged" on a straight line and consecutive observations are either in the same
## cluster or on the boundaries of clusters, i.e. 1,1,1,1|2,2|3,3,3,3 where the number i means that
## the observation is in cluster i.
## Input parameters
## ordered_obs: a permutation of the observations.
## obs_clusters: a list with each field corresponding to a cluster, i.e. the observations inside each cluster
## obs_labels: if the values inside obs_clusters are not the observations themselves but their indices. Then,
##	obs_labels are the values of the observations that correspond to the values in ordered_obs.
## b_onlydiag: logical should the function only return a vector indicating how many observations were correctly
##	clustered together in each cluster (default TRUE) or a matrix as described in returned value (FALSE)
## Value
## either a vector showing for each cluster how many observations are correctly clustered together
## 	or a matrix showing for each cluster how many observations from this cluster were put into every cluster
## Note that the order of the clusters in the output does not need to correspond to the order of the provided clusters.
## The function finds the permutation of the numbers 1:number of clusters that will maximize the success of the seriation.
    v_final_order<-rep(1,length(ordered_obs))
    for (i in 2:length(obs_clusters)){
	v_final_order[which(is.element(ordered_obs,obs_labels[obs_clusters[[i]]]))]<-i
    }
    m_cluster_Assessment<-matrix(0,length(obs_clusters),length(obs_clusters))
    m_cluster_Assessment[1,]<-table(c(1:length(obs_clusters),v_final_order[1:length(obs_clusters[[1]])]))-1
    currlen<-length(obs_clusters[[1]])
    if (length(obs_clusters)>1){
	for (i in 2:length(obs_clusters)){
	    m_cluster_Assessment[i,]<-table(c(1:length(obs_clusters),v_final_order[(currlen+1):(currlen+length(obs_clusters[[i]]))]))-1
	    currlen<-currlen+length(obs_clusters[[i]])
	}
    }
    all_perms<-permut(1:length(obs_clusters))
    cluster_scores<-apply(all_perms,1,function(x,m_cluster_Assessment){
	sum((diag(m_cluster_Assessment[x,])/apply(m_cluster_Assessment,2,sum)-1)^2)
    },m_cluster_Assessment=m_cluster_Assessment)
    index_best_permut<-which(cluster_scores == min(cluster_scores), arr.ind=TRUE)
    if (length(index_best_permut)>1){
	## we have more than one best assignement, we make a random choice right now
	index_best_permut<-sample(index_best_permut,1)
    }
    m_cluster_Assessment<-m_cluster_Assessment[all_perms[index_best_permut,],]
    res<-NA
    if(b_onlydiag){res<-diag(m_cluster_Assessment)}
    else{res<-m_cluster_Assessment}
    res
}


f_doRMaCzekAnalysis<-function(b_doplot,b_dosave,runnum,random_shuffle_of_data,num_traits,mPsi,num_repeats=100,graphic_pars=list(xmargins=8,xoffset=4.1,ymargins=7,yoffset=4.9,pointsize=0.55,spacesize=0.55),v_clustsize=c(30,30,30)){
## The function does a simulation study of the RMaCzek package. Each rerun simulates a clustered phylogeny, then, using mvSLOUCH, simulates 
## a multidimensional Ornstein-Uhlenbeck (OU) process of top of it. Apart from the optimum, the parameters of the OU process
## are created randomly. Afterwords, it uses various seriation methods from the RMaCzek to see
## how well it can recover the clades on the phylogeny.
## Input parameters
## b_doplot: logical, if TRUE, then Czekanowski's diagrams of the  last rerun are also plotted
## b_dosave: logical, if TRUE, then the results of the simulation study are saved in an .RData file
## runnum: string defining the study, used to create the file name of the .RData file if b_dosave is TRUE
## random_shuffle_of_data: logical, should the simulated OU data rows be shuffled (TRUE) or provided in the correct order
##	(FALSE, observtions from the same clade in adjacent rows) to the RMaCzke::czek_matrix() function.
## num_traits: the dimension of the simulated OU process
## mPsi: the optimum value inside each clade of the phylogeny.
## num_repeats: how many times should the simulation be rerun (default 100)
## graphic_pars: a list of parameters that are to aid in creating Czekanowski's diagram. The list should contain the following fields
## 		xmargins (default 8): size of x margins
## 		xoffset (default 4.1): offset on the x axis 
##		ymargins (default 7): size of y margins
## 		yoffset (default 4.9): offset on the y axis 
##		pointsize (default 0.55): size of coloured point on the x-axis labels
##		spacesize (default 0.55): size of space between points on the x-axis labels
##	These parameters are passed to the RMaCzek:::plot.czek_matrix() function.
## v_clustsize: a vector of cluster/clade sizes (defualt c(30,30,30))
## Value
## a list of length equalling num_repeats. Each field is a rerun on the simulation and contains the fields:
##	phyltree: the simulated clustered phylogeny
##	regimes: regimes (branch indices) with different levels of the optimum vector
##	OUOUparams: parameters of the OU process
##	OUOUdata: simulated OU data (shuffled if random_shuffle_of_data is TRUE)
##	ga_res: result of RMaCzek::czek_matrix() with the ga seriation method (genetic algorith included in RMaCzek)
##	ga_Um: value of Um factor from the ga seriation method (genetic algorith included in RMaCzek)
##	ga_path_length: value of path length factor from the ga seriation method (genetic algorith included in RMaCzek)
##	ga_clusts: number of correct observations in each cluster under the ga seriation method (genetic algorith included in RMaCzek)
##	ga_time: running time under the ga seriation method (genetic algorith included in RMaCzek)
##	QAP_2SUM_res: result of RMaCzek::czek_matrix() with the QAP_2SUM seriation method (from package seriation)
##	QAP_2SUM_Um: value of Um factor from the QAP_2SUM seriation method (from package seriation)
##	QAP_2SUM_path_length: value of path length factor from the QAP_2SUM seriation method (from package seriation)
##	QAP_2SUM_clusts: number of correct observations in each cluster under the QAP_2SUM seriation method (from package seriation)
##	QAP_2SUM_time: running time under the QAP_2SUM seriation method (from package seriation)
##	OLO_res: result of RMaCzek::czek_matrix() with the OLO seriation method (from package seriation)
##	OLO_Um: value of Um factor from the OLO seriation method (from package seriation)
##	OLO_path_length: value of path length factor from the OLO seriation method (from package seriation)
##	OLO_clusts: number of correct observations in each cluster under the OLO seriation method (from package seriation)
##	OLO_time: running time under the OLO seriation method (from package seriation)
    RMaCzek_res<-vector("list",num_repeats)    
    numobs<-sum(v_clustsize)
    v_plot_methods<-c("ga","OLO","QAP_2SUM")
    for (j in 1:num_repeats){
	print(j)
	RMaCzek_res[[j]]<-list(phyltree=NULL, regimes=NULL, OUOUparams=NULL, OUOUdata=NULL, ga_Um=NULL, ga_res=NULL, ga_clusts=NULL, OLO_Um=NULL, OLO_res=NULL, OLO_clusts=NULL, QAP_2SUM_Um=NULL, QAP_2SUM_res=NULL, QAP_2SUM_clusts=NULL)

	phylophyltree<-mvSLOUCH::simulate_clustered_phylogeny(v_clustsize,joining_branchlengths=c(20,NA),f_simclustphyl="sim.bd.taxa_Yule1",joiningphyl="sim.bd.taxa_Yule1",b_change_joining_branches=TRUE)
	v_regimes<-rep("reg1",sum(sapply(phylophyltree$edges_clusters,function(x){length(x)},simplify=TRUE)))
	v_regimes[phylophyltree$edges_clusters[[3]]]<-"reg2"
	v_regimes[phylophyltree$edges_clusters[[4]]]<-"reg3"
	Mtmp<-matrix(rnorm(num_traits*num_traits),num_traits,num_traits)
	Mtmp<-Mtmp%*%t(Mtmp)
	mSyy<-t(chol(Mtmp))
	OUOUparameters<-list(vY0=matrix(runif(num_traits),nrow=num_traits,ncol=1),
	A=diag(rexp(num_traits)),mPsi=mPsi,Syy=mSyy)
	OUOUdata<-simulOUCHProcPhylTree(phylophyltree,OUOUparameters,regimes=v_regimes)

	if(random_shuffle_of_data){
	    OUOUdata<-OUOUdata[sample(1:nrow(OUOUdata),nrow(OUOUdata)),,drop=FALSE]
	}

	RMaCzek_res[[j]]$phyltree<-phylophyltree
	RMaCzek_res[[j]]$OUOUparams<-OUOUparameters
	RMaCzek_res[[j]]$OUOUdata<-OUOUdata
	RMaCzek_res[[j]]$regimes<-v_regimes

	s_anal_time<-Sys.time()
	CzkMat_QAP_2SUM<-czek_matrix(OUOUdata,n_classes = 5,order="QAP_2SUM")
	e_anal_time<-Sys.time()
	RMaCzek_res[[j]]$QAP_2SUM_res<-CzkMat_QAP_2SUM
	RMaCzek_res[[j]]$QAP_2SUM_clusts<-f_cluster_assess(rownames(CzkMat_QAP_2SUM)[attr(CzkMat_QAP_2SUM,"order")],phylophyltree$tips_clusters,phylophyltree$tip.label,TRUE)
	RMaCzek_res[[j]]$QAP_2SUM_Um<-RMaCzek::Um_factor(dist(scale(OUOUdata)),order=attr(CzkMat_QAP_2SUM,"order"),inverse_um=FALSE)
	RMaCzek_res[[j]]$QAP_2SUM_path_length<-seriation::criterion(dist(scale(OUOUdata)),order=seriation::ser_permutation(attr(CzkMat_QAP_2SUM,"order")),method="Path_length")
	RMaCzek_res[[j]]$QAP_2SUM_time<-e_anal_time-s_anal_time
	print("done QPA_2SUM in time")
	print(RMaCzek_res[[j]]$QAP_2SUM_time)
	

	s_anal_time<-Sys.time()
	CzkMat_OLO<-czek_matrix(OUOUdata,n_classes = 5,order="OLO")
	e_anal_time<-Sys.time()
	RMaCzek_res[[j]]$OLO_res<-CzkMat_OLO
	RMaCzek_res[[j]]$OLO_clusts<-f_cluster_assess(rownames(CzkMat_QAP_2SUM)[attr(CzkMat_OLO,"order")],phylophyltree$tips_clusters,phylophyltree$tip.label,TRUE)
	RMaCzek_res[[j]]$OLO_Um<-RMaCzek::Um_factor(dist(scale(OUOUdata)),order=attr(CzkMat_OLO,"order"),inverse_um=FALSE)
	RMaCzek_res[[j]]$OLO_path_length<-seriation::criterion(dist(scale(OUOUdata)),order=seriation::ser_permutation(attr(CzkMat_OLO,"order")),method="Path_length")
	RMaCzek_res[[j]]$OLO_time<-e_anal_time-s_anal_time
	print("done OLO in time")
	print(RMaCzek_res[[j]]$OLO_time)

	s_anal_time<-Sys.time()
	CzkMat_ga<-czek_matrix(OUOUdata,n_classes = 5,order="ga")
	e_anal_time<-Sys.time()
	RMaCzek_res[[j]]$ga_res<-CzkMat_ga
	RMaCzek_res[[j]]$ga_clusts<-f_cluster_assess(rownames(CzkMat_ga)[attr(CzkMat_ga,"order")],phylophyltree$tips_clusters,phylophyltree$tip.label,TRUE)
	RMaCzek_res[[j]]$ga_Um<-RMaCzek::Um_factor(dist(scale(OUOUdata)),order=attr(CzkMat_ga,"order"),inverse_um=FALSE)
	RMaCzek_res[[j]]$ga_path_length<-seriation::criterion(dist(scale(OUOUdata)),order=seriation::ser_permutation(attr(CzkMat_ga,"order")),method="Path_length")
	RMaCzek_res[[j]]$ga_time<-e_anal_time-s_anal_time
	print("done GA in time")
	print(RMaCzek_res[[j]]$ga_time)
	##units(res)<-"weeks"
    }

    if (b_dosave){
	filename<-"RMaCzek_"
	if (random_shuffle_of_data){
	    filename<-paste0(filename,"shuffledorder")
	}else{filename<-paste0(filename,"correctorder")}
	filename<-paste0(filename,"_test_",runnum,".RData")
	save(RMaCzek_res,file=filename)
    }
    if (b_doplot){
	xwidth<-(graphic_pars$spacesize*(numobs-1)+graphic_pars$xmargins+numobs*graphic_pars$pointsize)
	ywidth<-(graphic_pars$spacesize*(numobs-1)+graphic_pars$ymargins+numobs*graphic_pars$pointsize)

	for (meth in v_plot_methods){
	    CzMat_res<-switch(meth,
		ga=CzkMat_ga,
		OLO=CzkMat_OLO,
		QAP_2SUM=CzkMat_QAP_2SUM
	    )
	    v_ordered_obs<-rownames(CzMat_res)[attr(CzMat_res,"order")]
	    axiscolors <- rep(gray(0),numobs)
	    axiscolors[which(is.element(v_ordered_obs,phylophyltree$tip.label[phylophyltree$tips_clusters[[2]]]))]<-gray(0.5)
	    axiscolors[which(is.element(v_ordered_obs,phylophyltree$tip.label[phylophyltree$tips_clusters[[3]]]))]<-gray(0.75)
	    pdf(paste0("phylRMaCzek_",meth,".pdf"))
	    plot(CzMat_res,axis=TRUE,plot_title=meth)
	    points(graphic_pars$xoffset/xwidth+((xwidth-graphic_pars$xmargins)/xwidth)*((1:numobs)-1)/(numobs),rep(0.01,numobs),col=axiscolors,pch=19,cex=0.4)
    	    points(rep(0.01,numobs),graphic_pars$yoffset/ywidth+((ywidth-graphic_pars$ymargins)/ywidth)*((1:numobs)-1)/(numobs),col=rev(axiscolors),pch=19,cex=0.4)
	    dev.off()
	}
	
    }
    RMaCzek_res
}

fprres<-function(filename,N=30,correct_Um=FALSE,calc_path_len=FALSE,btoprint=c(clusts=TRUE,Um=TRUE,path_len=TRUE,time=TRUE)){
## Function prints tables that summarize (report mean and variance) the cluster correctness Um, path_length criteria, running times for the 
## ga, QAP_2SUM and OLO seriation methods from the output of the f_doRMaCzekAnalysis() function
## Input parameters
## filename: name of file from which to read in the object
## N: size of each cluster (default 30), assumes that all the clusters have the same size
## correct_Um: logical, should the Um factor be recalculated (default FALSE, do not)
## calc_path_len: logical, should the path length factor be recalculated (default FALSE, do not)
## btoprint: logical, named, vector which tables to print out, by default all TRUE, entries:
##	clusts: should cluster correctness be summarized
##	Um: should Um factors be summarized
##	path_len: should path length factors be summarized
##	time: should running times be summarized
## Value
##	Nothing
    load(filename)

    print(filename)
    if (btoprint["clusts"]){
	print("clustering")
	M<-rbind(
	    ga_mean=apply(sapply(RMaCzek_res,function(x){x$ga_clusts}),1,mean),
	    ga_sd=apply(sapply(RMaCzek_res,function(x){x$ga_clusts}),1,sd),
	    QAP_2SUM_mean=apply(sapply(RMaCzek_res,function(x){x$QAP_2SUM_clusts}),1,mean),
	    QAP_2SUM_sd=apply(sapply(RMaCzek_res,function(x){x$QAP_2SUM_clusts}),1,sd),
    	    OLO_mean=apply(sapply(RMaCzek_res,function(x){x$OLO_clusts}),1,mean),
	    OLO_sd=apply(sapply(RMaCzek_res,function(x){x$OLO_clusts}),1,sd)
	)/N
	colnames(M)<-c("cluster_1","cluster_2","cluster_3")
	print(M)
    }


    if (btoprint["Um"]){
	if (correct_Um){
	    RMaCzek_res<-sapply(RMaCzek_res,function(x){
		x$ga_Um<-RMaCzek::Um_factor(dist(scale(x$OUOUdata)),order=attr(x$ga_res,"order"),inverse_um=FALSE)
		x$QAP_2SUM_Um<-RMaCzek::Um_factor(dist(scale(x$OUOUdata)),order=attr(x$QAP_2SUM_res,"order"),inverse_um=FALSE)
		x$OLO_Um<-RMaCzek::Um_factor(dist(scale(x$OUOUdata)),order=attr(x$OLO_res,"order"),inverse_um=FALSE)		
		x
	    },simplify=FALSE)
	    save(RMaCzek_res,file=filename)
	}
	print("Um")
	print(rbind(
	    ga_mean=mean(sapply(RMaCzek_res,function(x){x$ga_Um})),
	    ga_sd=sd(sapply(RMaCzek_res,function(x){x$ga_Um})),
	    QAP_2SUM_mean=mean(sapply(RMaCzek_res,function(x){x$QAP_2SUM_Um})),
	    QAP_2SUM_sd=sd(sapply(RMaCzek_res,function(x){x$QAP_2SUM_Um})),
	    OLO_mean=mean(sapply(RMaCzek_res,function(x){x$OLO_Um})),
	    OLO_sd=sd(sapply(RMaCzek_res,function(x){x$OLO_Um}))
	))
    }

    if (btoprint["path_len"]){
	if (calc_path_len){
	    RMaCzek_res<-sapply(RMaCzek_res,function(x){
		x$ga_path_length<-seriation::criterion(dist(scale(x$OUOUdata)),order=seriation::ser_permutation(attr(x$ga_res,"order")),method="Path_length")
		x$QAP_2SUM_path_length<-seriation::criterion(dist(scale(x$OUOUdata)),order=seriation::ser_permutation(attr(x$QAP_2SUM_res,"order")),method="Path_length")
		x$OLO_path_length<-seriation::criterion(dist(scale(x$OUOUdata)),order=seriation::ser_permutation(attr(x$OLO_res,"order")),method="Path_length")
		x
	    },simplify=FALSE)
	    save(RMaCzek_res,file=filename)
	}
	print("OLO path length")
	print(rbind(
	    ga_mean=mean(sapply(RMaCzek_res,function(x){x$ga_path_length})),
	    ga_sd=sd(sapply(RMaCzek_res,function(x){x$ga_path_length})),
	    QAP_2SUM_mean=mean(sapply(RMaCzek_res,function(x){x$QAP_2SUM_path_length})),
	    QAP_2SUM_sd=sd(sapply(RMaCzek_res,function(x){x$QAP_2SUM_path_length})),
	    OLO_mean=mean(sapply(RMaCzek_res,function(x){x$OLO_path_length})),
	    OLO_sd=sd(sapply(RMaCzek_res,function(x){x$OLO_path_length}))
	))
    }

    if (btoprint["time"]){
	print("timing")
	print(rbind(
	    ga_mean=mean(sapply(RMaCzek_res,function(x){as.numeric(x$ga_time,units="secs")})),
	    ga_sd=sd(sapply(RMaCzek_res,function(x){as.numeric(x$ga_time,units="secs")})),
	    QAP_2SUM_mean=mean(sapply(RMaCzek_res,function(x){as.numeric(x$QAP_2SUM_time,units="secs")})),
	    QAP_2SUM_mean=sd(sapply(RMaCzek_res,function(x){as.numeric(x$QAP_2SUM_time,units="secs")})),
	    OLO_mean=mean(sapply(RMaCzek_res,function(x){as.numeric(x$OLO_time,units="secs")})),
	    OLO_sd=sd(sapply(RMaCzek_res,function(x){as.numeric(x$OLO_time,units="secs")}))
	))
    }
    print("================================================")
    invisible(NA)
}

numreps<-100

num_traits<-10 ## number of traits
runnum<-"1" ## index of analysis

v_clustsize<-c(30,30,30)
## defining the random number generator
RNGversion("3.6.1") 
rexp(1)
load(file=paste0("Random_seed_",runnum,".RData"))

## do six different analyses 
## shuffled/correct order of data
## optima for each cluster are same/similar/very different
f_doRMaCzekAnalysis(FALSE,TRUE,paste0(runnum,"_distinctmeans"),TRUE,num_traits,mPsi=cbind("reg1"=rnorm(num_traits,mean=-50),"reg2"=rnorm(num_traits,mean=0),"reg3"=rnorm(num_traits,mean=50)),num_repeats=numreps,v_clustsize=v_clustsize)
f_doRMaCzekAnalysis(FALSE,TRUE,paste0(runnum,"_distinctmeans"),FALSE,num_traits,mPsi=cbind("reg1"=rnorm(num_traits,mean=-50),"reg2"=rnorm(num_traits,mean=0),"reg3"=rnorm(num_traits,mean=50)),num_repeats=numreps,v_clustsize=v_clustsize)

f_doRMaCzekAnalysis(FALSE,TRUE,paste0(runnum,"_similarmeans"),TRUE,num_traits,mPsi=cbind("reg1"=rnorm(num_traits,mean=-2),"reg2"=rnorm(num_traits,mean=0),"reg3"=rnorm(num_traits,mean=2)),num_repeats=numreps,v_clustsize=v_clustsize)
f_doRMaCzekAnalysis(FALSE,TRUE,paste0(runnum,"_similarmeans"),FALSE,num_traits,mPsi=cbind("reg1"=rnorm(num_traits,mean=-2),"reg2"=rnorm(num_traits,mean=0),"reg3"=rnorm(num_traits,mean=2)),num_repeats=numreps,v_clustsize=v_clustsize)

f_doRMaCzekAnalysis(FALSE,TRUE,paste0(runnum,"_samemeans"),TRUE,num_traits,mPsi=cbind("reg1"=rnorm(num_traits,mean=0),"reg2"=rnorm(num_traits,mean=0),"reg3"=rnorm(num_traits,mean=0)),num_repeats=numreps,v_clustsize=v_clustsize)
f_doRMaCzekAnalysis(FALSE,TRUE,paste0(runnum,"_samemeans"),FALSE,num_traits,mPsi=cbind("reg1"=rnorm(num_traits,mean=0),"reg2"=rnorm(num_traits,mean=0),"reg3"=rnorm(num_traits,mean=0)),num_repeats=numreps,v_clustsize=v_clustsize)

## summarize the cluster quality
options(scipen=999)
sink("ClusterQuality.txt")
fprres(paste0("RMaCzek_correctorder_test_",runnum,"_distinctmeans.RData"),correct_Um=FALSE,calc_path_len=FALSE,btoprint=c(clusts=TRUE,Um=TRUE,path_len=TRUE,time=TRUE))
fprres(paste0("RMaCzek_shuffledorder_test_",runnum,"_distinctmeans.RData"),correct_Um=FALSE,calc_path_len=FALSE,btoprint=c(clusts=TRUE,Um=TRUE,path_len=TRUE,time=TRUE))
fprres(paste0("RMaCzek_shuffledorder_test_",runnum,"_similarmeans.RData"),correct_Um=FALSE,calc_path_len=FALSE,btoprint=c(clusts=TRUE,Um=TRUE,path_len=TRUE,time=TRUE))
fprres(paste0("RMaCzek_correctorder_test_",runnum,"_similarmeans.RData"),correct_Um=FALSE,calc_path_len=FALSE,btoprint=c(clusts=TRUE,Um=TRUE,path_len=TRUE,time=TRUE))
fprres(paste0("RMaCzek_correctorder_test_",runnum,"_samemeans.RData"),correct_Um=FALSE,calc_path_len=FALSE,btoprint=c(clusts=TRUE,Um=TRUE,path_len=TRUE,time=TRUE))
fprres(paste0("RMaCzek_shuffledorder_test_",runnum,"_samemeans.RData"),correct_Um=FALSE,calc_path_len=FALSE,btoprint=c(clusts=TRUE,Um=TRUE,path_len=TRUE,time=TRUE))
sink()

## plot an example clustered phylogeny
## The exact plot here will depend on the random seed
load(file="Random_seed_plot.RData")
output<-f_doRMaCzekAnalysis(TRUE,FALSE,NULL,FALSE,num_traits,mPsi=cbind("reg1"=rnorm(num_traits,mean=-2),"reg2"=rnorm(num_traits,mean=0),"reg3"=rnorm(num_traits,mean=2)),num_repeats=1)
pdf("phyltree.pdf");plot(output[[1]]$phyltree,clust_cols=c(gray(0),gray(0.75),gray(0.5)),clust_edge.width=3,joiningphylo_edge.width=3,show.tip.label=FALSE);dev.off()






