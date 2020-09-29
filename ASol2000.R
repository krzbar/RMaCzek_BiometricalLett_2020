## Analysis of Soltysiak (2000)'s seals data (Table 10), found in RMaCzek as seals_similarities.
## The produced results (ordering of observations) can depend on the random seed.
## Here, we load the saved random seed that resulted in the best orderings.

library(RMaCzek)
RNGversion("3.6.1") 
rexp(1)
load("ASol2000_randomseed.RData")


diag(seals_similarities)<-100
dist_soltysiak_seals<-as.dist(100-seals_similarities)
int_breaks<-c(0,40,60,80,100)

czkm_seals_orig<-czek_matrix(dist_soltysiak_seals,original_diagram=TRUE,order=NA)
czkm_seals_orig_sym<-czek_matrix(dist_soltysiak_seals,original_diagram=FALSE,order=NA,n_classes=length(int_breaks)-1,interval_breaks=int_breaks)
czkm_seals_OLO<-czek_matrix(dist_soltysiak_seals,original_diagram=TRUE,order="OLO")
czkm_seals_OLO_sym<-czek_matrix(dist_soltysiak_seals,original_diagram=FALSE,order="OLO",n_classes=length(int_breaks)-1,interval_breaks=int_breaks)
czkm_seals_qap2sum<-czek_matrix(dist_soltysiak_seals,original_diagram=TRUE,order="QAP_2SUM")
czkm_seals_qap2sum_sym<-czek_matrix(dist_soltysiak_seals,original_diagram=FALSE,order="QAP_2SUM",n_classes=length(int_breaks)-1,interval_breaks=int_breaks)


pdf("Soltysiak_seals_orig.pdf");plot(czkm_seals_orig,plot_title="",label.cex=0.7);dev.off()
pdf("Soltysiak_seals_OLO.pdf");plot(czkm_seals_OLO,plot_title="",label.cex=0.7);dev.off()
pdf("Soltysiak_seals_qap2sum.pdf");plot(czkm_seals_qap2sum,plot_title="",label.cex=0.7);dev.off()
pdf("Soltysiak_seals_orig_sym.pdf");plot(czkm_seals_orig_sym,plot_title="",label.cex=0.7);dev.off()
pdf("Soltysiak_seals_OLO_sym.pdf");plot(czkm_seals_OLO_sym,plot_title="",label.cex=0.7);dev.off()
pdf("Soltysiak_seals_qap2sum_sym.pdf");plot(czkm_seals_qap2sum_sym,plot_title="",label.cex=0.7);dev.off()


## if it turned out that the Um values are not the same
if (attr(czkm_seals_qap2sum_sym,"Um")<attr(czkm_seals_qap2sum,"Um")){
    print("Symmetric diagram had better Um factor.")
    czkm_seals_qap2sum<-czek_matrix(dist_soltysiak_seals,original_diagram=TRUE,order=attr(czkm_seals_qap2sum_sym,"order"))
    attr(czkm_seals_qap2sum,"criterion_value")<-attr(czkm_seals_qap2sum_sym,"criterion_value")
    pdf("Soltysiak_seals_qap2sum.pdf");plot(czkm_seals_qap2sum,plot_title="",label.cex=0.7);dev.off()
}else{
    if (attr(czkm_seals_qap2sum,"Um")<attr(czkm_seals_qap2sum_sym,"Um")){
        print("Asymmetric diagram had better Um factor.")
	czkm_seals_qap2sum_sym<-czek_matrix(dist_soltysiak_seals,original_diagram=FALSE,order=attr(czkm_seals_qap2sum,"order"))
	attr(czkm_seals_qap2sum_sym,"criterion_value")<-attr(czkm_seals_qap2sum,"criterion_value")
	pdf("Soltysiak_seals_qap2sum_sym.pdf");plot(czkm_seals_qap2sum_sym,plot_title="",label.cex=0.7);dev.off()
    }
}
##=================================================================


## if it turned out that the Path length values are not the same
if (attr(czkm_seals_OLO_sym,"Path_length")<attr(czkm_seals_OLO,"Path_length")){
     print("Symmetric diagram had better OLO factor.")
    czkm_seals_OLO<-czek_matrix(dist_soltysiak_seals,original_diagram=TRUE,order=attr(czkm_seals_OLO_sym,"order"))
    attr(czkm_seals_OLO,"criterion_value")<-attr(czkm_seals_OLO_sym,"criterion_value")
    pdf("Soltysiak_seals_OLO.pdf");plot(czkm_seals_OLO,plot_title="",label.cex=0.7);dev.off()
}else{
    if (attr(czkm_seals_OLO,"Path_length")<attr(czkm_seals_OLO_sym,"Path_length")){
	print("Asymmetric diagram had better OLO factor.")
	czkm_seals_OLO_sym<-czek_matrix(dist_soltysiak_seals,original_diagram=FALSE,order=attr(czkm_seals_OLO,"order"))
	attr(czkm_seals_OLO_sym,"criterion_value")<-attr(czkm_seals_OLO,"criterion_value")
	pdf("Soltysiak_seals_OLO_sym.pdf");plot(czkm_seals_OLO_sym,plot_title="",label.cex=0.7);dev.off()
    }
}
##=================================================================

## compare the two OLO arrangements
print("Ordering and factor values under original OLO diagram")
print(czkm_seals_OLO)
print("Ordering and factor values under symmetric OLO diagram")
print(czkm_seals_OLO_sym)
