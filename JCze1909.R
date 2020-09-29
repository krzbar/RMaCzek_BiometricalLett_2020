## Analysis of Czekanowski (1909)'s skulls data (Table II), found in RMaCzek as skulls_distances.
## The produced results (ordering of observations) can depend on the random seed.
## Here, we load the saved random seed that resulted in the best orderings.


library(RMaCzek)

RNGversion("3.6.1") 
rexp(1)
load("JCze1909_randomseed.RData")

## correcting typo rowname
rownames(skulls_distances)[3]<-"Krapina C"

## we need to symmetrize the distance matrix as in the original
## paper, Czekanowski (1909), there is a minor typographical error
## in the provided distance matrix, based on Stolyhwa (1908)'s original data
sym_skulls_distances<-skulls_distances
sym_skulls_distances[5,9]<-10.504

czek_matrix_skulls<-czek_matrix(as.dist(sym_skulls_distances),order=NA,original_diagram=TRUE)
czek_matrix_skulls_sym<-czek_matrix(as.dist(sym_skulls_distances),order=NA,original_diagram=FALSE)
czek_matrix_skulls_OLO<-czek_matrix(as.dist(sym_skulls_distances),order="OLO",original_diagram=TRUE)
czek_matrix_skulls_OLO_sym<-czek_matrix(as.dist(sym_skulls_distances),order="OLO",original_diagram=FALSE)
czek_matrix_skulls_qap2sum_sym<-czek_matrix(as.dist(sym_skulls_distances),order="QAP_2SUM",original_diagram=FALSE)
czek_matrix_skulls_qap2sum<-czek_matrix(as.dist(sym_skulls_distances),order="QAP_2SUM",original_diagram=TRUE)


pdf("Czekanowski_skulls_orig.pdf");plot(czek_matrix_skulls,plot_title="",label.cex=0.5);dev.off()
pdf("Czekanowski_skulls_sym.pdf");plot(czek_matrix_skulls_sym,plot_title="",label.cex=0.5);dev.off()
pdf("Czekanowski_skulls_OLO.pdf");plot(czek_matrix_skulls_OLO,plot_title="",label.cex=0.5);dev.off()
pdf("Czekanowski_skulls_OLO_sym.pdf");plot(czek_matrix_skulls_OLO_sym,plot_title="",label.cex=0.5);dev.off()
pdf("Czekanowski_skulls_qap2sum_sym.pdf");plot(czek_matrix_skulls_qap2sum_sym,plot_title="",label.cex=0.5);dev.off()
pdf("Czekanowski_skulls_qap2sum.pdf");plot(czek_matrix_skulls_qap2sum,plot_title="",label.cex=0.5);dev.off()


## if it turned out that the Um values are not the same
if (attr(czek_matrix_skulls_qap2sum_sym,"Um")<attr(czek_matrix_skulls_qap2sum,"Um")){
    print("Symmetric diagram had better Um factor.")
    czek_matrix_skulls_qap2sum<-czek_matrix(as.dist(sym_skulls_distances),order=attr(czek_matrix_skulls_qap2sum_sym,"order"),original_diagram=TRUE)
    attr(czek_matrix_skulls_qap2sum,"criterion_value")<-attr(czek_matrix_skulls_qap2sum_sym,"criterion_value")
    pdf("Czekanowski_skulls_qap2sum.pdf");plot(czek_matrix_skulls_qap2sum,plot_title="",label.cex=0.5);dev.off()
}else{
    if (attr(czek_matrix_skulls_qap2sum,"Um")<attr(czek_matrix_skulls_qap2sum_sym,"Um")){
	print("Asymmetric diagram had better Um factor.")
        czek_matrix_skulls_qap2sum_sym<-czek_matrix(as.dist(sym_skulls_distances),order=attr(czek_matrix_skulls_qap2sum,"order"),original_diagram=FALSE)
        attr(czek_matrix_skulls_qap2sum_sym,"criterion_value")<-attr(czek_matrix_skulls_qap2sum,"criterion_value")
	pdf("Czekanowski_skulls_qap2sum_sym.pdf");plot(czek_matrix_skulls_qap2sum_sym,plot_title="",label.cex=0.5);dev.off()
    }
}
##=================================================================


## if it turned out that the Path length values are not the same
if (attr(czek_matrix_skulls_OLO_sym,"Path_length")<attr(czek_matrix_skulls_OLO,"Path_length")){
    print("Symmetric diagram had better OLO factor.")
    czek_matrix_skulls_OLO<-czek_matrix(as.dist(sym_skulls_distances),order=attr(czek_matrix_skulls_OLO_sym,"order"),original_diagram=TRUE)
    attr(czek_matrix_skulls_OLO,"criterion_value")<-attr(czek_matrix_skulls_OLO_sym,"criterion_value")
    pdf("Czekanowski_skulls_OLO.pdf");plot(czek_matrix_skulls_OLO,plot_title="",label.cex=0.5);dev.off()
}else{
    if (attr(czek_matrix_skulls_OLO,"Path_length")<attr(czek_matrix_skulls_OLO_sym,"Path_length")){
	print("Asymmetric diagram had better OLO factor.")    
        czek_matrix_skulls_OLO_sym<-czek_matrix(as.dist(sym_skulls_distances),order=attr(czek_matrix_skulls_OLO,"order"),original_diagram=FALSE)
        attr(czek_matrix_skulls_OLO_sym,"criterion_value")<-attr(czek_matrix_skulls_OLO,"criterion_value")
	pdf("Czekanowski_skulls_OLO_sym.pdf");plot(czek_matrix_skulls_OLO_sym,plot_title="",label.cex=0.5);dev.off()
    }
}
##=================================================================