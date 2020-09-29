## Analysis of Soltysiak, Jaskulski (1999)'s urns data (Table p. 181), found in RMaCzek as urns.
## The produced results (ordering of observations) can depend on the random seed.
## Here, we load the saved random seed that resulted in the best orderings.


library(RMaCzek)
RNGversion("3.6.1") 
rexp(1)
load("ASolPJas1999_randomseed.RData")


soltysiakjaskulski_urns_use<-urns[-9,] ## removed as too many missing values, observation "gr.52.1"

soltysiakjaskulski_urns_use_proposedorder<-c(1,3,5,7,2,4,8,11,10,12,6,9,14,13)
czkm_urns_orig<-czek_matrix(soltysiakjaskulski_urns_use,original_diagram=TRUE,order=soltysiakjaskulski_urns_use_proposedorder,scale_data=TRUE)
czkm_urns_orig_sym<-czek_matrix(soltysiakjaskulski_urns_use,original_diagram=FALSE,order=soltysiakjaskulski_urns_use_proposedorder,scale_data=TRUE)
czkm_urns_OLO<-czek_matrix(soltysiakjaskulski_urns_use,original_diagram=TRUE,order="OLO",scale_data=TRUE)
czkm_urns_OLO_sym<-czek_matrix(soltysiakjaskulski_urns_use,original_diagram=FALSE,order="OLO",scale_data=TRUE)
czkm_urns_qap2sum<-czek_matrix(soltysiakjaskulski_urns_use,original_diagram=TRUE,order="QAP_2SUM",scale_data=TRUE)
czkm_urns_qap2sum_sym<-czek_matrix(soltysiakjaskulski_urns_use,original_diagram=FALSE,order="QAP_2SUM",scale_data=TRUE)


pdf("SoltysiakJaskulski_urns_orig.pdf");plot(czkm_urns_orig,plot_title="",label.cex=0.9);dev.off()
pdf("SoltysiakJaskulski_urns_OLO.pdf");plot(czkm_urns_OLO,plot_title="",label.cex=0.9);dev.off()
pdf("SoltysiakJaskulski_urns_qap2sum.pdf");plot(czkm_urns_qap2sum,plot_title="",label.cex=0.9);dev.off()
pdf("SoltysiakJaskulski_urns_orig_sym.pdf");plot(czkm_urns_orig_sym,plot_title="",label.cex=0.9);dev.off()
pdf("SoltysiakJaskulski_urns_OLO_sym.pdf");plot(czkm_urns_OLO_sym,plot_title="",label.cex=0.9);dev.off()
pdf("SoltysiakJaskulski_urns_qap2sum_sym.pdf");plot(czkm_urns_qap2sum_sym,plot_title="",label.cex=0.9);dev.off()

## if it turned out that the Um values are not the same
if (attr(czkm_urns_qap2sum_sym,"Um")<attr(czkm_urns_qap2sum,"Um")){
    print("Symmetric diagram had better Um factor.")
    czkm_urns_qap2sum<-czek_matrix(soltysiakjaskulski_urns_use,original_diagram=TRUE,order=attr(czkm_urns_qap2sum_sym,"order"))
    attr(czkm_urns_qap2sum,"criterion_value")<-attr(czkm_urns_qap2sum_sym,"criterion_value")
    pdf("SoltysiakJaskulski_urns_qap2sum.pdf");plot(czkm_urns_qap2sum,plot_title="",label.cex=0.9);dev.off()
}else{
    if (attr(czkm_urns_qap2sum,"Um")<attr(czkm_urns_qap2sum_sym,"Um")){
        print("Asymmetric diagram had better Um factor.")
	czkm_urns_qap2sum_sym<-czek_matrix(soltysiakjaskulski_urns_use,original_diagram=FALSE,order=attr(czkm_urns_qap2sum,"order"))
        attr(czkm_urns_qap2sum_sym,"criterion_value")<-attr(czkm_urns_qap2sum,"criterion_value")
	pdf("SoltysiakJaskulski_urns_qap2sum_sym.pdf");plot(czkm_urns_qap2sum_sym,plot_title="",label.cex=0.9);dev.off()
    }
}
##=================================================================


## if it turned out that the Path length values are not the same
if (attr(czkm_urns_OLO_sym,"Path_length")<attr(czkm_urns_OLO,"Path_length")){
    print("Symmetric diagram had better OLO factor.")
    czkm_urns_OLO<-czek_matrix(soltysiakjaskulski_urns_use,original_diagram=TRUE,order=attr(czkm_urns_OLO_sym,"order"))
    attr(czkm_urns_OLO,"criterion_value")<-attr(czkm_urns_OLO_sym,"criterion_value")
    pdf("SoltysiakJaskulski_urns_OLO.pdf");plot(czkm_urns_OLO,plot_title="",label.cex=0.9);dev.off()
}else{
    if (attr(czkm_urns_OLO,"Path_length")<attr(czkm_urns_OLO_sym,"Path_length")){
	print("Asymmetric diagram had better OLO factor.")
	czkm_urns_OLO_sym<-czek_matrix(soltysiakjaskulski_urns_use,original_diagram=FALSE,order=attr(czkm_urns_OLO,"order"))
        attr(czkm_urns_OLO_sym,"criterion_value")<-attr(czkm_urns_OLO,"criterion_value")
	pdf("SoltysiakJaskulski_urns_OLO_sym.pdf");plot(czkm_urns_OLO_sym,plot_title="",label.cex=0.9);dev.off()
    }
}
##=================================================================
