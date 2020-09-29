## Analysis of Warzecha (2015)'s Internet availability data for some Silesian Voivodeship 
## counties data (basis for Fig. 3) found in RMaCzek as internet_availability.
## The produced results (ordering of observations) can depend on the random seed.
## Here, we load the saved random seed that resulted in the best orderings.


library(RMaCzek)
RNGversion("3.6.1") 
rexp(1)
#save(.Random.seed,file="KWar2015_randomseed.RData")
load("KWar2015_randomseed.RData")


## Internet availability data found in RMaCzek are obtained directly from Katarzyna Warzecha 
distancematrix_warzecha_internetavailability<-internet_availability$internet_availability_distances
warzecha_internetavailability_countynames_dataorder<-rownames(distancematrix_warzecha_internetavailability)
dist_warzecha_internetavailability<-as.dist(distancematrix_warzecha_internetavailability)

warzecha_internetavailability_countynames_proposedorder<-
c("bielski","Gliwice","Rybnik","pszczynski","Bielsko-Biala","Dabrowa_Gor.",
"bedzinski", "Chorzow","cieszynski", "Myslowice","Bytom",
"wodzislawski", "Czestochowa", "raciborski", "zywiecki", "zawiercianski",
"klobucki", "lubliniecki", "Tychy", "tarnogorski", "mikolowski",
"Katowice", "Jaworzno", "Zory", "Siemianowice_S.", "bierunsko-ledz.",
"Swietochlowice", "gliwicki", "czestochowski", "myszkowski",
"Piekary_Slaskie", "rybnicki", "Zabrze", "Jastrzebie-Zdroj",
"Ruda_Slaska", "Sosnowiec")

warzecha_internetavailability_proposedorder<-sapply(warzecha_internetavailability_countynames_proposedorder,function(x,warzecha_internetavailability_countynames_dataorder){which(x==warzecha_internetavailability_countynames_dataorder)},warzecha_internetavailability_countynames_dataorder=warzecha_internetavailability_countynames_dataorder,simplify=TRUE)
names(warzecha_internetavailability_proposedorder)<-NULL

czkm_internetavailability_orig<-czek_matrix(dist_warzecha_internetavailability,original_diagram=TRUE,order=warzecha_internetavailability_proposedorder)
czkm_internetavailability_orig_sym<-czek_matrix(dist_warzecha_internetavailability,original_diagram=FALSE,order=warzecha_internetavailability_proposedorder)
czkm_internetavailability_OLO<-czek_matrix(dist_warzecha_internetavailability,original_diagram=TRUE,order="OLO")
czkm_internetavailability_OLO_sym<-czek_matrix(dist_warzecha_internetavailability,original_diagram=FALSE,order="OLO")
czkm_internetavailability_qap2sum<-czek_matrix(dist_warzecha_internetavailability,original_diagram=TRUE,order="QAP_2SUM")
czkm_internetavailability_qap2sum_sym<-czek_matrix(dist_warzecha_internetavailability,original_diagram=FALSE,order="QAP_2SUM")



pdf("Warzecha_internetavailability_orig.pdf");plot(czkm_internetavailability_orig,plot_title="",label.cex=0.5);dev.off()
pdf("Warzecha_internetavailability_OLO.pdf");plot(czkm_internetavailability_OLO,plot_title="",label.cex=0.5);dev.off()
pdf("Warzecha_internetavailability_qap2sum.pdf");plot(czkm_internetavailability_qap2sum,plot_title="",label.cex=0.5);dev.off()
pdf("Warzecha_internetavailability_orig_sym.pdf");plot(czkm_internetavailability_orig_sym,plot_title="",label.cex=0.5);dev.off()
pdf("Warzecha_internetavailability_OLO_sym.pdf");plot(czkm_internetavailability_OLO_sym,plot_title="",label.cex=0.5);dev.off()
pdf("Warzecha_internetavailability_qap2sum_sym.pdf");plot(czkm_internetavailability_qap2sum_sym,plot_title="",label.cex=0.5);dev.off()

## if it turned out that the Um values are not the same
if (attr(czkm_internetavailability_qap2sum_sym,"Um")<attr(czkm_internetavailability_qap2sum,"Um")){
    print("Symmetric diagram had better Um factor.")
    czkm_internetavailability_qap2sum<-czek_matrix(dist_warzecha_internetavailability,original_diagram=TRUE,order=attr(czkm_internetavailability_qap2sum_sym,"order"))
    attr(czkm_internetavailability_qap2sum,"criterion_value")<-attr(czkm_internetavailability_qap2sum_sym,"criterion_value")
    pdf("Warzecha_internetavailability_qap2sum.pdf");plot(czkm_internetavailability_qap2sum,plot_title="",label.cex=0.5);dev.off()
}else{
    if (attr(czkm_internetavailability_qap2sum,"Um")<attr(czkm_internetavailability_qap2sum_sym,"Um")){
        print("Asymmetric diagram had better Um factor.")
        czkm_internetavailability_qap2sum_sym<-czek_matrix(dist_warzecha_internetavailability,original_diagram=FALSE,order=attr(czkm_internetavailability_qap2sum,"order"))
        attr(czkm_internetavailability_qap2sum_sym,"criterion_value")<-attr(czkm_internetavailability_qap2sum,"criterion_value")
        pdf("Warzecha_internetavailability_qap2sum_sym.pdf");plot(czkm_internetavailability_qap2sum_sym,plot_title="",label.cex=0.5);dev.off()
    }
}
##=================================================================


## if it turned out that the Path length values are not the same
if (attr(czkm_internetavailability_OLO_sym,"Path_length")<attr(czkm_internetavailability_OLO,"Path_length")){
    print("Symmetric diagram had better OLO factor.")
    czkm_internetavailability_OLO<-czek_matrix(dist_warzecha_internetavailability,original_diagram=TRUE,order=attr(czkm_internetavailability_OLO_sym,"order"))
    attr(czkm_internetavailability_OLO,"criterion_value")<-attr(czkm_internetavailability_OLO_sym,"criterion_value")
    pdf("Warzecha_internetavailability_OLO.pdf");plot(czkm_internetavailability_OLO,plot_title="",label.cex=0.5);dev.off()
}else{
    if (attr(czkm_internetavailability_OLO,"Path_length")<attr(czkm_internetavailability_OLO_sym,"Path_length")){
        print("Asymmetric diagram had better OLO factor.")
        czkm_internetavailability_OLO_sym<-czek_matrix(dist_warzecha_internetavailability,original_diagram=FALSE,order=attr(czkm_internetavailability_OLO,"order"))
        attr(czkm_internetavailability_OLO_sym,"criterion_value")<-attr(czkm_internetavailability_OLO,"criterion_value")
        pdf("Warzecha_internetavailability_OLO_sym.pdf");plot(czkm_internetavailability_OLO_sym,plot_title="",label.cex=0.5);dev.off()
    }
}
##=================================================================

## print out the ordering of the counties
print("Original ordering by Warzecha (2015) with factor values")
print(czkm_internetavailability_orig)
print("Ordering and factor values under OLO")
print(czkm_internetavailability_OLO)
print("Ordering and factor values under QAP_2SUM")
print(czkm_internetavailability_qap2sum)


column_order_stat_grouping<-c(8,10,12,16)
czkm_internetavailability_orig_largergroups<-czek_matrix(dist_warzecha_internetavailability,original_diagram=TRUE,order=warzecha_internetavailability_proposedorder,column_order_stat_grouping=column_order_stat_grouping)
czkm_internetavailability_OLO_largergroups<-czek_matrix(dist_warzecha_internetavailability,original_diagram=TRUE,order=attr(czkm_internetavailability_OLO,"order"),column_order_stat_grouping=column_order_stat_grouping)
czkm_internetavailability_qap2sum_largergroups<-czek_matrix(dist_warzecha_internetavailability,original_diagram=TRUE,order=attr(czkm_internetavailability_qap2sum,"order"),column_order_stat_grouping=column_order_stat_grouping)
pdf("Warzecha_internetavailability_orig_largergroups.pdf");plot(czkm_internetavailability_orig_largergroups,plot_title="",label.cex=0.5);dev.off()
pdf("Warzecha_internetavailability_OLO_largergroups.pdf");plot(czkm_internetavailability_OLO_largergroups,plot_title="",label.cex=0.5);dev.off()
pdf("Warzecha_internetavailability_qap2sum_largergroups.pdf");plot(czkm_internetavailability_qap2sum_largergroups,plot_title="",label.cex=0.5);dev.off()
