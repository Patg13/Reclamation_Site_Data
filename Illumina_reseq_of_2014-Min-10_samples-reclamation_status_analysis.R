#############################################################################
#############################################################################
#############################################################################
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("phyloseq")

getwd()
setwd("W:/Soil_Fungi_PES/FieldTrip2014/DNAsoil2014/Illumina reseq of 2014-Min-10 and 2016-natural-white-spruce samples")

	library(phyloseq)####
	library("ggplot2")####
	library("plyr")####
	library("vegan")####

#Env.#
env.mo<-read.csv("Illumina reseq of 2014-Min-10_data table.csv",header=T)
names(env.mo)
rownames(env.mo)<-env.mo$X
env.allsample<-env.mo
env.mo<-sample_data(env.mo[env.mo$TreeSp=="Jack Pine",])
env.mo$SampleShort

	distmat2<-cbind(env.mo$North_m,env.mo$East_m)
	#distmat2<-cbind(sample_data(GP1)$edgedis_m,sample_data(GP1)$y.dist)
	distmat2<-as.matrix(distmat2)
	distmat2
	pcnm2 <- pcnm(dist(distmat2))  #using easting/ northing

env.mo<-cbind(env.mo,pcnm2$vectors)	#launch code at very end
env.mo<-sample_data(env.mo)

#OTUs#
datamat<-read.csv("Analysis19471_ITSSoil_20Dec2017_15h23.new.besthits.with_biof.NO_SINGLETONS.blast.csv",header=T)
names(datamat)
datamat$Size<-rowSums(datamat[,rownames(env.mo)])
	#range(rowSums(datamat[,rownames(env.mo)]))
	#range(rowSums(datamat[,14:42]))
datamat<-datamat[datamat$Size>10,]
#datamat<-subset(datamat,datamat$Size>10)

		##Latin names##
		ALLspecies.dat<-aggregate(datamat[,14:42], list(datamat$ID,datamat$Trophic.function,datamat$Trophic.status,datamat$Lifestyle,datamat$Growth.form),sum)
		rownames(ALLspecies.dat)<-ALLspecies.dat[,1]
		ALLspecies.sub<-ALLspecies.dat[,rownames(env.mo)]
			species1314.otu<-otu_table(ALLspecies.sub, taxa_are_rows=TRUE, errorIfNULL=TRUE)

		#FOR Presence/absence#
		species1314.pres<-ALLspecies.sub
		species1314.pres[species1314.pres>0] <-1 
			species1314.otupres<-otu_table(species1314.pres, taxa_are_rows=TRUE, errorIfNULL=TRUE)

	#Trophic##for latin.names
	troph<-data.frame(ALLspecies.dat[,1:5])  
	colnames(troph)<-c("ID" ,"Trophic_function","Trophic_status","Lifestyle","Growth_form")
	tax1 <- tax_table(as.matrix(troph))
	colnames(tax1)
	taxa_names(tax1)<-ALLspecies.dat[,1]

colnames(species1314.otu)
rownames(env.mo)
head(species1314.otu)
tail(species1314.otu)

###
phydat<-merge_phyloseq(species1314.otu,env.mo,tax1)
	colnames(tax_table(phydat))
	colnames(sample_data(phydat))
	str(otu_table(phydat))
phydat.pres<-merge_phyloseq(species1314.otupres,env.mo,tax1)


		##Trophic status##
		ALLspecies.dat<-aggregate(datamat[,14:42], list(datamat$Trophic.function),sum)
		rownames(ALLspecies.dat)<-ALLspecies.dat[,1]
		ALLspecies.sub<-ALLspecies.dat[,rownames(env.mo)]
			species1314.otu<-otu_table(ALLspecies.sub, taxa_are_rows=TRUE, errorIfNULL=TRUE)
				#func.names<-c("noID","Animal.p","Endophyte"
				# ,"Myco.p","Myco.p.yeast","Path.yeast","Plant.p"
				# ,"Plant.p.yeast","Sap.","Sap.brown.rot"
				# ,"Sap.fac.yeast","Sap.p","Sap.white.rot"
				# ,"Sap.yeast","AM","EM","EM.white.rot","Ericoid","Lichen","Unknown")
				func.names<-c("noID","AP","EP"
				 ,"MP","MPY","PY","PP"
				 ,"PPY","ST","SBR"
				 ,"SFY","SP","SWR"
				 ,"SY","AM","EM","EMWR","ER","L","UK")
				cbind(func.names,rownames(species1314.otu))
				rownames(species1314.otu)<-func.names
			phydat.troph<-merge_phyloseq(species1314.otu,env.mo)


		##Trophic status: SUM of SPECIES presence / absence##
		ALLspecies.dat<-aggregate(datamat[,14:42], list(datamat$ID,datamat$Trophic.function,datamat$Trophic.status,datamat$Lifestyle,datamat$Growth.form),sum)
		rownames(ALLspecies.dat)<-ALLspecies.dat[,1]
		ALLspecies.sub<-ALLspecies.dat[,rownames(env.mo)]
		 #FOR Presence/absence#
		 species1314.pres<-ALLspecies.sub
		 species1314.pres[species1314.pres>0] <-1 
				#Trophic stat summed by species#
				ALLspecies.dat<-ALLspecies.dat[c(match(rownames(species1314.pres), ALLspecies.dat$Group.1)),]
				troph.stat<-ALLspecies.dat$Group.2	#for 'grouped-troph functions'
				species1314.foragg<-cbind(species1314.pres,troph.stat)
				trophbysumsp<-aggregate(species1314.foragg[,1:12], list(species1314.foragg$troph.stat),sum)
				trophbysumsp$Group.1
				species1314.mo<-data.frame(trophbysumsp[,2:13])
			levels(droplevels(trophbysumsp$Group.1))	#droplevels(trophbysumsp)$Group.1
			#rownames(species1314.mo)<-c("noID","Animal.p","Endophyte"
			# ,"Myco.p","Myco.p.yeast","Path.yeast","Plant.p"
			# ,"Plant.p.yeast","Sap.","Sap.brown.rot"
			# ,"Sap.fac.yeast","Sap.p","Sap.white.rot"
			# ,"Sap.yeast","AM","EM","EM.white.rot","Ericoid","Lichen","Unknown")
			rownames(species1314.mo)<-c("noID","AP","EP"
			 ,"MP","MPY","PY","PP"
			 ,"PPY","ST","SBR"
			 ,"SFY","SP","SWR"
			 ,"SY","AM","EM","EMWR","ER","L","UK")
			cbind(levels(droplevels(trophbysumsp$Group.1)),rownames(species1314.mo))
				#colnames(species1314.mo)<-rownames(env.mo)
			species1314.otu<-otu_table(species1314.mo, taxa_are_rows=TRUE, errorIfNULL=TRUE)
			phydat.trophpres<-merge_phyloseq(species1314.otu,env.mo)


#################################################################3
#TROPHdat: Pres/Abs or NOTtransformed dataset#
GP1.trophpres <- prune_taxa(taxa_sums(phydat.trophpres )>0, phydat.trophpres )

	####TROPHdat: transform dataset###	1914 taxa ##
	GP2.troph <- prune_taxa(taxa_sums(phydat.troph )>0, phydat.troph )
	GP1.troph = transform_sample_counts(GP2.troph, function(x) 1e+06 * x/sum(x))	

#################################################################3
#Pres/Abs or NOTtransformed dataset#
GP1.pres <- prune_taxa(taxa_sums(phydat.pres )>0, phydat.pres )

	###transform dataset###	1914 taxa ##
	GP2 <- prune_taxa(taxa_sums(phydat )>0, phydat )
	GP1 = transform_sample_counts(GP2, function(x) 1e+06 * x/sum(x))	

################################
GPr = transform_sample_counts(phydat , function(x) x/sum(x))
GPf = filter_taxa(GPr, function(x) var(x) > 1e-05, TRUE)
GPg = transform_sample_counts(GPf , function(x) 1e+06 * x)

#############################
	rec.site <- prune_samples(sample_data(GP1.pres)$ReclaimStat == "Reclaim", GP1.pres)
	rec.site <- prune_taxa(taxa_sums(rec.site )>0, rec.site )

	burn.site <- prune_samples(sample_data(GP1.pres)$ReclaimStat == "Forestry", GP1.pres)
	burn.site <- prune_taxa(taxa_sums(burn.site )>0, burn.site )

	intersect(rownames(otu_table(burn.site)),rownames(otu_table(rec.site)))

Notrec <- prune_samples(sample_data(GP1)$ReclaimStat != "Reclaim", GP1)
		rec2 <- prune_samples(sample_data(GP1)$ReclaimStat == "Reclaim", GP1)
		rec3 <- prune_taxa(taxa_sums(Notrec )==0, rec2 )
Notburn <- prune_samples(sample_data(GP1)$ReclaimStat != "Forestry", GP1)
		burn2 <- prune_samples(sample_data(GP1)$ReclaimStat == "Forestry", GP1)
		burn3 <- prune_taxa(taxa_sums(Notburn )==0, burn2 )

###############################################
levels(as.factor(tax_table(GP1)[,2]))
paste( rownames(otu_table(subset_taxa(GP2, Trophic_function %in% c("Biotroph/Animal parasite")))),
  rowSums(otu_table(subset_taxa(GP2, Trophic_function %in% c("Biotroph/Animal parasite")))))
paste( rownames(otu_table(subset_taxa(GP2, Trophic_function %in% c("Biotroph/Endophyte")))),
  rowSums(otu_table(subset_taxa(GP2, Trophic_function %in% c("Biotroph/Endophyte")))))
paste( rownames(otu_table(subset_taxa(GP2, Trophic_function %in% c("Biotroph/Mycoparasite")))),
  rowSums(otu_table(subset_taxa(GP2, Trophic_function %in% c("Biotroph/Mycoparasite")))))
paste( rownames(otu_table(subset_taxa(GP2, Trophic_function %in% c("Biotroph/Mycoparasite/yeast")))),
  rowSums(otu_table(subset_taxa(GP2, Trophic_function %in% c("Biotroph/Mycoparasite/yeast")))))
paste( rownames(otu_table(subset_taxa(GP2, Trophic_function %in% c("Biotroph/Pathogen/yeast")))),
  rowSums(otu_table(subset_taxa(GP2, Trophic_function %in% c("Biotroph/Pathogen/yeast")))))
paste( rownames(otu_table(subset_taxa(GP2, Trophic_function %in% c("Biotroph/Plant pathogen/yeast")))),
  rowSums(otu_table(subset_taxa(GP2, Trophic_function %in% c("Biotroph/Plant pathogen/yeast")))))

paste( rownames(otu_table(subset_taxa(GP2, Trophic_function %in% c("Saprotroph/Brown rot")))),
  rowSums(otu_table(subset_taxa(GP2, Trophic_function %in% c("Saprotroph/Brown rot")))))
paste( rownames(otu_table(subset_taxa(GP2, Trophic_function %in% c("Saprotroph/Pathogen")))),
  rowSums(otu_table(subset_taxa(GP2, Trophic_function %in% c("Saprotroph/Pathogen")))))
paste( rownames(otu_table(subset_taxa(GP2, Trophic_function %in% c("Saprotroph/White rot")))),
  rowSums(otu_table(subset_taxa(GP2, Trophic_function %in% c("Saprotroph/White rot")))))

paste( rownames(otu_table(subset_taxa(GP2, Trophic_function %in% c("Symbiotroph/Ectomycorrhizal/white rot")))),
  rowSums(otu_table(subset_taxa(GP2, Trophic_function %in% c("Symbiotroph/Ectomycorrhizal/white rot")))))
paste( rownames(otu_table(subset_taxa(GP2, Trophic_function %in% c("Symbiotroph/Lichenized")))),
  rowSums(otu_table(subset_taxa(GP2, Trophic_function %in% c("Symbiotroph/Lichenized")))))

	get_sample(GP2,i="Ramaria abietina 18S rRNA gene, 5.8S rRN")

rownames(otu_table(subset_taxa(GP2, Trophic_function %in% c("Saprotroph/Facultative yeast"))))
rownames(otu_table(subset_taxa(GP2, Trophic_function %in% c("Saprotroph/Yeast"))))

rownames(otu_table(subset_taxa(GP2, Trophic_function %in% c("Symbiotroph/Arbuscular mycorrhizal"))))
a<-(subset_taxa(GP2, Trophic_function %in% c("")))
b<-(subset_taxa(GP2, Trophic_function %in% c("Unknown Function")))

################################################
################################################
###Taxa unique to a tax_table label###
	#str(tax_table(GP1)[,3])
	get_taxa_unique(GP1, "Trophic_function")
	get_taxa_unique(GP1, "Trophic_status")
	get_taxa_unique(GP1, "Lifestyle")
		subset_taxa(GP1, Trophic_status %in% c("Saprotroph"))
		  subset_taxa(GP1, Lifestyle %in% c("Saprotroph"))
		  subset_taxa(GP1, Lifestyle %in% c("Facultative yeast"))
		  subset_taxa(GP1, Lifestyle %in% c("White rot"))
		  subset_taxa(GP1, Lifestyle %in% c("Brown rot"))
		  subset_taxa(GP1, Lifestyle %in% c("Yeast"))

		subset_taxa(GP1, Trophic_status %in% c("Symbiotroph"))  
		  g<-subset_taxa(GP1, Lifestyle %in% c("Arbuscular mycorrhizal"))
		  subset_taxa(burn.site, Lifestyle %in% c("Ectomycorrhizal"))
		  subset_taxa(GP1, Lifestyle %in% c("Ericoid"))
		  subset_taxa(GP1, Lifestyle %in% c("Lichenized"))

		subset_taxa(GP1, Trophic_status %in% c("Biotroph"))
		  subset_taxa(GP1, Lifestyle %in% c("Animal parasite"))
		  subset_taxa(GP1, Lifestyle %in% c("Plant pathogen"))
		  subset_taxa(GP1, Lifestyle %in% c("Mycoparasite"))
		  subset_taxa(GP1, Lifestyle %in% c("Pathogen"))
		  subset_taxa(GP1, Lifestyle %in% c("Endophyte"))

		  subset_taxa(GP1, Lifestyle %in% c(""))
		  subset_taxa(GP1, Lifestyle %in% c("Unknown Function"))

	get_sample(GP1,i="Placynthiella icmalea 18S ribosomal RNA ")
grep("Cenococcum*",rownames(otu_table(GP1)),value=T)

grep("Amphinema*",rownames(otu_table(burn.site)),value=T)
grep("Amphinema*",rownames(otu_table(rec.site)),value=T)


################################################
	  ###Ubiquitous taxa across all samples###
	  GPall.pres <- prune_taxa(taxa_sums(GP1.pres )>=12, GP2 )
	  paste(rownames(otu_table(GPall.pres)),rowSums(otu_table(GPall.pres)),tax_table(GPall.pres)[,2])

################################################
#############################################################################
#######Correlation matrix######
#############################################################################
corProb <- function(X, dfr = nrow(X) - 2) { 
    R <- cor(X) 
    above <- row(R) < col(R) 
    r2 <- R[above]^2 
    Fstat <- r2 * dfr / (1 - r2) 
    R[above] <- 1 - pf(Fstat, 1, dfr) 
    class(R) <- "corProb" 
    R } 
print.corProb <- function(x, digits = getOption("digits"), quote = FALSE, na.print = "", 
    justify = "none") { 
    xx <- format(unclass(round(x, digits = 4)), digits = digits, justify = justify) 
    if (any(ina <- is.na(x))) 
        xx[ina] <- na.print 
    cat("\nCorrelations are shown below the diagonal\n") 
    cat("P-values are shown above the diagonal\n\n") 
    print(xx, quote = quote) 
    invisible(x) } 

corProb(cbind(sample_data(GP1)$OrgDep_cm
	,sample_data(GP1)$min_NO3
	,sample_data(GP1)$min_NH4
	,sample_data(GP1)$min_totalN
	,sample_data(GP1)$min_P
	,sample_data(GP1)$PCNM1
	,sample_data(GP1)$PCNM2
	,sample_data(GP1)$PCNM3))


####################################
####################
##Network analyses##= "Facultative yeast" = "darkgoldenrod4",
####################
collifestyle <- c(" "="pink","Biotroph/Animal parasite" = "deepskyblue2","Symbiotroph/Arbuscular mycorrhizal" = "coral","Symbiotroph/Ectomycorrhizal"  = "coral3","Biotroph/Endophyte"="bisque1", "Symbiotroph/Ericoid" = "coral2", "Symbiotroph/Lichenized"="green","Biotroph/Mycoparasite"="dodgerblue1","Saprotroph/Pathogen"="blue",
	"Biotroph/Plant pathogen"="cornflowerblue","Saprotroph"="darkgoldenrod1","Unknown Function"="purple","Saprotroph/Brown rot"="bisque3","Saprotroph/Facultative yeast"="darkgoldenrod4","Symbiotroph/Ectomycorrhizal/white rot"="white","Saprotroph/White rot"="white","Biotroph/Mycoparasite/yeast"="grey10","Biotroph/Pathogen/yeast"="grey50", "Biotroph/Plant pathogen/yeast"="grey80","Saprotroph/Yeast"="darkgoldenrod3")
#symbvec<- c("Min10"=22,"Org"=21,"Root"=23,"Min35"=24)#use symbols for soil_layer (Min=circle vs Org=square vs Seedling= diamond vs Soil=triangle)
symbvec<- c("Min10"=15,"Org"=16,"Root"=18,"Min35"=17)#use symbols for soil_layer (Min=circle vs Org=square vs Seedling= diamond vs Soil=triangle)
#symbvec<- c("Jack Pine"=19,"Wt. Spruce"=15,"Sib.Larch"=17)#use symbols for soil_layer (JPine=circle vs Wtspruce=square vs sib.larch=triangle)
#cols <- c("-5" = "grey50", "0" = "grey80","5" = "#FEE5D9","15" = "#FCAE91","35" = "#FB6A4A","50" = "#CB181D")
colvec <- c("Forestry" = "grey20","Reclaim"= "grey80")

	#plot_net#
	GP100 = prune_taxa(names(sort(taxa_sums(GP1), TRUE))[1:100], GP1)
	plot_net(GP100 ,type="taxa",distance = "bray",color="Trophic_function",point_label="ID",laymeth="fruchterman.reingold", maxdist=0.2)#
	p<-plot_net(GP1 ,type="samples",distance = "bray",point_label="TreeSp",color= "DistType", maxdist=0.7,laymeth="fruchterman.reingold",hjust=0.5)#
	p + scale_colour_manual(values = colvec, name = "DistType")+scale_shape_manual(values = symbvec, name = "Type", label=c("Min10","Org","Root","Min35"))


	#plotting_network#
	tg <- make_network( rec.site  , "taxa", "bray", 0.1)
	p<-plot_network(tg,  rec.site , "taxa", line_weight = 0.4,hjust=1.15,color="Trophic_function",label =NULL)
	p + scale_colour_manual(values = collifestyle)
		tg <- make_network( burn.site  , "taxa", "bray", 0.1)
		p<-plot_network(tg,  burn.site, "taxa", line_weight = 0.4,hjust=1.15,color="Trophic_function",label =NULL)
		p + scale_colour_manual(values = collifestyle)

	sg <- make_network(GP1 , "samples", "bray", 0.71,laymeth="fruchterman.reingold")
	p<-plot_network(sg, GP1 , "samples",label="NA",color="ReclaimStat",line_weight = 1,hjust=-0.5, point_size=5, line_alpha =0.5)	#shape="Type"
	p +theme(legend.position="right", legend.title = element_text(colour="black", size=25,face="bold"), legend.text = element_text(colour="black", size=23))+geom_text(label=p$data$SiteFull,size=4,vjust=1.7)+ scale_colour_manual(values = colvec, name = "Stand type",label=c("Forestry","Reclamation"))
		#+scale_shape_manual(values = symbvec, name = "Soil fraction", label=c("Min.soil & roots","Forest floor","Roots","Mineral soil"))


#########################################################################
##BEST HEATMAPs##########################################################
#########################################################################
	p<-plot_heatmap(GPg , "NMDS", "bray", "SiteFull")
		p$scales$scales[[2]]$name <- "Counts per million reads" #"Normalized counts" #
		p$scales$scales[[1]]$name <- "Coarse soils (Illumina)"
			print(p+ylab("Fungal taxa"))

	###Most abundant 250 subset###
	gp1heat <- prune_taxa(names(sort(taxa_sums(GP1 ),TRUE)[1:250]), GP1 )
	p <- plot_heatmap(GP1 , "NMDS", "bray", "SiteFull")
		p$scales$scales[[2]]$name <- "Counts per million reads"
		p$scales$scales[[1]]$name <- "All soil fractions"
			print(p+ylab("Fungal taxa"))

	#############################
	###Trophic status: fungi>0###
	p <- plot_heatmap(GP1.troph , "NMDS", "bray", "SiteFull")
		p$scales$scales[[2]]$name <- "Counts per million reads"
		p$scales$scales[[1]]$name <- "Coarse soils (Illumina)"
			print(p+ylab("Functional groups"))



#########################################################################
##CCA and RDAs##########################################################
#########################################################################
	dca.litter.comp<-decorana(t(otu_table(GP1.trophpres)))
	dca.litter.comp	#library norm data -->gradient more than 2.5 therefore use CCA
	plot(dca.litter.comp

#full/reduced models: TAXA: rel.abun##
mod<-cca(t(otu_table(GP1))~sample_data(GP1)$PCNM1+sample_data(GP1)$PCNM2+sample_data(GP1)$PCNM3+sample_data(GP1)$PCNM4+sample_data(GP1)$PCNM5+sample_data(GP1)$PCNM6,scale=TRUE)
mod<-cca(t(otu_table(GP1))~sample_data(GP1)$min_NO3+sample_data(GP1)$min_NH4+sample_data(GP1)$min_totalN+sample_data(GP1)$min_P,scale=TRUE)
  mod<-cca(t(otu_table(GP1))~sample_data(GP1)$PCNM1+sample_data(GP1)$PCNM3,scale=TRUE)	
  mod<-cca(t(otu_table(GP1))~sample_data(GP1)$min_NO3+sample_data(GP1)$min_P,scale=TRUE)	#NO3 borderline sig by marginal test

	#full/reduced models: TAXA: pres/abs##
	mod<-cca(t(otu_table(GP1.pres))~sample_data(GP1.pres)$PCNM1+sample_data(GP1.pres)$PCNM2+sample_data(GP1.pres)$PCNM3+sample_data(GP1.pres)$PCNM4+sample_data(GP1.pres)$PCNM5+sample_data(GP1.pres)$PCNM6,scale=TRUE)
	mod<-cca(t(otu_table(GP1.pres))~sample_data(GP1.pres)$min_NO3+sample_data(GP1.pres)$min_NH4+sample_data(GP1.pres)$min_totalN+sample_data(GP1.pres)$min_P,scale=TRUE)
	  mod<-cca(t(otu_table(GP1.pres))~sample_data(GP1.pres)$PCNM1,scale=TRUE)
	  mod<-cca(t(otu_table(GP1.pres))~sample_data(GP1.pres)$min_P,scale=TRUE)

#full/reduced models: Trophic function: rel.abun##
mod<-rda(t(otu_table(GP1.troph))~sample_data(GP1.troph)$PCNM1+sample_data(GP1.troph)$PCNM2+sample_data(GP1.troph)$PCNM3+sample_data(GP1.troph)$PCNM4+sample_data(GP1.troph)$PCNM5+sample_data(GP1.troph)$PCNM6,scale=TRUE)
mod<-rda(t(otu_table(GP1.troph))~sample_data(GP1.troph)$min_NO3+sample_data(GP1.troph)$min_NH4+sample_data(GP1.troph)$min_totalN+sample_data(GP1.troph)$min_P,scale=TRUE)
  mod<-rda(t(otu_table(GP1.troph))~sample_data(GP1.troph)$PCNM1+sample_data(GP1.troph)$PCNM2,scale=TRUE)
  mod<-rda(t(otu_table(GP1.troph))~sample_data(GP1.troph)$min_NH4+sample_data(GP1.troph)$min_P,scale=TRUE)

	#full/reduced models: Trophic function: pres/abs##
	mod<-rda(t(otu_table(GP1.trophpres))~sample_data(GP1.trophpres)$PCNM1+sample_data(GP1.trophpres)$PCNM2+sample_data(GP1.trophpres)$PCNM3+sample_data(GP1.trophpres)$PCNM4+sample_data(GP1.trophpres)$PCNM5+sample_data(GP1.trophpres)$PCNM6,scale=TRUE)
	mod<-rda(t(otu_table(GP1.trophpres))~sample_data(GP1.trophpres)$min_NO3+sample_data(GP1.trophpres)$min_NH4+sample_data(GP1.trophpres)$min_totalN+sample_data(GP1.trophpres)$min_P,scale=TRUE)
	  mod<-rda(t(otu_table(GP1.trophpres))~sample_data(GP1.trophpres)$PCNM1,scale=TRUE)
	  mod<-rda(t(otu_table(GP1.trophpres))~sample_data(GP1.trophpres)$min_totalN,scale=TRUE)

	ordistep(mod, direction =  "both", Pin = 0.05, Pout = 0.1, pstep = 100, perm.max = 1000, steps = 50, trace = TRUE)

		anova.cca(mod)
		anova.cca(mod,by="margin",permu=999)
		anova.cca(mod,by="term",permu=999)
		mod
		#summary(mod,by="term",permu=999)


#############
#TAXA: rel.abun##
mod1<-cca(t(otu_table(GP1))~sample_data(GP1)$ReclaimStat,scale=TRUE)
mod<-cca(t(otu_table(GP1))~sample_data(GP1)$RegenYr,scale=TRUE)
mod<-cca(t(otu_table(GP1))~sample_data(GP1)$OrgDep_cm,scale=TRUE)
  #mod<-cca(t(otu_table(GP1))~sample_data(GP1)$PCNM1+sample_data(GP1)$PCNM3,scale=TRUE)	
  mod<-cca(t(otu_table(GP1))~sample_data(GP1)$pcnm1_jp,scale=TRUE)	
  mod<-cca(t(otu_table(GP1))~sample_data(GP1)$min_NO3+sample_data(GP1)$min_P,scale=TRUE)	#NO3 borderline sig by marginal test

	#TAXA: pres/abs##
	mod2<-cca(t(otu_table(GP1.pres))~sample_data(GP1.pres)$ReclaimStat,scale=TRUE)
	mod<-cca(t(otu_table(GP1.pres))~sample_data(GP1.pres)$RegenYr,scale=TRUE)
	mod<-cca(t(otu_table(GP1.pres))~sample_data(GP1.pres)$OrgDep_cm,scale=TRUE)
	  #mod<-cca(t(otu_table(GP1.pres))~sample_data(GP1.pres)$PCNM1,scale=TRUE)
	  mod<-cca(t(otu_table(GP1.pres))~sample_data(GP1.pres)$pcnm1_jp,scale=TRUE)
	  mod<-cca(t(otu_table(GP1.pres))~sample_data(GP1.pres)$min_P,scale=TRUE)

#Trophic function: rel.abun##
mod3<-rda(t(otu_table(GP1.troph))~sample_data(GP1.troph)$ReclaimStat,scale=TRUE)
mod<-rda(t(otu_table(GP1.troph))~sample_data(GP1.troph)$RegenYr,scale=TRUE)
mod<-rda(t(otu_table(GP1.troph))~sample_data(GP1.troph)$OrgDep_cm,scale=TRUE)
  #mod<-rda(t(otu_table(GP1.troph))~sample_data(GP1.troph)$PCNM1+sample_data(GP1.troph)$PCNM2,scale=TRUE)
  mod<-rda(t(otu_table(GP1.troph))~sample_data(GP1.troph)$pcnm1_jp+sample_data(GP1.troph)$pcnm2_jp,scale=TRUE)
  mod<-rda(t(otu_table(GP1.troph))~sample_data(GP1.troph)$min_NH4+sample_data(GP1.troph)$min_P,scale=TRUE)

	#Trophic function: pres/abs##
	mod4<-rda(t(otu_table(GP1.trophpres))~sample_data(GP1.trophpres)$ReclaimStat,scale=TRUE)
	mod<-rda(t(otu_table(GP1.trophpres))~sample_data(GP1.trophpres)$RegenYr,scale=TRUE)
	mod<-rda(t(otu_table(GP1.trophpres))~sample_data(GP1.trophpres)$OrgDep_cm,scale=TRUE)
	  #mod<-rda(t(otu_table(GP1.trophpres))~sample_data(GP1.trophpres)$PCNM1,scale=TRUE)
	  mod<-rda(t(otu_table(GP1.trophpres))~sample_data(GP1.trophpres)$pcnm1_jp,scale=TRUE)
	  mod<-rda(t(otu_table(GP1.trophpres))~sample_data(GP1.trophpres)$min_totalN,scale=TRUE)

		anova.cca(mod)
		anova.cca(mod,by="margin",permu=999)
		anova.cca(mod,by="term",permu=999)
		summary(mod)


###Conditioned models### #remove ffdepth, min_P, PCNM1 since significantly different between reclamation status
#TAXA: rel.abun##
#mod<-cca(t(otu_table(GP1))~sample_data(GP1)$ReclaimStat+Condition(sample_data(GP1)$OrgDep_cm),scale=TRUE)
#mod<-cca(t(otu_table(GP1))~sample_data(GP1)$ReclaimStat+Condition(sample_data(GP1)$pcnm1_jp),scale=TRUE)	
#mod<-cca(t(otu_table(GP1))~sample_data(GP1)$ReclaimStat+Condition(sample_data(GP1)$min_NO3+sample_data(GP1)$min_P),scale=TRUE)	#NO3 borderline sig by marginal test
 mod1<-cca(t(otu_table(GP1))~sample_data(GP1)$ReclaimStat+Condition(sample_data(GP1)$min_NO3),scale=TRUE)	#NO3 borderline sig by marginal test

	#TAXA: pres/abs##
	#mod<-cca(t(otu_table(GP1.pres))~sample_data(GP1.pres)$ReclaimStat+Condition(sample_data(GP1.pres)$OrgDep_cm),scale=TRUE)
	#mod<-cca(t(otu_table(GP1.pres))~sample_data(GP1.pres)$ReclaimStat+Condition(sample_data(GP1.pres)$pcnm1_jp),scale=TRUE)
	#mod<-cca(t(otu_table(GP1.pres))~sample_data(GP1.pres)$ReclaimStat+Condition(sample_data(GP1.pres)$min_P),scale=TRUE)
	 mod2<-cca(t(otu_table(GP1.pres))~sample_data(GP1.pres)$ReclaimStat,scale=TRUE)

#Trophic function: rel.abun##
#mod<-rda(t(otu_table(GP1.troph))~sample_data(GP1.troph)$ReclaimStat+Condition(sample_data(GP1.troph)$OrgDep_cm),scale=TRUE)
#mod<-rda(t(otu_table(GP1.troph))~sample_data(GP1.troph)$ReclaimStat+Condition(sample_data(GP1.troph)$pcnm1_jp+sample_data(GP1.troph)$pcnm2_jp),scale=TRUE)
#mod<-rda(t(otu_table(GP1.troph))~sample_data(GP1.troph)$ReclaimStat+Condition(sample_data(GP1.troph)$min_NH4+sample_data(GP1.troph)$min_P),scale=TRUE)
 mod3<-rda(t(otu_table(GP1.troph))~sample_data(GP1.troph)$ReclaimStat+Condition(sample_data(GP1.troph)$min_NH4+sample_data(GP1.troph)$pcnm2_jp),scale=TRUE)

	#Trophic function: pres/abs##
	#mod<-rda(t(otu_table(GP1.trophpres))~sample_data(GP1.trophpres)$ReclaimStat+Condition(sample_data(GP1.trophpres)$OrgDep_cm),scale=TRUE)
	#mod<-rda(t(otu_table(GP1.trophpres))~sample_data(GP1.trophpres)$ReclaimStat+Condition(sample_data(GP1.trophpres)$pcnm1_jp),scale=TRUE)
	 mod4<-rda(t(otu_table(GP1.trophpres))~sample_data(GP1.trophpres)$ReclaimStat+Condition(sample_data(GP1.trophpres)$min_totalN),scale=TRUE)


		anova.cca(mod)
		anova.cca(mod,by="axis",permu=999)
		mod4

########################################
###plotting ordinations#################
########################################
scl <- 0 ## scl=0 (unscaled raw scores)#scl=2 (species scores scaled by eigenvalue) #scaling = 3 (both scecies and site scores are scaled symmetrically by square root of eigenvalue)
colvec <- c("black")
symbvec<- c("Min10"=21,"Root"=24)#use symbols for soil_layer (Min=circle vs Org=square vs Seedling= diamond vs Soil=triangle)
#symbvec<- c("Min10"=21,"Org"=22,"Min35"=23,"Root"=24)#use symbols for soil_layer (Min=circle vs Org=square vs Seedling= diamond vs Soil=triangle)
#cexvec<-c("none"=2.22,"Reclaim"=1.8) # reclaim versus natural stand 
#cexvec<-1
colorSp <- c("Forestry" = "black", "Reclaim" = "grey")
#colorSp <- c("Jack Pine" = "black", "Sib. Larch" = "grey","Wt. Spruce" = "white")

par(mfrow=c(2,2))	
#par(oma=c(0,0,0,0))
par(mar=c(5.1,4.1,0.5,2.1))

mintotal.mod<-mod4

minfit<-scores(mintotal.mod, display = "sites", shrink = FALSE, choices=c(1)) 
		wss <- (nrow(minfit)-1)*sum(apply(minfit,2,var))
		for (i in 2:10) wss[i] <- sum(kmeans(minfit, 
	  		 centers=i)$withinss)
		#plot(1:10, wss, type="b", xlab="Number of Clusters",
	  	#	ylab="Within groups sum of squares")
	mintotal.kclust<-kmeans(minfit, centers=2, iter.max = 10, nstart = 25)
		mintotal.kclust$cluster
	
	plot(mintotal.mod, type = "n", scaling = scl ,xlab="",ylab="")
		    #text(mintotal.mod,choices=c(1,2), display = "spec", cex = 0.8, scaling = scl, pch=3, col="red")
		    #text(mintotal.mod,choices=c(1,2), display = "spec", cex = 1, scaling = scl, pch=3, col="red")
		    points(mintotal.mod,choices=c(1,2), display = "spec", cex = 0.8, scaling = scl, pch=3, col="red")
		  	#text(mintotal.mod,choices=c(1,2), display = "cn", cex=1, col="darkred")

	ordihull(mintotal.mod, mintotal.kclust$cluster, draw = "lines",	##plotting clusters  
     		    col = "black",scaling = scl )

	    #GP1#color=reclaim, size = regen yr, pch = stand stg
	    with(sample_data(GP1), points(mintotal.mod,choices=c(1,2), display = "sites", col = "black",
	      scaling = scl , pch = symbvec[Type], bg = colorSp[ReclaimStat],cex=2))
			#transect.lab#
	   		# with(sample_data(GP1),points(mintotal.mod,choices=c(1,2),display="sites",col="black",
			#  scaling = scl , pch = c("14","2","6","7","8")[SiteShort],cex=2))

	#legend("bottomright",pch=c(21,21,  22,21,23,24)
	#,col="black",pt.bg=c("black","grey"	,"black","black","black","black"),
	#legend=c("Forestry","Reclamation",		"Forest floor","Coarse soil","Fine soil","Roots"),cex=1
	#,pt.cex=c(1.9,1.9,  1.9,1.9,1.9,1.9),bty="t",bg="white")  

	legend("bottomright",pch=c(21,21)#,  21,24)
	,col="black",pt.bg=c("black","grey")#,"black","black")
	,legend=c("Forestry","Reclamation")#,"Coarse soils","Roots")
	,cex=1,pt.cex=c(1.9,1.9)#,  1.9,1.9)
	,bty="t",bg="white") 


  specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
#title(main="",xlab = paste("CA1 ",specify_decimal(mintotal.mod$CA$eig[1]/sum(mintotal.mod$CA$eig)*100,2),"%"), ylab = paste("CA2 ",specify_decimal(mintotal.mod$CA$eig[2]/sum(mintotal.mod$CA$eig)*100,2),"%"),cex.lab=1.2,line=2.5)  
title(main="",xlab = paste("CCA1 ",specify_decimal(mintotal.mod$CCA$tot.chi/mintotal.mod$tot.chi*100,2),"%"), ylab = paste("CA1 ",specify_decimal(mintotal.mod$CA$eig[1]/sum(mintotal.mod$CA$eig)*100,2),"%"),cex.lab=1.2,line=2.5)  
#title(main=paste(specify_decimal(mintotal.mod$CCA$tot.chi/mintotal.mod$tot.chi*100,2),"% of total inertia"),xlab = paste("CCA1 ",specify_decimal(mintotal.mod$CCA$eig[1]/sum(mintotal.mod$CCA$eig)*100,2),"%"), ylab = paste("CCA2 ",specify_decimal(mintotal.mod$CCA$eig[2]/sum(mintotal.mod$CCA$eig)*100,2),"%"),cex.lab=1.2,line=2.5)  
#for RDA
#title(main="",xlab = paste("PC1 ",specify_decimal(mintotal.mod$CA$eig[1]/sum(mintotal.mod$CA$eig)*100,2),"%"), ylab = paste("PC2 ",specify_decimal(mintotal.mod$CA$eig[2]/sum(mintotal.mod$CA$eig)*100,2),"%"),cex.lab=1.2,line=2.5)  
title(main="",xlab = paste("RDA1 ",specify_decimal(mintotal.mod$CCA$tot.chi/mintotal.mod$tot.chi*100,2),"%"), ylab = paste("PC1 ",specify_decimal(mintotal.mod$CA$eig[1]/sum(mintotal.mod$CA$eig)*100,2),"%"),cex.lab=1.2,line=2.5)  
#title(main=paste(specify_decimal(mintotal.mod$CCA$tot.chi/mintotal.mod$tot.chi*100,2),"% of total inertia"),xlab = paste("RDA1 ",specify_decimal(mintotal.mod$CCA$eig[1]/sum(mintotal.mod$CCA$eig)*100,2),"%"), ylab = paste("RDA2 ",specify_decimal(mintotal.mod$CCA$eig[2]/sum(mintotal.mod$CCA$eig)*100,2),"%"),cex.lab=1.2,line=2.5)  

	#reclaimStat +cond.(non-corr. factors)	#scl=0 
	 text(2.5,-2.4,"A",cex=4)	
	 text(-2.4,-3,"B",cex=4)	
	 text(0.45,0.39,"C",cex=4)	
	 text(-0.69,0.37,"D",cex=4)	

		#reclaimStat [NO conditions]	#scl=0 
		 text(1.9,1.6,"A",cex=4)	
		 text(-2.4,2.6,"B",cex=4)	
		 text(0.41,0.34,"C",cex=4)	
		 text(-0.73,0.47,"D",cex=4)	



########################################################################
########################################################################
###upregulated taxa###
###Up/down reg by species###
mod.tax<-cca(t(otu_table(GP1))~sample_data(GP1)$ReclaimStat+Condition(sample_data(GP1)$min_NO3),scale=TRUE)	#NO3 borderline sig by marginal test
mod.pres<-cca(t(otu_table(GP1.pres))~sample_data(GP1.pres)$ReclaimStat,scale=TRUE)
mod.tax1<-cca(t(otu_table(GP1))~sample_data(GP1)$ReclaimStat,scale=TRUE)	

 
 	tax.reclaim<-attributes(mod.tax$CCA$v[mod.tax$CCA$v[,1]< -0.5,])	#need space b/w < and -
	  GP1.tax.reclaim<- prune_taxa( unlist(tax.reclaim) , GP1.pres)
	  ems.tax.reclaim<-subset_taxa(GP1.tax.reclaim, Lifestyle %in% c("Ectomycorrhizal"))
 		tax1.reclaim<-attributes(mod.tax1$CCA$v[mod.tax1$CCA$v[,1]< -0.5,])	#need space b/w < and -
		  GP1.tax1.reclaim<- prune_taxa( unlist(tax1.reclaim) , GP1.pres)
	 	  ems.tax1.reclaim<-subset_taxa(GP1.tax1.reclaim, Lifestyle %in% c("Ectomycorrhizal"))
 	  match.tax.reclaim<- prune_taxa( unlist(tax1.reclaim) , GP1.tax.reclaim)

	tax.forestry<-attributes(mod.tax$CCA$v[mod.tax$CCA$v[,1]> 0.5,])	#need space b/w < and -
	  GP1.tax.forestry<- prune_taxa( unlist(tax.forestry) , GP1.pres)
	  ems.tax.forestry<-subset_taxa(GP1.tax.forestry, Lifestyle %in% c("Ectomycorrhizal"))
		tax1.forestry<-attributes(mod.tax1$CCA$v[mod.tax1$CCA$v[,1]> 0.5,])	#need space b/w < and -
		  GP1.tax1.forestry<- prune_taxa( unlist(tax1.forestry) , GP1.pres)
		  ems.tax1.forestry<-subset_taxa(GP1.tax1.forestry, Lifestyle %in% c("Ectomycorrhizal"))
 	  match.tax.forestry<- prune_taxa( unlist(tax1.forestry) , GP1.tax.forestry)

		pres.reclaim<-attributes(mod.pres$CCA$v[mod.pres$CCA$v[,1]< -0.5,])	#need space b/w < and -
		  GP1.pres.reclaim<- prune_taxa( unlist(pres.reclaim) , GP1.pres)
		  ems.pres.reclaim<-subset_taxa(GP1.pres.reclaim, Lifestyle %in% c("Ectomycorrhizal"))
	#		pres1.reclaim<-attributes(mod.pres1$CCA$v[mod.pres1$CCA$v[,1]< -0.5,])	#need space b/w < and -
	#		  GP1.pres1.reclaim<- prune_taxa( unlist(pres1.reclaim) , GP1.pres)
	#	 	  ems.pres1.reclaim<-subset_taxa(GP1.pres1.reclaim, Lifestyle %in% c("Ectomycorrhizal"))
	# 	  match.pres.reclaim<- prune_taxa( unlist(pres1.reclaim) , GP1.pres.reclaim)

		pres.forestry<-attributes(mod.pres$CCA$v[mod.pres$CCA$v[,1]> 0.5,])	#need space b/w < and -
		  GP1.pres.forestry<- prune_taxa( unlist(pres.forestry) , GP1.pres)
		  ems.pres.forestry<-subset_taxa(GP1.pres.forestry, Lifestyle %in% c("Ectomycorrhizal"))
	#		pres1.forestry<-attributes(mod.pres1$CCA$v[mod.pres1$CCA$v[,1]> 0.5,])	#need space b/w < and -
	#		  GP1.pres1.forestry<- prune_taxa( unlist(pres1.forestry) , GP1.pres)
	#		  ems.pres1.forestry<-subset_taxa(GP1.pres1.forestry, Lifestyle %in% c("Ectomycorrhizal"))
	# 	  match.pres.forestry<- prune_taxa( unlist(pres1.forestry) , GP1.pres.forestry)

	 	  match.both.reclaim <- prune_taxa( rownames(otu_table( GP1.pres.reclaim )) ,  match.tax.reclaim )
	 	  match.both.forestry <- prune_taxa( rownames(otu_table(GP1.pres.forestry )) ,  match.tax.forestry )

	 	  #match.both.reclaim <- prune_taxa( rownames(otu_table( GP1.pres.reclaim )) ,  GP1.tax.reclaim )
	 	  #match.both.forestry <- prune_taxa( rownames(otu_table(GP1.pres.forestry )) ,  GP1.tax.forestry )

	taxaSize<-aggregate(datamat[,2], list(datamat$ID),sum)
	rownames(taxaSize)<-taxaSize[,1]
reclaim <-taxaSize[rownames(otu_table(match.both.reclaim )),]
	trim.reclaim <-subset(reclaim ,reclaim $x>60)
	tax_table(match.both.reclaim )[rownames(trim.reclaim ),2]
forestry <-taxaSize[rownames(otu_table(match.both.forestry )),]
	trim.forestry <-subset(forestry ,forestry $x>60)
	tax_table(match.both.forestry )[rownames(trim.forestry ),2]


#####################################################################
########Diversity INDEXES############################################
#####################################################################
div.dat<-estimate_richness(GP2)
div.dat
rownames(div.dat)
rownames(env.mo)
colnames(env.mo)
div.dat<-cbind(div.dat,env.mo)
colnames(div.dat)

t.test(div.dat$Observed~div.dat$ReclaimStat)
t.test(div.dat$Chao1~div.dat$ReclaimStat)
t.test(div.dat$ACE~div.dat$ReclaimStat)
t.test(div.dat$Shannon~div.dat$ReclaimStat)
t.test(div.dat$Simpson~div.dat$ReclaimStat)
t.test(div.dat$InvSimpson~div.dat$ReclaimStat)
t.test(div.dat$Fisher~div.dat$ReclaimStat)



