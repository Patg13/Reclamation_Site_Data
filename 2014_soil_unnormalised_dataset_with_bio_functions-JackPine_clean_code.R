setwd("W:/Soil_Fungi_PES/FieldTrip2014/DNAsoil2014")


	library(phyloseq)####
	library("ggplot2")####
	library("plyr")####
	library("vegan")####
###################################################
env.mat<-read.csv("Stat-tests_2013-14 site descriptionForestReclam.csv",header=T)
names(env.mat)
env.mat$Sample
sitemerge.env1<-c(paste(env.mat$Sample,sep=""))

	env.mo<-subset(env.mat, env.mat$SiteShort=="M2"|env.mat$SiteShort=="M7"&!env.mat$Sample=="M07_1LM10_2014"|env.mat$SiteShort=="TR1P"|env.mat$SiteShort=="TR1B")
	head(env.mo)
	sitemerge.env2<-c(paste(env.mo$Sample, sep=""))
	titles.env<-c(paste(env.mo$SiteShort,env.mo$Type,env.mo$Rep, sep=""))
		rownames(env.mo)<-titles.env

#	distmat2<-cbind(env.mo$North_m,env.mo$East_m)
#	#distmat2<-cbind(sample_data(GP1)$edgedis_m,sample_data(GP1)$y.dist)
#	distmat2<-as.matrix(distmat2)
#	distmat2
#	pcnm2 <- pcnm(dist(distmat2))  #using easting/ northing

#env.mo<-cbind(env.mo,pcnm2$vectors)	#launch code at very end
env.mo<-sample_data(env.mo)


		##Latin names##
	datamat<-read.csv("Stat-tests_2013-14 soil unnormalised data set with biological functionsCopy.csv",header=T)
	ALLspecies.dat<-aggregate(datamat[,8:96], list(datamat$Latin.name),sum)
 			dat.2<-data.frame(cbind(env.mat[,2:5],env.mat$Loc))
			colnames(dat.2)<-c("SiteShort","Rep","Type","Type_MO","Loc")
			rownames(dat.2)<-(env.mat[,1])
		tallspecies.dat<-t(ALLspecies.dat[,2:90])
		colnames(tallspecies.dat)<-ALLspecies.dat[,1]
			rownames(tallspecies.dat)
		allspecies.test<-cbind(tallspecies.dat,dat.2)
			rownames(allspecies.test)
			colnames(allspecies.test)
			#For all soil layers##
			tspecies1314<-subset(allspecies.test, allspecies.test$SiteShort=="M2"|allspecies.test$SiteShort=="M7"|allspecies.test$SiteShort=="TR1P"|allspecies.test$SiteShort=="TR1B")
			length(tspecies1314)
			species1314.mo<-data.frame(t(tspecies1314[,1:347]))
			titles1314.mo<-c(paste(tspecies1314$SiteShort,tspecies1314$Type,tspecies1314$Rep, sep=""))
			colnames(species1314.mo)<-titles1314.mo
		     species1314.mo<-species1314.mo[,c(match(rownames(env.mo), colnames(species1314.mo)))] 
			taxa.rel<-otu_table(species1314.mo, taxa_are_rows=TRUE, errorIfNULL=TRUE)

		#FOR Presence/absence#
		species1314.pres<-species1314.mo
		species1314.pres[species1314.pres>0] <-1 
		taxa.pres<-otu_table(species1314.pres, taxa_are_rows=TRUE, errorIfNULL=TRUE)

###
phydat<-merge_phyloseq(taxa.rel,env.mo)
	colnames(sample_data(phydat))
	str(otu_table(phydat))
phydat.pres<-merge_phyloseq(taxa.pres,env.mo)


		##Trophic status##
	ALLspecies.dat<-aggregate(datamat[,8:96], list(datamat$Trophic_Status),sum)
 			dat.2<-data.frame(cbind(env.mat[,2:5],env.mat$Loc))
			colnames(dat.2)<-c("SiteShort","Rep","Type","Type_MO","Loc")
			rownames(dat.2)<-(env.mat[,1])
		tallspecies.dat<-t(ALLspecies.dat[,2:90])
		colnames(tallspecies.dat)<-ALLspecies.dat[,1]
			rownames(tallspecies.dat)
		allspecies.test<-cbind(tallspecies.dat,dat.2)
			rownames(allspecies.test)
			colnames(allspecies.test)
			#For all soil layers##
			tspecies1314<-subset(allspecies.test, allspecies.test$SiteShort=="M2"|allspecies.test$SiteShort=="M7"|allspecies.test$SiteShort=="TR1P"|allspecies.test$SiteShort=="TR1B")
			length(tspecies1314)
			species1314.mo<-data.frame(t(tspecies1314[,1:15]))
			titles1314.mo<-c(paste(tspecies1314$SiteShort,tspecies1314$Type,tspecies1314$Rep, sep=""))
			colnames(species1314.mo)<-titles1314.mo
		     species1314.mo<-species1314.mo[,c(match(rownames(env.mo), colnames(species1314.mo)))] 
			#rownames(species1314.mo)<-c("Animal.p","Myco.p","Plant.p","Sap.","Sap.brown.rot","Sap.fac.yeast","Sap.p","Sap.white.rot","Sap.yeast","AM","EM","EM.white.rot","Ericoid","Lichen","Unknown")
			rownames(species1314.mo)<-c("AP","MP","PP","ST","SBR","SFY","SP","SWR","SY","AM","EM","EMWR","ER","L","UK")
		troph.rel<-otu_table(species1314.mo, taxa_are_rows=TRUE, errorIfNULL=TRUE)
			phydat.troph<-merge_phyloseq(troph.rel,env.mo)


	##Trophic status: SUM of SPECIES presence / absence##
	ALLspecies.dat<-aggregate(datamat[,8:96], list(datamat$Latin.name,datamat$Trophic_Status),sum)
 		dat.2<-data.frame(cbind(env.mat[,2:5],env.mat$Loc))
		colnames(dat.2)<-c("SiteShort","Rep","Type","Type_MO","Loc")
		rownames(dat.2)<-(env.mat[,1])
		tallspecies.dat<-t(ALLspecies.dat[,3:91])
		colnames(tallspecies.dat)<-ALLspecies.dat[,1]
		rownames(tallspecies.dat)
		allspecies.test<-cbind(tallspecies.dat,dat.2)
			rownames(allspecies.test)
			colnames(allspecies.test)
			#For all soil layers##
			tspecies1314<-subset(allspecies.test, allspecies.test$SiteShort=="M2"|allspecies.test$SiteShort=="M7"|allspecies.test$SiteShort=="TR1P"|allspecies.test$SiteShort=="TR1B")
			length(tspecies1314)
			species1314.mo<-data.frame(t(tspecies1314[,1:347]))
			titles1314.mo<-c(paste(tspecies1314$SiteShort,tspecies1314$Type,tspecies1314$Rep, sep=""))
			colnames(species1314.mo)<-titles1314.mo
		     species1314.mo<-species1314.mo[,c(match(rownames(env.mo), colnames(species1314.mo)))] 
			#FOR Presence/absence#
			species1314.mo[species1314.mo>0] <-1
				#Trophic stat summed by species#
				ALLspecies.dat<-ALLspecies.dat[c(match(rownames(species1314.mo), ALLspecies.dat$Group.1)),]
				troph.stat<-ALLspecies.dat$Group.2
				species1314.foragg<-cbind(species1314.mo,troph.stat)
				trophbysumsp<-aggregate(species1314.foragg[,1:43], list(species1314.foragg$troph.stat),sum)
				trophbysumsp$Group.1
				species1314.mo<-data.frame(trophbysumsp[,2:44])
				#rownames(species1314.mo)<-c("Animal.p","Myco.p","Plant.p","Sap.","Sap.brown.rot","Sap.fac.yeast","Sap.p","Sap.white.rot","Sap.yeast","AM","EM","EM.white.rot","Ericoid","Lichen","Unknown")
				rownames(species1314.mo)<-c("AP","MP","PP","ST","SBR","SFY","SP","SWR","SY","AM","EM","EMWR","ER","L","UK")
		troph.pres<-otu_table(species1314.mo, taxa_are_rows=TRUE, errorIfNULL=TRUE)
			phydat.trophpres<-merge_phyloseq(troph.pres,env.mo)

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

#######################################################################
		###Fungal taxa: rel.abund.###
		min10 <- prune_samples(sample_data(GP1)$Type == "Min10", GP1)
		min10 <- prune_taxa(taxa_sums(min10 )>0, min10 )

		min35 <- prune_samples(sample_data(GP1)$Type == "Min35", GP1)
		min35 <- prune_taxa(taxa_sums(min35 )>0, min35 )

		org <- prune_samples(sample_data(GP1)$Type_MO == "org", GP1)
		org <- prune_taxa(taxa_sums(org )>0, org )

		roots <- prune_samples(sample_data(GP1)$Type_MO == "root", GP1)
		roots <- prune_taxa(taxa_sums(roots )>0, roots )
	
		bothrm <- prune_samples(sample_data(GP1)$Type == "Root"|sample_data(GP1)$Type == "Min10", GP1)
		bothrm <- prune_taxa(taxa_sums(bothrm )>0, bothrm )

			###Fungal taxa: pres.abs###
			min10.pres <- prune_samples(sample_data(GP1.pres)$Type == "Min10", GP1.pres)
			min10.pres <- prune_taxa(taxa_sums(min10.pres )>0, min10.pres )

			min35.pres <- prune_samples(sample_data(GP1.pres)$Type == "Min35", GP1.pres)
			min35.pres <- prune_taxa(taxa_sums(min35.pres )>0, min35.pres )
	
			org.pres <- prune_samples(sample_data(GP1.pres)$Type_MO == "org", GP1.pres)
			org.pres <- prune_taxa(taxa_sums(org.pres )>0, org.pres )
	
			roots.pres <- prune_samples(sample_data(GP1.pres)$Type_MO == "root", GP1.pres)
			roots.pres <- prune_taxa(taxa_sums(roots.pres )>0, roots.pres )

			bothrm.pres <- prune_samples(sample_data(GP1.pres)$Type == "Root"|sample_data(GP1.pres)$Type == "Min10", GP1.pres)
			bothrm.pres <- prune_taxa(taxa_sums(bothrm.pres)>0, bothrm.pres)

		###Trophic function: rel.abund.###
		min10.troph <- prune_samples(sample_data(GP1.troph)$Type == "Min10", GP1.troph)
		min10.troph <- prune_taxa(taxa_sums(min10.troph )>0, min10.troph )

		min35.troph <- prune_samples(sample_data(GP1.troph)$Type == "Min35", GP1.troph)
		min35.troph <- prune_taxa(taxa_sums(min35.troph )>0, min35.troph )

		org.troph <- prune_samples(sample_data(GP1.troph)$Type_MO == "org", GP1.troph)
		org.troph <- prune_taxa(taxa_sums(org.troph )>0, org.troph )

		roots.troph <- prune_samples(sample_data(GP1.troph)$Type_MO == "root", GP1.troph)
		roots.troph <- prune_taxa(taxa_sums(roots.troph )>0, roots.troph )

		bothrm.troph <- prune_samples(sample_data(GP1.troph)$Type == "Root"|sample_data(GP1.troph)$Type == "Min10", GP1.troph)
		bothrm.troph <- prune_taxa(taxa_sums(bothrm.troph )>0, bothrm.troph )

			###Trophic function: pres.abs###
			min10.trophpres<- prune_samples(sample_data(GP1.trophpres)$Type == "Min10", GP1.trophpres)
			min10.trophpres<- prune_taxa(taxa_sums(min10.trophpres)>0, min10.trophpres)

			min35.trophpres<- prune_samples(sample_data(GP1.trophpres)$Type == "Min35", GP1.trophpres)
			min35.trophpres<- prune_taxa(taxa_sums(min35.trophpres)>0, min35.trophpres)

			org.trophpres<- prune_samples(sample_data(GP1.trophpres)$Type_MO == "org", GP1.trophpres)
			org.trophpres<- prune_taxa(taxa_sums(org.trophpres)>0, org.trophpres)

			roots.trophpres<- prune_samples(sample_data(GP1.trophpres)$Type_MO == "root", GP1.trophpres)
			roots.trophpres<- prune_taxa(taxa_sums(roots.trophpres)>0, roots.trophpres)

			bothrm.trophpres <- prune_samples(sample_data(GP1.trophpres )$Type == "Root"|sample_data(GP1.trophpres )$Type == "Min10", GP1.trophpres )
			bothrm.trophpres <- prune_taxa(taxa_sums(bothrm.trophpres )>0, bothrm.trophpres )


	rec.site <- prune_samples(sample_data(GP1.pres)$ReclaimStat == "Reclamation", GP1.pres)
	rec.site <- prune_taxa(taxa_sums(rec.site )>0, rec.site )

	burn.site <- prune_samples(sample_data(GP1.pres)$ReclaimStat == "Forestry", GP1.pres)
	burn.site <- prune_taxa(taxa_sums(burn.site )>0, burn.site )

Notrec <- prune_samples(sample_data(GP1)$ReclaimStat != "Reclamation", GP1)
		rec2 <- prune_samples(sample_data(GP1)$ReclaimStat == "Reclamation", GP1)
		rec3 <- prune_taxa(taxa_sums(Notrec )==0, rec2 )
Notburn <- prune_samples(sample_data(GP1)$ReclaimStat != "Forestry", GP1)
		burn2 <- prune_samples(sample_data(GP1)$ReclaimStat == "Forestry", GP1)
		burn3 <- prune_taxa(taxa_sums(Notburn )==0, burn2 )

###################################################
####Unique taxa####################################
###################################################
Notallmin <- prune_samples(sample_data(GP1)$Type == "Root"|sample_data(GP1)$Type == "Org", GP1)
		allmin2 <- prune_samples(sample_data(GP1)$Type == "Min10"|sample_data(GP1)$Type == "Min35", GP1)
		allmin3 <- prune_taxa(taxa_sums(Notallmin )==0, allmin2  )

levels(datamat$Trophic_Status)
datamat$Latin.name[datamat$Trophic_Status=="Biotroph/Animal parasite"]
datamat$Latin.name[datamat$Trophic_Status=="Biotroph/Mycoparasite"]
datamat$Latin.name[datamat$Trophic_Status=="Saprotroph/Yeast"]
datamat$Latin.name[datamat$Trophic_Status=="Symbiotroph/Arbuscular mycorrhizal"]
datamat$Latin.name[datamat$Trophic_Status=="Symbiotroph/Ericoid"]

################################################
###########Ubiquitous taxa across all samples###
	  GPall.pres <- prune_taxa(taxa_sums(GP1.pres )>=43, GP2 )
	  paste(rownames(otu_table(GPall.pres)),rowSums(otu_table(GPall.pres)),tax_table(GPall.pres)[,2])

#####################################################################
########Diversity INDEXES############################################
#####################################################################
sample_data(GP2)
	#good method-fill#	could make color= darker version of fill color
p<-plot_richness(GP2, x = "ReclaimStat",measures=c("Observed", "Chao1","ACE","Fisher", "Shannon", "InvSimpson")) + geom_boxplot()
	print(p+xlab("Reclamation Status")+
		aes(fill=Type)+
		scale_fill_discrete(name="Soil fractions",labels=c("Coarse soils","Fine soils","Forest floors","Roots"))
		)


####################################
####################
##Network analyses##
####################
symbvec<- c("Min10"=16,"Min35"=18,"Org"=15,"Root"=17)
colvec <- c("Forestry" = "grey20","Reclamation"= "grey80")

sg <- make_network(GP1 , "samples", "bray", 0.71,laymeth="fruchterman.reingold")
	p<-plot_network(sg, GP1 , "samples",label="NA",color="ReclaimStat",line_weight = 1,hjust=-0.5, point_size=5, line_alpha =0.5,shape="Type")	
	p +theme(legend.position="right", legend.title = element_text(colour="black", size=25,face="bold"), legend.text = element_text(colour="black", size=23) )+geom_text(label=p$data$SiteFull,size=4,vjust=1.7)+ scale_colour_manual(values = colvec, name = "Stand type",label=c("Forestry","Reclamation"))+scale_shape_manual(values = symbvec, name = "Soil fraction", label=c("Coarse soil","Fine soil","Forest floor","Roots"))

	#original# size=3
 
#########################################################################
##BEST HEATMAPs##########################################################
#########################################################################
	p<-plot_heatmap(GPg , "NMDS", "bray", "SiteFull")
		p$scales$scales[[2]]$name <- "Counts per million reads" #"Normalized counts" #
		p$scales$scales[[1]]$name <- "Soil samples (Pyrosequencing)"
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
		p$scales$scales[[1]]$name <- "Soil samples (Pyrosequencing)"
			print(p+ylab("Functional groups"))


#######################################################################
####Constrained ordinations############################################
#######################################################################

##TAXA:rel. abund.
mod1<-cca(t(otu_table(min10))~sample_data(min10)$ReclaimStat+Condition(sample_data(min10)$min_NO3),scale=TRUE)
mod<-cca(t(otu_table(min35))~sample_data(min35)$ReclaimStat+Condition(sample_data(min35)$min_NO3),scale=TRUE)
mod<-cca(t(otu_table(org))~sample_data(org)$ReclaimStat+Condition(sample_data(org)$RegenYr),scale=TRUE)
mod1<-cca(t(otu_table(roots))~sample_data(roots)$ReclaimStat,scale=TRUE)
mod<-cca(t(otu_table(bothrm))~sample_data(bothrm)$ReclaimStat+Condition(sample_data(bothrm)$RegenYr+sample_data(bothrm)$min_NO3+sample_data(bothrm)$pcnm2_jp+sample_data(bothrm)$Type),scale=TRUE)
mod<-cca(t(otu_table(GP1))~sample_data(GP1)$ReclaimStat+Condition(sample_data(GP1)$RegenYr+sample_data(GP1)$min_NO3+sample_data(GP1)$min_NH4+sample_data(GP1)$pcnm2_jp+sample_data(GP1)$Type),scale=TRUE)

##TAXA:pres/abs
mod2<-cca(t(otu_table(min10.pres))~sample_data(min10)$ReclaimStat,scale=TRUE)
mod<-cca(t(otu_table(min35.pres))~sample_data(min35)$ReclaimStat+Condition(sample_data(min35)$min_NO3),scale=TRUE)
mod<-cca(t(otu_table(org.pres))~sample_data(org)$ReclaimStat+Condition(sample_data(org)$RegenYr+sample_data(org)$pcnm2_jp),scale=TRUE)
mod<-cca(t(otu_table(roots.pres))~sample_data(roots)$ReclaimStat+Condition(sample_data(roots)$min_NO3),scale=TRUE)
mod<-cca(t(otu_table(bothrm.pres))~sample_data(bothrm)$ReclaimStat+Condition(sample_data(bothrm)$min_NH4+sample_data(bothrm)$min_NO3+sample_data(bothrm)$pcnm2_jp+sample_data(bothrm)$Type),scale=TRUE)
mod<-cca(t(otu_table(GP1.pres))~sample_data(GP1)$ReclaimStat+Condition(sample_data(GP1)$RegenYr+sample_data(GP1)$min_NO3+sample_data(GP1)$min_NH4+sample_data(GP1)$min_totalN+sample_data(GP1)$pcnm2_jp+sample_data(GP1)$Type),scale=TRUE)

##FUNCTION: rel. abund.
mod<-rda(t(otu_table(min10.troph))~sample_data(min10)$ReclaimStat+Condition(sample_data(min10)$min_NO3),scale=TRUE)
mod<-rda(t(otu_table(min35.troph))~sample_data(min35)$ReclaimStat+Condition(sample_data(min35)$min_NH4+sample_data(min35)$pcnm2_jp),scale=TRUE)
mod<-rda(t(otu_table(org.troph))~sample_data(org)$ReclaimStat+Condition(sample_data(org)$pcnm2_jp),scale=TRUE)
mod<-rda(t(otu_table(roots.troph))~sample_data(roots)$ReclaimStat+Condition(sample_data(roots)$RegenYr),scale=TRUE)
mod<-rda(t(otu_table(bothrm.troph))~sample_data(bothrm)$ReclaimStat+Condition(sample_data(bothrm)$RegenYr+sample_data(bothrm)$pcnm2_jp+sample_data(bothrm)$Type),scale=TRUE)
mod<-rda(t(otu_table(GP1.troph))~sample_data(GP1)$ReclaimStat+Condition(sample_data(GP1)$min_NO3+sample_data(GP1)$pcnm2_jp+sample_data(GP1)$Type),scale=TRUE)

##FUNCTION: pres/abs
mod<-rda(t(otu_table(min10.trophpres))~sample_data(min10)$ReclaimStat,scale=TRUE)
mod<-rda(t(otu_table(min35.trophpres))~sample_data(min35)$ReclaimStat+Condition(sample_data(min35)$min_NH4),scale=TRUE)
mod<-rda(t(otu_table(org.trophpres))~sample_data(org)$ReclaimStat+Condition(sample_data(org)$min_NO3),scale=TRUE)
mod<-rda(t(otu_table(roots.trophpres))~sample_data(roots)$ReclaimStat,scale=TRUE)
mod<-rda(t(otu_table(bothrm.trophpres))~sample_data(bothrm)$ReclaimStat+Condition(sample_data(bothrm)$min_NH4+sample_data(bothrm)$pcnm2_jp+sample_data(bothrm)$Type),scale=TRUE)
mod<-rda(t(otu_table(GP1.trophpres))~sample_data(GP1)$ReclaimStat+Condition(sample_data(GP1)$min_NH4+sample_data(GP1)$pcnm2_jp+sample_data(GP1)$Type),scale=TRUE)

	permutest(mod,permu=999)
	anova(mod,by="term",perm.max=999)
	summary(mod)


####################################################################
####################################################################
########################################
###plotting ordinations#################
########################################
mod1<-cca(t(otu_table(bothrm))~sample_data(bothrm)$ReclaimStat+Condition(sample_data(bothrm)$RegenYr+sample_data(bothrm)$min_NO3+sample_data(bothrm)$pcnm2_jp+sample_data(bothrm)$Type),scale=TRUE)
mod2<-cca(t(otu_table(bothrm.pres))~sample_data(bothrm)$ReclaimStat+Condition(sample_data(bothrm)$min_NH4+sample_data(bothrm)$min_NO3+sample_data(bothrm)$pcnm2_jp+sample_data(bothrm)$Type),scale=TRUE)
mod3<-rda(t(otu_table(bothrm.troph))~sample_data(bothrm)$ReclaimStat+Condition(sample_data(bothrm)$RegenYr+sample_data(bothrm)$pcnm2_jp+sample_data(bothrm)$Type),scale=TRUE)
mod4<-rda(t(otu_table(bothrm.trophpres))~sample_data(bothrm)$ReclaimStat+Condition(sample_data(bothrm)$min_NH4+sample_data(bothrm)$pcnm2_jp+sample_data(bothrm)$Type),scale=TRUE)


scl <- 0 ## scl=0 (unscaled raw scores)#scl=2 (species scores scaled by eigenvalue) #scaling = 3 (both scecies and site scores are scaled symmetrically by square root of eigenvalue)
colvec <- c("black")
#symbvec<- c("Min10"=21,"Org"=22,"Min35"=23,"Root"=24)#use symbols for soil_layer (Min=circle vs Org=square vs Seedling= diamond vs Soil=triangle)
symbvec<- c("Min10"=21,"Root"=24)#use symbols for soil_layer (Min=circle vs Org=square vs Seedling= diamond vs Soil=triangle)
#cexvec<-c("none"=2.22,"Reclaim"=1.8) # reclaim versus natural stand 
#cexvec<-1
colorSp <- c("Forestry" = "black", "Reclamation" = "grey")
#colorSp <- c("Jack Pine" = "black", "Sib. Larch" = "grey","Wt. Spruce" = "white")

par(mfrow=c(2,2))	
#par(oma=c(0,0,0,0))
par(mar=c(5.1,4.1,0.5,2.1))

mintotal.mod<-mod2

minfit<-scores(mintotal.mod, display = "sites", shrink = FALSE, choices=c(1)) 
		wss <- (nrow(minfit)-1)*sum(apply(minfit,2,var))
		for (i in 2:10) wss[i] <- sum(kmeans(minfit, 
	  		 centers=i)$withinss)
		#plot(1:10, wss, type="b", xlab="Number of Clusters",
	  	#	ylab="Within groups sum of squares")
	mintotal.kclust<-kmeans(minfit, centers=2, iter.max = 10, nstart = 25)
		mintotal.kclust$cluster
	
	plot(mintotal.mod, type = "n", scaling = scl ,xlab="",ylab="")
		    #text(mintotal.mod,choices=c(1,2), display = "spec", cex = 1, scaling = scl, pch=3, col="red")
		    points(mintotal.mod,choices=c(1,2), display = "spec", cex = 0.8, scaling = scl, pch=3, col="red")
		  	#text(mintotal.mod,choices=c(1,2), display = "cn", cex=1, col="darkred")

	ordihull(mintotal.mod, mintotal.kclust$cluster, draw = "lines",	##plotting clusters  
     		    col = "black",scaling = scl )

	    #bothrm#color=reclaim, size = regen yr, pch = stand stg
	    with(sample_data(bothrm), points(mintotal.mod,choices=c(1,2), display = "sites", col = "black",
	      scaling = scl , pch = symbvec[Type], bg = colorSp[ReclaimStat],cex=2))
	
#	legend("bottomright",pch=c(21,21,  22,21,23,24)
#	,col="black",pt.bg=c("black","grey"	,"black","black","black","black"),
#	legend=c("Forestry","Reclamation","Forest floor","Coarse soil","Fine soil","Roots"),cex=1
#	,pt.cex=c(1.9,1.9,  1.9,1.9,1.9,1.9),bty="t",bg="white")  

	legend("bottomleft",pch=c(21,21,  21,24)
	,col="black",pt.bg=c("black","grey","black","black"),
	legend=c("Forestry","Reclamation","Coarse soil","Roots"),cex=1
	,pt.cex=c(1.9,1.9,  1.9,1.9),bty="t",bg="white") 

  specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
#title(main="",xlab = paste("CA1 ",specify_decimal(mintotal.mod$CA$eig[1]/sum(mintotal.mod$CA$eig)*100,2),"%"), ylab = paste("CA2 ",specify_decimal(mintotal.mod$CA$eig[2]/sum(mintotal.mod$CA$eig)*100,2),"%"),cex.lab=1.2,line=2.5)  
title(main="",xlab = paste("CCA1 ",specify_decimal(mintotal.mod$CCA$tot.chi/mintotal.mod$tot.chi*100,2),"%"), ylab = paste("CA1 ",specify_decimal(mintotal.mod$CA$eig[1]/sum(mintotal.mod$CA$eig)*100,2),"%"),cex.lab=1.2,line=2.5)  
#title(main=paste(specify_decimal(mintotal.mod$CCA$tot.chi/mintotal.mod$tot.chi*100,2),"% of total inertia"),xlab = paste("CCA1 ",specify_decimal(mintotal.mod$CCA$eig[1]/sum(mintotal.mod$CCA$eig)*100,2),"%"), ylab = paste("CCA2 ",specify_decimal(mintotal.mod$CCA$eig[2]/sum(mintotal.mod$CCA$eig)*100,2),"%"),cex.lab=1.2,line=2.5)  
#for RDA
#title(main="",xlab = paste("PC1 ",specify_decimal(mintotal.mod$CA$eig[1]/sum(mintotal.mod$CA$eig)*100,2),"%"), ylab = paste("PC2 ",specify_decimal(mintotal.mod$CA$eig[2]/sum(mintotal.mod$CA$eig)*100,2),"%"),cex.lab=1.2,line=2.5)  
title(main="",xlab = paste("RDA1 ",specify_decimal(mintotal.mod$CCA$tot.chi/mintotal.mod$tot.chi*100,2),"%"), ylab = paste("PC1 ",specify_decimal(mintotal.mod$CA$eig[1]/sum(mintotal.mod$CA$eig)*100,2),"%"),cex.lab=1.2,line=2.5)  
#title(main=paste(specify_decimal(mintotal.mod$CCA$tot.chi/mintotal.mod$tot.chi*100,2),"% of total inertia"),xlab = paste("RDA1 ",specify_decimal(mintotal.mod$CCA$eig[1]/sum(mintotal.mod$CCA$eig)*100,2),"%"), ylab = paste("RDA2 ",specify_decimal(mintotal.mod$CCA$eig[2]/sum(mintotal.mod$CCA$eig)*100,2),"%"),cex.lab=1.2,line=2.5)  

	#reclaimStat +cond.(non-corr. factors)	#scl=0 
	 text(2.3,4,"A",cex=4)	
	 text(-5.4,4.3,"B",cex=4)	
	 text(0.47,-0.64,"C",cex=4)	
	 text(-0.49,-0.36,"D",cex=4)	

		#reclaimStat [NO conditions]	#scl=0 
		 text(2.9,-4.4,"A",cex=4)	
		 text(-3,-3.4,"B",cex=4)	
		 text(0.44,-0.48,"C",cex=4)	
		 text(-0.5,-0.43,"D",cex=4)	



########################################################################
###upregulated taxa###
###Up/down reg by species###
mod.tax<-cca(t(otu_table(bothrm))~sample_data(bothrm)$ReclaimStat+Condition(sample_data(bothrm)$RegenYr+sample_data(bothrm)$min_NO3+sample_data(bothrm)$pcnm2_jp+sample_data(bothrm)$Type),scale=TRUE)
mod.tax1<-cca(t(otu_table(bothrm))~sample_data(bothrm)$ReclaimStat+Condition(sample_data(bothrm)$Type),scale=TRUE)
mod.pres<-cca(t(otu_table(bothrm.pres))~sample_data(bothrm)$ReclaimStat+Condition(sample_data(bothrm)$min_NH4+sample_data(bothrm)$min_NO3+sample_data(bothrm)$pcnm2_jp+sample_data(bothrm)$Type),scale=TRUE)
mod.pres1<-cca(t(otu_table(bothrm.pres))~sample_data(bothrm)$ReclaimStat+Condition(sample_data(bothrm)$Type),scale=TRUE)

mod.tax<-mod1
mod.pres<-mod2
 
 	tax.reclaim<-attributes(mod.tax$CCA$v[mod.tax$CCA$v[,1]< -0.5,])	#need space b/w < and -
	  GP1.tax.reclaim<- prune_taxa( unlist(tax.reclaim) , GP2)
 		tax1.reclaim<-attributes(mod.tax1$CCA$v[mod.tax1$CCA$v[,1]< -0.5,])	#need space b/w < and -
		  GP1.tax1.reclaim<- prune_taxa( unlist(tax1.reclaim) , GP2)
 	  match.tax.reclaim<- prune_taxa( unlist(tax1.reclaim) , GP1.tax.reclaim)

	tax.forestry<-attributes(mod.tax$CCA$v[mod.tax$CCA$v[,1]> 0.5,])	#need space b/w < and -
	  GP1.tax.forestry<- prune_taxa( unlist(tax.forestry) , GP2)
		tax1.forestry<-attributes(mod.tax1$CCA$v[mod.tax1$CCA$v[,1]> 0.5,])	#need space b/w < and -
		  GP1.tax1.forestry<- prune_taxa( unlist(tax1.forestry) , GP2)
 	  match.tax.forestry<- prune_taxa( unlist(tax1.forestry) , GP1.tax.forestry)

		pres.reclaim<-attributes(mod.pres$CCA$v[mod.pres$CCA$v[,1]< -0.5,])	#need space b/w < and -
		  GP1.pres.reclaim<- prune_taxa( unlist(pres.reclaim) , GP2)
			pres1.reclaim<-attributes(mod.pres1$CCA$v[mod.pres1$CCA$v[,1]< -0.5,])	#need space b/w < and -
			  GP1.pres1.reclaim<- prune_taxa( unlist(pres1.reclaim) , GP2)
	 	  match.pres.reclaim<- prune_taxa( unlist(pres1.reclaim) , GP1.pres.reclaim)

		pres.forestry<-attributes(mod.pres$CCA$v[mod.pres$CCA$v[,1]> 0.5,])	#need space b/w < and -
		  GP1.pres.forestry<- prune_taxa( unlist(pres.forestry) , GP2)
			pres1.forestry<-attributes(mod.pres1$CCA$v[mod.pres1$CCA$v[,1]> 0.5,])	#need space b/w < and -
			  GP1.pres1.forestry<- prune_taxa( unlist(pres1.forestry) , GP2)
	 	  match.pres.forestry<- prune_taxa( unlist(pres1.forestry) , GP1.pres.forestry)

	 	  #match.both.reclaim <- prune_taxa( rownames(otu_table( match.pres.reclaim )) ,  match.tax.reclaim )
	 	  #match.both.forestry <- prune_taxa( rownames(otu_table(match.pres.forestry )) ,  match.tax.forestry )

	 	  match.both.reclaim <- prune_taxa( rownames(otu_table( GP1.pres.reclaim )) ,  GP1.tax.reclaim )
	 	  match.both.forestry <- prune_taxa( rownames(otu_table(GP1.pres.forestry )) ,  GP1.tax.forestry )

	datamat<-read.csv("Stat-tests_2013-14 soil unnormalised data set with biological functionsCopy.csv",header=T)
datamat[match(rownames(otu_table(match.both.reclaim )),datamat$Latin.name),3:4]
datamat[match(rownames(otu_table(match.both.forestry )),datamat$Latin.name),3:4]

rowSums(otu_table(match.both.reclaim ))
rowSums(otu_table(match.both.forestry ))


length(datamat$Latin.name[datamat$Trophic_Status=="Symbiotroph/Arbuscular mycorrhizal"]
)


#################
 reclaimPyro <- prune_taxa( rownames(otu_table( reclaim_1314 )) ,  match.both.reclaim )
 forestryPyro <- prune_taxa( rownames(otu_table(forestry_1314 )) ,  match.both.forestry )

