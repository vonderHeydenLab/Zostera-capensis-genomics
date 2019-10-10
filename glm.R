#read datasets
env_status <- read.csv(file="env_status.csv", header = TRUE,sep=";", dec=".")
AR <- read.csv(file="AR_env_loci_edit.csv", header = TRUE,sep=";", dec=".")
He <- read.csv(file="he_env_loci_edit1.csv", header = TRUE,sep=";", dec=".")
Nu_div <- read.csv(file="nucleotide_div_env_edit.csv", header = TRUE,sep=";", dec=".")
area<- read.csv(file="diversity_indices_status.csv", header = TRUE,sep=";", dec=".")

#split the genetic diversity data into 2 clusters. Run on the complete dataset as well as
#2 clusters seperately to see if they respond differently to disturbances

#keep row 1 to row 5 for cluster 1
AR_1<- rbind(AR[1:5,])
He_1<- rbind(He[1:5,])

#keep row 6 to 8 for cluster 2
AR_2<- rbind(AR[6:8,])
He_2<- rbind(He[6:8,])
#Nu-div is in a different format- edit this by hand
Nu_div_1 <- read.csv(file="nucleotide_div_env_edit_cluster1.csv", header = TRUE,sep=";", dec=".")
Nu_div_2 <- read.csv(file="nucleotide_div_env_edit_cluster2.csv", header = TRUE,sep=";", dec=".")

cor.test(area$AR,area$submerged_macrophyte_area,type="spearman")
cor.test(area$He,area$submerged_macrophyte_area,type="spearman")
cor.test(area$nucleotide_diversity,area$submerged_macrophyte_area,type="spearman")

#use reshape 'melt' function to 'stack' the data 
library(reshape2)

AR_melt <- melt(AR, id =c('site','change_in_flow',	'pollution',	'habitat_loss',	
                          'mining',	'artificial_breaching',	'fishing_effort',	'bait_collection',	'ecological_category', 'submerged_macrophyte_area'))
He_melt <- melt(He, id =c('site','change_in_flow',	'pollution',	'habitat_loss',	
                          'mining',	'artificial_breaching',	'fishing_effort',	'bait_collection',	'ecological_category','submerged_macrophyte_area'))

just_He<-cbind(He_melt$value)
just_AR<-cbind(AR_melt$value)
just_Nu_div<-cbind(Nu_div$nucleotide_diversity)

# Use the glm because you can see the effect of each level of each facor
#####GLM Allelic richness
AR_glm<- glm(just_AR~ pollution + habitat_loss + mining + fishing_effort + change_in_flow , data = AR_melt, family=gaussian)
#removed ecological_category, 'step' still results in the same final model
#a_glm<-anova(AR_glm)
#anova(AR_glm, test="LRT")
summary(AR_glm)
#alias(AR_glm, complete=TRUE)

step(AR_glm)
AR_glm_reduced<-glm(formula = just_AR ~ habitat_loss + fishing_effort, family = gaussian, 
                  data = AR_melt)
summary(AR_glm_reduced)

#######################################################
#significant:fishing_effort, negative relationship; non-sig:habitat loss
#######################################################

#####GLM Expected heterozygosity
He_glm<- glm(just_He~ pollution + habitat_loss + mining + 
               fishing_effort + change_in_flow , data = He_melt, family=gaussian)
#a_He_glm<-anova(He_glm)
summary(He_glm)
step(He_glm)

He_glm_reduced<-glm(just_He ~ mining + fishing_effort + change_in_flow, 
                    family = gaussian, data = He_melt)
summary(He_glm_reduced)

summary(aov(just_He~ pollution + habitat_loss + mining + fishing_effort + change_in_flow , data = He_melt))

#################################################################################
#mining, negative; fishing effort, negative; change in flow, positive(non-sig)
#################################################################################

#####GLM Nucleotide diversity
Nu_div_glm<- glm(just_Nu_div~ pollution + habitat_loss + mining + 
                   fishing_effort + change_in_flow + ecological_category+submerged_macrophyte_area, data = Nu_div, family=gaussian)
Nu_div_glm<- glm(just_Nu_div~ pollution + habitat_loss + mining + 
                   fishing_effort + change_in_flow , data = Nu_div, family=gaussian)
#a_Nu_div_glm<-anova(Nu_div_glm)
summary(Nu_div_glm)
step(Nu_div_glm)
Nu_div_glm_reduced<- glm(just_Nu_div ~ habitat_loss + mining + fishing_effort, 
                         family = gaussian, data = Nu_div)
summary(Nu_div_glm_reduced)

summary(aov(just_Nu_div~ pollution + habitat_loss + mining + fishing_effort + change_in_flow , data = Nu_div))

#############################################################################
#habitat_loss,negative;mining,
#############################################################################

###################################################
#Post hoc Tukey HSD 
###################################################

library(vegan)

He_melt$site <- as.factor(He_melt$site)
He.aov <- aov(He_melt$value ~ He_melt$site, data = He_melt)
summary(He.aov)
He.tukey<-TukeyHSD(He.aov)

# Group the treatments that are not different from each other together.
library(multcompView)
generate_label_df <- function(He.tukey, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- He.tukey[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}

# Apply the function on my dataset
LABELS=generate_label_df(He.tukey , "He_melt$site")


# A panel of colors to draw each group with the same color :
my_colors=c( rgb(143,199,74,maxColorValue = 255),rgb(242,104,34,maxColorValue = 255), rgb(111,145,202,maxColorValue = 255),rgb(254,188,18,maxColorValue = 255) , rgb(74,132,54,maxColorValue = 255),rgb(236,33,39,maxColorValue = 255),rgb(165,103,40,maxColorValue = 255))

# Draw the basic boxplot
a=boxplot(He_melt$value ~ He_melt$site , ylim=c(min(He_melt$value) , 1.1*max(He_melt$value)) , col=my_colors[as.numeric(LABELS[,1])] , ylab="He" , xlab="site", main="He per locus across sites", outline=FALSE)

# I want to write the letter over each box. Over is how high I want to write it.
over=0.1*max( a$stats[nrow(a$stats),] )

#Add the labels
text( c(1:nlevels(He_melt$site)) , a$stats[nrow(a$stats),]+over , LABELS[,1]  , col=my_colors[as.numeric(LABELS[,1])] ,pos=3, offset = 0.5)

AR_melt$site <- as.factor(AR_melt$site)
AR.aov <- aov(AR_melt$value ~ AR_melt$site, data = AR_melt)
summary(AR.aov)
TukeyHSD(AR.aov)

Nu_div$site <- as.factor(Nu_div$site)
Nu_div.aov <- aov(Nu_div$nucleotide_diversity ~ Nu_div$site, data = Nu_div)
summary(Nu_div.aov)
Nu.tukey<-TukeyHSD(Nu_div.aov)

####plotting
generate_label_df.nu.div <- function(Nu.tukey, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- Nu.tukey[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}

# Apply the function on my dataset
LABELS=generate_label_df(Nu.tukey , "Nu_div$site")


# A panel of colors to draw each group with the same color :
my_colors=c( rgb(143,199,74,maxColorValue = 255),rgb(242,104,34,maxColorValue = 255), rgb(111,145,202,maxColorValue = 255),rgb(254,188,18,maxColorValue = 255) , rgb(74,132,54,maxColorValue = 255),rgb(236,33,39,maxColorValue = 255),rgb(165,103,40,maxColorValue = 255))

# Draw the basic boxplot
a=boxplot(Nu_div$nucleotide_diversity ~ Nu_div$site , ylim=c(min(Nu_div$nucleotide_diversity) , 1.1*max(Nu_div$nucleotide_diversity)) , col=my_colors[as.numeric(LABELS[,1])] , ylab="nucleotide diversity", xlab="site" , main="Nucleotide diversity per locus across sites", outline=FALSE)

# I want to write the letter over each box. Over is how high I want to write it.
over=0.1*max( a$stats[nrow(a$stats),] )

#Add the labels
text( c(1:nlevels(Nu_div$site)) , a$stats[nrow(a$stats),]+over , LABELS[,1]  , col=my_colors[as.numeric(LABELS[,1])] ,pos=3, offset = 0.5)

