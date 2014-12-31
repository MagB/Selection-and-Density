##DATA  ANALYSIS STARTS HERE
#1.  How are phenotypes & absolute seed # associated with density?
#  •  Table of all correlations, including density as a character. Actually, should probably make:
#  •	1st table for all trait correlations (z) 
#  •	2nd table for corr between each trait and density of each scale. 
# I need teh Hmisc package to obtain pvalues for the correlations.

#Here is a helpful quick R tutorial that covers some basic concepts.
#library("swirl")
#swirl()

install.packages("Hmisc")
library("Hmisc")
options(digits=3)




#2011 datasets are split into site 1 (this includes the 800s plants) and site 2 and site4. Note that site 4 has different data that the other 2 sites.
s1_traits_density=read.csv(file="./output/s1_traits_density_2011.csv", header=T)
s2_traits_density=read.csv(file="./output/s2_traits_density_2011.csv", header=T)
s4_traits_density=read.csv(file="./output/s4_traits_density_2011.csv", header=T)

##20009 data
s1_traits_2009=read.csv(file="./output/s1_traits_density_2009.csv", header=T)


# data sets:s1_traits_density
head(s1_traits_density)
head(s2_traits_density)
head(s1_traits_2009)
head(s4_traits_density)

#step 1 adjust density!!!!! AT this point December 12 2014, we have not adjusted the density data for site 2 to correct for plants that were present on the site but not tagged


names(s1_traits_density)
ncol(s1_traits_density)


##Step 2 all pairwise correlations:


names(s1_traits_density)
traits=(13:27)

##Correlation among traits for s1
cor(s1_traits_density[,13:27], use="complete.obs", method="pearson")
rcorr(as.matrix(s1_traits_density[,13:27]))

##Correlation among traits and density for s1
dens_trts=c(6:11, 13:27)
cor(s1_traits_density[,dens_trts], use="complete.obs", method="pearson")
rcorr(as.matrix(s1_traits_density[,dens_trts]))

#Correlation among traits for 2009
names(s1_2009_tr_density)

trs=c(4:19,21,26)
both=c(trs, 32:37)

#Correlation among traits for s1 in 2009
cor(s1_2009_tr_density[,trs],use="complete.obs", method="pearson")
rcorr(as.matrix(s1_2009_tr_density[,trs]))

#Correlation among traits and density for s1 in 2009
cor(s1_2009_tr_density[,both], use="complete.obs", method="pearson")
rcorr(as.matrix(s1_2009_tr_density[,both]))
#prints out multiple empty lines to make some space in the console
#cat('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n')
    
    
    
    
##########REGRESSIONS:
names(s1_2009_tr_density)
s1_2009_tr_density$stem_diam_st=(s1_2009_tr_density$stem_diam-mean(s1_2009_tr_density$stem_diam, na.rm=TRUE))/sd(s1_2009_tr_density$stem_diam, na.rm=TRUE) 
s1_2009_tr_density$dens_Rad_200_st=(s1_2009_tr_density$dens_Rad_200-mean(s1_2009_tr_density$dens_Rad_200, na.rm=TRUE))/sd(s1_2009_tr_density$dens_Rad_200, na.rm=TRUE) 

s1_2009_tr_density$rel_seed=s1_2009_tr_density$seed/mean(s1_2009_tr_density$seed, na.rm=TRUE)

reg1=lm(s1_2009_tr_density$rel_seed~s1_2009_tr_density$stem_diam_st + s1_2009_tr_density$dens_Rad_10_st, data=s4_traits_density)
summary(reg1)

reg1=lm(s1_2009_tr_density$rel_seed~s1_2009_tr_density$stem_diam_st + s1_2009_tr_density$dens_Rad_200_st + s1_2009_tr_density$stem_diam_st*s1_2009_tr_density$dens_Rad_200_st, data=s4_traits_density)
summary(reg1)

library(scatterplot3d)
#THIS WILL MAKE A 3D SCATTER PLOT
s3d=scatterplot3d(s1_2009_tr_density$dens_Rad_200_st , s1_2009_tr_density$stem_diam_st, s1_2009_tr_density$rel_seed,pch=16, highlight.3d=TRUE, type="h", main="3D Scatterplot")

#THIS 2 LINES OF CODE MAKE A PLANE RUN THROUH THE 3D SCATTER PLOT
fit <- lm(s1_2009_tr_density$rel_seed ~ s1_2009_tr_density$stem_diam_st+s1_2009_tr_density$dens_Rad_200_st) 
s3d$plane3d(fit)




reg1=lm(seed~X10cm+X30cm+X200cm + Flr.num, data=s4_traits_density)
    
reg1=lm(seed~X10cm+ length.to.last.flr + Flr.num, data=s4_traits_density)
summary(reg1)
    
reg2=lm(seed~X30cm+ length.to.last.flr + Flr.num, data=s4_traits_density)
summary(reg2)
    
reg3=lm(seed~X200cm+ length.to.last.flr + Flr.num, data=s4_traits_density)
summary(reg3)
    
    
    #what about just density
reg1=lm(seed~X10cm, data=s4_traits_density, na.rm=T)
    plot(s4_traits_density$X10cm, s4_traits_density$seed)
    summary(reg1)
    
    reg2=lm(seed~X30cm, data=s4_traits_density)
    plot(s4_traits_density$X30cm, s4_traits_density$seed)
    summary(reg2)
    
    reg3=lm(seed~X200cm, data=s4_traits_density)
    plot(s4_traits_density$X200cm, s4_traits_density$seed)
    summary(reg3)
    
    #what about seed per fruit?
    s4_traits_density$seed_fruit=s4_traits_density$seed/s4_traits_density$fruit
    reg1=lm(seed[s4_traits_density$fruit_01==1 & s4_traits_density$seed<10000]~X10cm[s4_traits_density$fruit_01==1  & s4_traits_density$seed<10000], data=s4_traits_density)
    plot(s4_traits_density$X30cm[s4_traits_density$fruit_01==1 & s4_traits_density$seed<10000], s4_traits_density$seed[s4_traits_density$fruit_01==1 & s4_traits_density$seed<10000])
    summary(reg1)
    
    reg2=lm(seed_fruit~X30cm, data=s4_traits_density)
    plot(s4_traits_density$X30cm, s4_traits_density$seed_fruit)
    summary(reg2)
    
    reg3=lm(seed_fruit~X200cm, data=s4_traits_density)
    plot(s4_traits_density$X200cm, s4_traits_density$seed_fruit)
    summary(reg3)
    
    
    #what about just density and fruit number as a logistic model
    s4_traits_density$fruit[is.na(s4_traits_density$fruit)]==0
    s4_traits_density$fruit_01=0
    s4_traits_density$fruit_01[s4_traits_density$fruit >0] = 1
    s4_traits_density$fruit_01[s4_traits_density$fruit==0] =0
    
    s4_traits_density$fruit_01 = factor(s4_traits_density$fruit_01)
    mylogit = glm(fruit_01~X10cm, data=s4_traits_density, family="binomial")
    summary(mylogit)
    str(s4_traits_density$fruit_01)

    MnewData$NHets<-apply(newData[, GenoColumns], 1,GetNHet)
    
    
    
    
    reg1=lm(fruit~X10cm, data=s4_traits_density, na.rm=T)
    
    plot(s4_traits_density$X10cm, s4_traits_density$fruit)
    summary(reg1)
    
    reg2=lm(fruit~X30cm, data=s4_traits_density)
    plot(s4_traits_density$X30cm, s4_traits_density$fruit)
    summary(reg2)
    
    reg3=lm(fruit~X200cm, data=s4_traits_density)
    plot(s4_traits_density$X200cm, s4_traits_density$fruit)
    summary(reg3)
