setwd("/Users/Maggie/Documents/Lobelia cardinalis/Selection and Density 2014/")
install.packages("chron")
library(chron)
chron(s1_traits_2009$med_jule_f)
chron(s1_traits_density_2011$med_flr_date)

#open whatever datasets you want to work with
#Note that each dataset uses slightly different variable names, we can change the variable names to be consistent or just use whichever name is appropriate for each dataset.

s1_traits_density_2011=read.csv(file="./output/s1_traits_density_2011.csv", header=T)
names(s1_traits_density_2011)
s1_traits_density_2011$med_flr_date

s2_traits_density_2011=read.csv(file="./output/s2_traits_density_2011.csv", header=T)
names(s2_traits_density_2011)


s4_traits_density_2011=read.csv(file="./output/s4_traits_density_2011.csv", header=T)
names(s4_traits_density_2011)
traits_2011=c("total_seed","total_fruit", "mid_wid", "snd", "med_flr_date" ,"avg_display","flr_num", "height_last","dens_Rad_10","dens_Rad_20","dens_Rad_30","dens_Rad_50","dens_Rad_100","dens_Rad_200")


##20009 data
s1_traits_2009=read.csv(file="./output/s1_traits_density_2009.csv", header=T)

names(s1_traits_2009)
#for 2009 use these traits:
traits_2009=c("total_seed","fruit_num", "mid_wid", "snd", "med_jule_f" ,"avg_flr_num","flower_num", "ht_last_flr","dens_Rad_10","dens_Rad_20","dens_Rad_30","dens_Rad_50","dens_Rad_100","dens_Rad_200")

names(s1_traits_2009)[which(!names(s1_traits_2009) %in%  names(s1_traits_density_2011))]

#this function calculates the raw regression coefficients for whatever traits the user defines, the default is the 6 traits we typically use.
raw_regression_univariate <- function(data,trait=c("snd","mid_wid", "len_last", "flr_num","avg_display","height_last", "med_flr_date"),y="total_seed" ){

  if(!(is.data.frame(data))){stop("invalid dataset name")}  
  output_file=data.frame(traits_list=as.character(), intercept=numeric(0), estimate=numeric(0), stringsAsFactors=FALSE)
  colnames(output_file)=c("trait", "intercept","estimate")

  
  trait_columns=c(which(colnames(data) %in% trait))
  y_var_pos=which(colnames(data) %in% y)
  for (i in trait_columns) {
    
    myvars=c(y_var_pos, i)
    
    subdata=data[,myvars]
    subdata=subdata[complete.cases(subdata), ]
    print(nrow(subdata))
    #this standardizes the traits of interest
    subdata$standardized_trait=(subdata[,2]-mean(subdata[,2], na.rm=TRUE))/sd(subdata[,2], na.rm=TRUE) 
   
    #this relatives fitness
    subdata$rel_fitness=subdata[,1]/mean(subdata[,1])
    regress_raw<-as.numeric(lm(subdata$rel_fitness~   subdata$standardized_trait , data=subdata)$coef)

    newrow=c(colnames(data[i]), regress_raw)

    output_file=rbind(output_file, data.frame(traits_list=colnames(data[i]), intercept=regress_raw[1],estimate=regress_raw[2]))
  
    
  }
  #colnames(output_file)=c("trait", "intercept","estimate")
  return(output_file)
}

#the function raw_regression_univariate will also prints out the sample sizes
univariate_estimates=raw_regression_univariate(s1_traits_density_2011,c("snd","mid_wid", "len_last", "flr_num","avg_display", "med_flr_date"))
univariate_estimates=raw_regression_univariate(s1_traits_density_2011,c("snd","mid_wid", "len_last", "flr_num","avg_display", "med_flr_date", "dens_Rad_10", "dens_Rad_20", "dens_Rad_30", "dens_Rad_50", "dens_Rad_100", "dens_Rad_200"))





##FUNCTION FOR BOOTSTRAPPING THE UNIVARIATE ESTIMATES STARTS HERE
#This function calculates the bootstrap coefficients
output_boot_estimates=data.frame()
boot_uni= function(data,trait,y="total_seed"){
  set.seed(2)

  trait_column=c(which(colnames(data) %in% trait))
  y_var_pos=which(colnames(data) %in% y)

  myvars=c(y_var_pos,  trait_column)
  subdata=data[,myvars]
  subdata=subdata[complete.cases(subdata), ]
  
  output_boot_estimates=data.frame(intercept=numeric(0), estimate=numeric(0), stringsAsFactors=FALSE)
 
  
for(i in c(1:10000)){
  #this creates a randomized dataset
  subdata_boot <- subdata[sample(1:nrow(subdata), nrow(subdata), replace=TRUE),]#this samples with replacement within this data set that consists of only the trait and fitness value for the desired trait
 
  #this standardizes the traits of interest
  subdata_boot$standardized_trait=(subdata_boot[,2]-mean(subdata_boot[,2], na.rm=TRUE))/sd(subdata_boot[,2], na.rm=TRUE) 
  
  #this relatives fitness
  subdata_boot$rel_fitness=subdata_boot[,1]/mean(subdata_boot[,1])
  
  #this calculates the single univariate regression for the randomized data set
  regress_boot=as.numeric(lm(subdata_boot$rel_fitness~   subdata_boot$standardized_trait , data=subdata_boot)$coef)
  
  output_boot_estimates=rbind(output_boot_estimates, (regress_boot))

}
  output_boot_estimates$trait=trait
  colnames(output_boot_estimates)=c("intercept","estimate","trait")
  
  
  return(output_boot_estimates)
}


get_CIs = function(data){

  quants=data.frame()
  q1=quantile(data$estimate, c(0.025, 0.975))
  pval=(length(data$estimate[data$estimate<=0])/length(data$estimate))
  quants=data.frame(q1[1],q1[2])
  quants$p_val=pval  
  quants$trait=unique(data$trait)
  rownames(quants) <- 1:nrow(quants)
  
  colnames(quants)=c("2.5%","97.5%","p-val", "trait")
  return(quants)

}


univariate_boot_data=boot_uni(s1_traits_density_2011,"mid_wid")
univariate_boot_data=boot_uni(s1_traits_density_2011,"snd")
univariate_boot_data=boot_uni(s1_traits_density_2011,"len_last")
univariate_boot_data=boot_uni(s1_traits_density_2011,"flr_num")
univariate_boot_data=boot_uni(s1_traits_density_2011,"med_flr_date")
univariate_boot_data=boot_uni(s1_traits_density_2011,"dens_Rad_10")
univariate_boot_data=boot_uni(s1_traits_density_2011,"dens_Rad_20")
univariate_boot_data=boot_uni(s1_traits_density_2011,"dens_Rad_30")
univariate_boot_data=boot_uni(s1_traits_density_2011,"dens_Rad_50")
univariate_boot_data=boot_uni(s1_traits_density_2011,"dens_Rad_100")
univariate_boot_data=boot_uni(s1_traits_density_2011,"dens_Rad_200")
#this will calculate the confidence intervals and p-values
univariate_CI_intervals=get_CIs(univariate_boot_data)
get_CIs(univariate_boot_data)
  
  

  
#next I can get multivariate data  

#add density to regressions
  
raw_regression_multivariate <- function(data,trait=c("snd","mid_wid", "len_last", "flr_num","avg_display", "med_flr_date"),y="total_seed" ){

  trait_columns=c(which(colnames(data) %in% trait))
  y_var_pos=which(colnames(data) %in% y)

  myvars=c(y_var_pos, trait_columns)
  #this makes a dataset with only the traits of interest
  subdata=data[,myvars]
  #this removes any incomplete data from the subset data
  subdata=subdata[complete.cases(subdata), ]
  print(nrow(subdata))
  #this finds the columns in the subset data that has the predictor variables
  predictor_positions=c(which(colnames(subdata) %in% trait))

  #this finds the number of predictors I could also just do the length of the trait list
  predictors=ncol(subdata)
    
  for (i in predictor_positions) { 

    #this standardizes the traits of interest,the original trait value is replaced by the standardized trait values
    #This will prepare the traits for regression but then I still need to relativize fitness, which will happen after this loop
    subdata[,i]=(subdata[,i]-mean(subdata[,i], na.rm=TRUE))/sd(subdata[,i], na.rm=TRUE) 
    
  }
  #this relatives fitness
  subdata$rel_fitness=subdata[,1]/mean(subdata[,1])
  
  #now I do the multivariate regression

    regress_raw<-lm(subdata$rel_fitness ~. , data=subdata[,predictor_positions])
    
  c<-capture.output(summary(regress_raw))
  
  return(c)
}


###boot strap the multi
raw_regression_multivariate_boot <- function(bootnum=10, data,trait=c("snd","mid_wid", "len_last", "flr_num","avg_display", "med_flr_date"),y="total_seed" ){  set.seed(2)
  newrow={}
  
  trait_columns=c(which(colnames(data) %in% trait))
  y_var_pos=which(colnames(data) %in% y)
  
  
  myvars=c(y_var_pos, trait_columns)
  #this makes a dataset with only the traits of interest
  subdata=data[,myvars]
  #this removes any incomplete data from the subset data
  subdata=subdata[complete.cases(subdata), ]
  
  #this finds the columns in the subset data that has the predictor variables
  predictor_positions=c(which(colnames(subdata) %in% trait))
  
  
  #this finds the columns in the subset data that has the response variable
  fitness_positions=c(which(colnames(subdata) %in% y))
  
  #this finds the number of predictors I could also just do the length of the trait list
  predictors=ncol(subdata)
  
  newrow={}
for(number in c(1:bootnum)){
      #this creates a randomized dataset
   
    subdata_boot <- subdata[sample(1:nrow(subdata), nrow(subdata), replace=TRUE),]#this samples with replacement within this data set that consists of only the trait and fitness value for the desired trait
      
    for (i in predictor_positions) { 
    
     #this standardizes the traits of interest,the original trait value is replaced by the standardized trait values
      #This will prepare the traits for regression but then I still need to relativize fitness, which will happen after this loop
      subdata_boot[,i]=( subdata_boot[,i]-mean(subdata_boot[,i], na.rm=TRUE))/sd(subdata_boot[,i], na.rm=TRUE) 
    }
       #this relatives fitness
    subdata_boot$rel_fitness=subdata_boot[,fitness_positions]/mean(subdata_boot[,fitness_positions])
      
      #now I do the multivariate regression
      
      regress_raw<-lm(subdata_boot$rel_fitness ~. , data=subdata_boot[,predictor_positions])

  
      names(coef(regress_raw)) 
      newrow=rbind(newrow, coef(regress_raw))
      #colnames(newrow)=c("intercept", paste())    
      }
  return(newrow)
  
}

#function for  CI and pvals
get_mult_CIs = function(data2){
  
  quants=data.frame(row.names=NULL)
  ncol(data2)
    for (i in c(2:ncol(data2))){
      
      q1=quantile(data2[,i], c(0.025, 0.975))
      pval=(length(data2[,i][data2[,i]<=0])/length(data2[,i]))
      lower_CI=as.numeric(q1[1])
      upper_CI=as.numeric(q1[2])
      p_val=pval  
      newrow=c(colnames(data2)[i],lower_CI, upper_CI, p_val)
      quants=as.matrix(rbind(quants,newrow),row.name=NULL)
     
    }  
    quants=as.data.frame(quants, row.names=c(1:length(quants)))
    colnames(quants)=c("trait","lower_CI", "upper_CI", "p_val")

    return(quants)
}


#call function to get multi trait beta estimates
#the default list is the 6 traits
raw_regression_multivariate(s1_traits_density_2011)
raw_regression_multivariate(s1_traits_density_2011,c("snd","mid_wid", "len_last", "flr_num","avg_display", "med_flr_date","dens_Rad_200"))

#call the bootstrap function
mult_dat=raw_regression_multivariate_boot(10000,s1_traits_density_2011,c("snd","mid_wid", "len_last", "flr_num","avg_display", "med_flr_date","dens_Rad_200"))

#call the CI and pvals
get_mult_CIs(mult_dat)






#Do quad terms:
### for single univariate polynomial just input 1 trait at a time


#quadratic estimates function
quad<- function(data,trait=c("snd","mid_wid", "len_last", "flr_num","avg_display", "med_flr_date"),y="total_seed" ){
  newrow={}
  
  trait_columns=c(which(colnames(data) %in% trait))
  y_var_pos=which(colnames(data) %in% y)
  
  
  myvars=c(y_var_pos, trait_columns)
  #this makes a dataset with only the traits of interest
  subdata=data[,myvars]
  #this removes any incomplete data from the subset data
  subdata=subdata[complete.cases(subdata), ]
  
  #this finds the columns in the subset data that has the predictor variables
  predictor_positions=c(which(colnames(subdata) %in% trait))
  
  
  #this finds the columns in the subset data that has the response variable
  fitness_positions=c(which(colnames(subdata) %in% y))
  
  #this finds the number of predictors I could also just do the length of the trait list
  predictors=ncol(subdata)

    for (i in predictor_positions) { 
      
      #this standardizes the traits of interest,the original trait value is replaced by the standardized trait values
      #This will prepare the traits for regression but then I still need to relativize fitness, which will happen after this loop
      subdata[,i]=( subdata[,i]-mean(subdata[,i], na.rm=TRUE))/sd(subdata[,i], na.rm=TRUE) 
    }
    
    #this relatives fitness
    subdata$rel_fitness=subdata[,fitness_positions]/mean(subdata[,fitness_positions])
    
    #now I do the multivariate regression
    
    if (length(predictor_positions)==7){regress_raw<-lm(subdata$rel_fitness ~ . + I(subdata[,predictor_positions[1]]^2)+I(subdata[,predictor_positions[2]]^2)+I(subdata[,predictor_positions[3]]^2)+I(subdata[,predictor_positions[4]]^2)+I(subdata[,predictor_positions[5]]^2)+I(subdata[,predictor_positions[6]]^2)+I(subdata[,predictor_positions[7]]^2), data=subdata[,predictor_positions])} else if(length(predictor_positions)==6) {regress_raw<-lm(subdata$rel_fitness ~. + I(subdata[,predictor_positions[1]]^2)+I(subdata[,predictor_positions[2]]^2)+I(subdata[,predictor_positions[3]]^2)+I(subdata[,predictor_positions[4]]^2)+I(subdata[,predictor_positions[5]]^2)+I(subdata[,predictor_positions[6]]^2), data=subdata[,predictor_positions])} else if(length(predictor_positions)==5) {regress_raw<-lm(subdata$rel_fitness ~. + I(subdata[,predictor_positions[1]]^2)+I(subdata[,predictor_positions[2]]^2)+I(subdata[,predictor_positions[3]]^2)+I(subdata[,predictor_positions[4]]^2)+I(subdata[,predictor_positions[5]]^2), data=subdata[,predictor_positions]) } else if(length(predictor_positions)==4) {regress_raw<-lm(subdata$rel_fitness ~. + I(subdata[,predictor_positions[1]]^2)+I(subdata[,predictor_positions[2]]^2)+I(subdata[,predictor_positions[3]]^2)+I(subdata[,predictor_positions[4]]^2), data=subdata[,predictor_positions])} else if(length(predictor_positions)==3) {regress_raw<-lm(subdata$rel_fitness ~. + I(subdata[,predictor_positions[1]]^2)+I(subdata[,predictor_positions[2]]^2)+I(subdata[,predictor_positions[3]]^2), data=subdata[,predictor_positions])} else if(length(predictor_positions)==2) {regress_raw<-lm(subdata$rel_fitness ~. + I(subdata[,predictor_positions[1]]^2)+I(subdata[,predictor_positions[2]]^2), data=subdata[,predictor_positions])} else if(length(predictor_positions)==1) {regress_raw<-lm(subdata$rel_fitness ~ subdata[,predictor_positions[1]] + I(subdata[,predictor_positions[1]]^2), data=subdata)}
    
  return(regress_raw)
  
}

#quadratic bootstrap
quad_boot<- function(bootnum=10000, data, trait=c("snd","mid_wid", "len_last", "flr_num","avg_display", "med_flr_date"),y="total_seed" ){
  newrow={}
  
  trait_columns=c(which(colnames(data) %in% trait))
  y_var_pos=which(colnames(data) %in% y)
  
  
  myvars=c(y_var_pos, trait_columns)
  #this makes a dataset with only the traits of interest
  subdata=data[,myvars]
  #this removes any incomplete data from the subset data
  subdata=subdata[complete.cases(subdata), ]
  
  #this finds the columns in the subset data that has the predictor variables
  predictor_positions=c(which(colnames(subdata) %in% trait))
  
  
  #this finds the columns in the subset data that has the response variable
  fitness_positions=c(which(colnames(subdata) %in% y))
  
  #this finds the number of predictors I could also just do the length of the trait list
  predictors=ncol(subdata)
  
  newrow={}
  for(i in c(1:bootnum)){
    #this creates a randomized dataset
    
    subdata_boot <- subdata[sample(1:nrow(subdata), nrow(subdata), replace=TRUE),]#this samples with replacement within this data set that consists of only the trait and fitness value for the desired trait
        #this standardizes the traits of interest,the original trait value is replaced by the standardized trait values
      #This will prepare the traits for regression but then I still need to relativize fitness, which will happen after this loop
   
    for (i in predictor_positions) { 
      
      #this standardizes the traits of interest,the original trait value is replaced by the standardized trait values
      #This will prepare the traits for regression but then I still need to relativize fitness, which will happen after this loop
      subdata_boot[,i]=( subdata_boot[,i]-mean(subdata_boot[,i], na.rm=TRUE))/sd(subdata_boot[,i], na.rm=TRUE) 
    }
    
    #this relatives fitness
    subdata_boot$rel_fitness=subdata_boot[,fitness_positions]/mean(subdata_boot[,fitness_positions])

    #now I do the multivariate regression

    if (length(predictor_positions)==7){regress_raw<-lm(subdata_boot$rel_fitness ~ . + I(subdata_boot[,predictor_positions[1]]^2)+I(subdata_boot[,predictor_positions[2]]^2)+I(subdata_boot[,predictor_positions[3]]^2)+I(subdata_boot[,predictor_positions[4]]^2)+I(subdata_boot[,predictor_positions[5]]^2)+I(subdata_boot[,predictor_positions[6]]^2)+I(subdata_boot[,predictor_positions[7]]^2), data=subdata_boot[,predictor_positions])} else if(length(predictor_positions)==6) {regress_raw<-lm(subdata_boot$rel_fitness ~. + I(subdata_boot[,predictor_positions[1]]^2)+I(subdata_boot[,predictor_positions[2]]^2)+I(subdata_boot[,predictor_positions[3]]^2)+I(subdata_boot[,predictor_positions[4]]^2)+I(subdata_boot[,predictor_positions[5]]^2)+I(subdata_boot[,predictor_positions[6]]^2), data=subdata_boot[,predictor_positions])} else if(length(predictor_positions)==5) {regress_raw<-lm(subdata_boot$rel_fitness ~. + I(subdata_boot[,predictor_positions[1]]^2)+I(subdata_boot[,predictor_positions[2]]^2)+I(subdata_boot[,predictor_positions[3]]^2)+I(subdata_boot[,predictor_positions[4]]^2)+I(subdata_boot[,predictor_positions[5]]^2), data=subdata_boot[,predictor_positions]) } else if(length(predictor_positions)==4) {regress_raw<-lm(subdata_boot$rel_fitness ~. + I(subdata_boot[,predictor_positions[1]]^2)+I(subdata_boot[,predictor_positions[2]]^2)+I(subdata_boot[,predictor_positions[3]]^2)+I(subdata_boot[,predictor_positions[4]]^2), data=subdata_boot[,predictor_positions])} else if(length(predictor_positions)==3) {regress_raw<-lm(subdata_boot$rel_fitness ~. + I(subdata_boot[,predictor_positions[1]]^2)+I(subdata_boot[,predictor_positions[2]]^2)+I(subdata_boot[,predictor_positions[3]]^2), data=subdata_boot[,predictor_positions])} else if(length(predictor_positions)==2) {regress_raw<-lm(subdata_boot$rel_fitness ~. + I(subdata_boot[,predictor_positions[1]]^2)+I(subdata_boot[,predictor_positions[2]]^2), data=subdata_boot[,predictor_positions])} else if(length(predictor_positions)==1) {regress_raw<-lm(subdata_boot$rel_fitness ~ subdata_boot[,predictor_positions[1]] + I(subdata_boot[,predictor_positions[1]]^2), data=subdata_boot)}
    
    

    names(coef(regress_raw)) 
    newrow=rbind(newrow, coef(regress_raw))
    #colnames(newrow)=c("intercept", paste())    
  
  }
  return(newrow)
  
}


#function for  CI and pvals
get_quads_CIs = function(data2){
  
  quants=data.frame(row.names=NULL)
  ncol(data2)
  for (i in c(2:ncol(data2))){
    
    q1=quantile(data2[,i], c(0.025, 0.975))
    pval=(length(data2[,i][data2[,i]<=0])/length(data2[,i]))
    lower_CI=as.numeric(q1[1])
    upper_CI=as.numeric(q1[2])
    p_val=pval  
    newrow=c(colnames(data2)[i],lower_CI, upper_CI, p_val)
    quants=as.matrix(rbind(quants,newrow),row.name=NULL)
    
  }  
  quants=as.data.frame(quants, row.names=c(1:length(quants)))
  colnames(quants)=c("trait","lower_CI", "upper_CI", "p_val")
  
  return(quants)
}

#call the quad function on whatever combo of traits you want
quad(s1_traits_density_2011,c("snd","mid_wid", "len_last", "flr_num","avg_display", "med_flr_date","dens_Rad_200"))
#call the bootstrp cfunction
quad_data=quad_boot(bootnum=10000, s1_traits_density_2011,c("snd","mid_wid", "len_last", "flr_num","avg_display", "med_flr_date","dens_Rad_200"))
#call the CI and p val function
get_quads_CIs(quad_data)
