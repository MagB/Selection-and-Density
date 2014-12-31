density interaction
data=s1_traits_density_2011


#this function takes a dataset, density at a particular scale, a trait and fitness info
density_interaction<- function(data,density="dens_Rad_10", trait=c("mid_wid"),y="total_seed" ){
  newrow={}
  if (!(trait %in% names(data))) {stop("trait is not in dataset")}  
  
  trait_columns=c(which(colnames(data) %in% trait))
  y_var_pos=which(colnames(data) %in% y)
  dens_column=which(colnames(data) %in% density)
  
  myvars=c(y_var_pos, dens_column, trait_columns)
  #this makes a dataset with only the traits of interest
  subdata=data[,myvars]
  #this removes any incomplete data from the subset data
  subdata=subdata[complete.cases(subdata), ]
  
  #this finds the columns in the subset data that has the predictor variables
  trait_positions=c(which(colnames(subdata) %in% trait))
  
  
  #this finds the columns in the subset data that has the response variable
  fitness_positions=1
  dens_position=2

  
  for (i in (dens_pos:max(trait_positions))) { 
  
    #this standardizes density and the traits of interest,the original trait value is replaced by the standardized trait values
    #This will prepare the traits for regression but then I still need to relativize fitness, which will happen after this loop
    subdata[,i]=( subdata[,i]-mean(subdata[,i], na.rm=TRUE))/sd(subdata[,i], na.rm=TRUE) 
  }
  
  #this relatives fitness
  subdata$rel_fitness=subdata[,fitness_positions]/mean(subdata[,fitness_positions])
  
  regress_raw<-lm(subdata$rel_fitness ~ subdata[,dens_position[1]] + subdata[,trait_positions[1]] + subdata[,trait_positions[1]]* subdata[,dens_position[1]], data=subdata)
  c<-capture.output(summary(regress_raw))  
  #now I do the multivariate regression
  
  #if (length(predictor_positions)==7){regress_raw<-lm(subdata$rel_fitness ~ . + I(subdata[,predictor_positions[1]]^2)+I(subdata[,predictor_positions[2]]^2)+I(subdata[,predictor_positions[3]]^2)+I(subdata[,predictor_positions[4]]^2)+I(subdata[,predictor_positions[5]]^2)+I(subdata[,predictor_positions[6]]^2)+I(subdata[,predictor_positions[7]]^2), data=subdata[,predictor_positions])} else if(length(predictor_positions)==6) {regress_raw<-lm(subdata$rel_fitness ~. + I(subdata[,predictor_positions[1]]^2)+I(subdata[,predictor_positions[2]]^2)+I(subdata[,predictor_positions[3]]^2)+I(subdata[,predictor_positions[4]]^2)+I(subdata[,predictor_positions[5]]^2)+I(subdata[,predictor_positions[6]]^2), data=subdata[,predictor_positions])} else if(length(predictor_positions)==5) {regress_raw<-lm(subdata$rel_fitness ~. + I(subdata[,predictor_positions[1]]^2)+I(subdata[,predictor_positions[2]]^2)+I(subdata[,predictor_positions[3]]^2)+I(subdata[,predictor_positions[4]]^2)+I(subdata[,predictor_positions[5]]^2), data=subdata[,predictor_positions]) } else if(length(predictor_positions)==4) {regress_raw<-lm(subdata$rel_fitness ~. + I(subdata[,predictor_positions[1]]^2)+I(subdata[,predictor_positions[2]]^2)+I(subdata[,predictor_positions[3]]^2)+I(subdata[,predictor_positions[4]]^2), data=subdata[,predictor_positions])} else if(length(predictor_positions)==3) {regress_raw<-lm(subdata$rel_fitness ~. + I(subdata[,predictor_positions[1]]^2)+I(subdata[,predictor_positions[2]]^2)+I(subdata[,predictor_positions[3]]^2), data=subdata[,predictor_positions])} else if(length(predictor_positions)==2) {regress_raw<-lm(subdata$rel_fitness ~. + I(subdata[,predictor_positions[1]]^2)+I(subdata[,predictor_positions[2]]^2), data=subdata[,predictor_positions])} else if(length(predictor_positions)==1) {regress_raw<-lm(subdata$rel_fitness ~ subdata[,predictor_positions[1]] + I(subdata[,predictor_positions[1]]^2), data=subdata)}
  print(names(subdata))
  return(c)
  
}


#this function does the bootstrapping
density_interaction_boot<- function(bootnum=10, data,density="dens_Rad_10", trait=c("mid_wid"),y="total_seed" ){
  newrow={}
  if (!(trait %in% names(data))) {stop("trait is not in dataset")}  
  
  trait_columns=c(which(colnames(data) %in% trait))
  y_var_pos=which(colnames(data) %in% y)
  dens_column=which(colnames(data) %in% density)
  
  myvars=c(y_var_pos, dens_column, trait_columns)
  #this makes a dataset with only the traits of interest
  subdata=data[,myvars]
  #this removes any incomplete data from the subset data
  subdata=subdata[complete.cases(subdata), ]
  
  #this finds the columns in the subset data that has the predictor variables
  trait_positions=c(which(colnames(subdata) %in% trait))
  
  
  #this finds the columns in the subset data that has the response variable
  fitness_positions=1
  dens_position=2
  
  newrow={}
  for(number in c(1:bootnum)){
    #this creates a randomized dataset
    
    subdata_boot <- subdata[sample(1:nrow(subdata), nrow(subdata), replace=TRUE),]#this samples with replacement within this data set that consists of only the trait and fitness value for the desired trait
    
  for (i in (dens_pos:max(trait_positions))) { 
    
    #this standardizes density and the traits of interest,the original trait value is replaced by the standardized trait values
    #This will prepare the traits for regression but then I still need to relativize fitness, which will happen after this loop
    subdata_boot[,i]=( subdata_boot[,i]-mean(subdata_boot[,i], na.rm=TRUE))/sd(subdata_boot[,i], na.rm=TRUE) 
  }
  
  #this relatives fitness
    subdata_boot$rel_fitness=subdata_boot[,fitness_positions]/mean(subdata_boot[,fitness_positions])
  
    
    #now do the regression
  regress_raw<-lm(subdata_boot$rel_fitness ~ subdata_boot[,dens_position[1]] + subdata_boot[,trait_positions[1]] + subdata_boot[,trait_positions[1]]* subdata_boot[,dens_position[1]], data=subdata_boot)
    names(coef(regress_raw)) 
    newrow=rbind(newrow, coef(regress_raw))
    #colnames(newrow)=c("intercept", paste())    
  }
  colnames(newrow)=c("intercept", c(density, trait,"interaction")) 
  return(newrow)
  
}

#this function gets CI and pvals
get_CIs = function(data2){
  
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

density_interaction(s1_traits_density_2011,density="dens_Rad_10", trait=c("snd"),y="total_seed" )
dens_boot_data=density_interaction_boot(10,s1_traits_density_2011,density="dens_Rad_10", trait=c("snd"),y="total_seed" )
get_CIs(dens_boot_data)
