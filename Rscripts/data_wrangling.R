##In these datasets we use density which does not include the focal plant in the density calculation.
save.image("density_selection2014.RData")
load("density_selection2014.RData")


#S4 i have trait data (height length + AND) I have density at 10 30 50 cm and seed data and fruit (but fruit number has to be calculated?


seed_s4=read.csv("./data/Seed_data_s4_2011.csv", header=TRUE)
density_s4=read.csv("./data/s4_density", header=TRUE, )
head(seed_s4)
str(seed_s4)
nrow(seed_s4)
nrow(density_s4)

seed_by_indiv=as.data.frame(aggregate(seed_s4$Seed, by=list(as.factor(seed_s4$ID)), sum))
seed_by_indiv[1:2]
colnames(seed_by_indiv)=c("ID", "seed")
fruit_by_indiv=as.data.frame(aggregate(seed_s4$Fruit.number, by=list(as.factor(seed_s4$ID)), sum))
colnames(fruit_by_indiv)=as.factor(c("ID", "fruit"))
fitness=as.data.frame(merge(fruit_by_indiv, seed_by_indiv, by="ID"))

traits_s4=read.csv("./data/s4_traits.csv", header=TRUE)
traits$ID <- row.names(df1) 

head(traits_s4)

temp=merge(fitness, density_s4, by=c("ID"), all=T)

s4_trht=merge(temp, traits_s4, by=c("ID"), all=T)
s4_trht[s4_trht$ID=="11.01",]
write.csv(s4_trht,file="./output/s4_trht.csv")
s4_trht=read.csv("./output/s4_trht.csv", header=T)

head(s4_trht)
#Let's convert the xy coordinates into denstiy
s4_trht$dens_Rad_10=(s4_trht$X10cm/(pi*0.10*0.10))
s4_trht$dens_Rad_30=(s4_trht$X30cm/(pi*0.30*0.30))
s4_trht$dens_Rad_200=(s4_trht$X10cm/pi*(2.00*2.00))


summary(s4_trht)
#NOW LETS MERGE MARKS DENSITY DATA WITH THE OTHER DATASETS
#s1_s2_s800 traits and merge with XY coordinates for these sties
#in this dataset len_first and len_last is the length of the stem not height from ground.

traits_s1s2=read.csv("./data/s1_s2_800_traits.csv", header=TRUE)
head(traits_s1s2)
height_s1s2=read.csv("./data/s1_s2_800_height.csv", header=TRUE)
head(height_s1s2)

#colnames(height_s1s2) = c( "Pop","ID","height_first","height_last","x","x1")
s1_s2_trht=merge(traits_s1s2, height_s1s2, by=c("Pop","ID"), all=T)
head(s1_s2_trht)
s1_s2_trht$total_seed[is.na(s1_s2_trht$total_seed)]==0

density=read.csv("./data/All_Sites_density2011__Radius_ALL.csv", header=TRUE)
# I went in to the file All_Sites_density2011__Radius_ALL.csv (Friday Dec19) and changed all the sample ID that had letters to lower case.
#s1_s2_ht=read.csv("./data/s1_s2_800_height.csv", header=T)
#s1_s2_tr= read.csv("./data/s1_s2_800_traits.csv", header=T)
#head(density)

density_s1=density[density$Pop=="s1",]
density_s2=density[density$Pop=="s2",]
nrow(density_s1)
nrow(density_s2[density_s2$ID==31,])

s1_s2_trht_s1=s1_s2_trht[s1_s2_trht$Pop=="s1",]
nrow(s1_s2_trht_s1)

s1_s2_trht_s2=s1_s2_trht[s1_s2_trht$Pop=="s2",]
nrow(s1_s2_trht_s2)


s1_traits_density=merge(density_s1, s1_s2_trht_s1, by="ID", all=T)
s2_traits_density=merge(density_s2, s1_s2_trht_s2, by="ID", all=T)
nrow(s1_traits_density)
nrow(s2_traits_density)

head(s1_traits_density)
head(s2_traits_density)

summary(s1_traits_density)
summary(s2_traits_density)

write.csv(s1_traits_density, file="./output/s1_traits_density.csv")
write.csv(s2_traits_density, file="./output/s2_traits_density.csv")

#I did some manual fixing of the s1 dataset dec12 2013
#I did some manual fixing of the s2 dataset in dec12 2013
s1_traits_density=read.csv(file="./output/s1_traits_density_2011.csv", header=T)

#Now I have the trait data and fitness data merged and ready to be merged with surrounding density data.
#s1_s2_trht has the data for s1 s2 and 800s this data needs to be merged with surrounding plant density
#s4_trht has the data for s4 ofr floral traits, seed and fruit number for 10 30 and 50 cm 


##20009 data
s1_traits_2009=read.csv(file="./data/2009/s1_2009_traits.csv", header=T)
s1_density_2009=read.csv(file="./data/2009/s1_density2009__Radius_ALL.csv", header=T)
head(s1_traits_2009)

s1_2009_tr_density=merge(s1_traits_2009,s1_density_2009, by="ID", all=T)

#names(s1_2009_tr_density)[names(s1_2009_tr_density)=="wid"] <- "mid_wid"
#in the output file s1_traits_density_2009, I deleted the middle lobe length and width data and changed the name of the width of the bottom 3 petals fro wid to mid_wid to be 
#consistent with the 2011 data
write.csv(s1_2009_tr_density,file="./output/s1_traits_density_2009.csv")    
