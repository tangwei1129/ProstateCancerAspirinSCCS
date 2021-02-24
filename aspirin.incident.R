# source('code.R')


library(survival)
library(knitr)
library(haven)


## data import ##
SASData <- read_sas("aspirin_prostatecan_aug2020.sas7bdat")
#write.csv(SASData,"aspirin_prostatecan.csv") 

d=SASData
#d=read.csv("aspirin_prostatecan.csv",header=T)

dat=d
dat[dat == "9999"]=NA
#dat[dat == "7777"]=0
dat[dat == "8888"]=NA
 

###################################################################################################
data=data.frame(dat)
# bmi convert #
bmi=cut(data$BMI, breaks=c(0,18.5,25,30,35,40,100),right=F)
data$BMI=as.numeric(bmi)

# combine fa.x and br.y prostate history, with Fa_Cancer.z #
history_combine=function(x,y,z){
  x=ifelse(z==0,0,x)
  #y=ifelse(z==0,0,y)
  
  x[x=='7777'] = 0
  y[y=='7777'] = 0
  
  e=ifelse(x==1|y==1,1,NA)
  e=ifelse(is.na(x) & is.na(y) ,-9,e)
  e[is.na(e)] = 0
  e[e == -9] = NA
  return(e)
}


data$family=history_combine(data$Fa_ProstateCancer,data$Br_ProstateCancer,data$Fa_Cancer)

# combine OTCNSAIDS.x and RxNSAIDS.y as nsaids #
nsaid_combine=function(x,y){
  z=ifelse(is.na(x),y,x)
  z[z=='7777']=NA
  return(z)
}
data$nsaids=nsaid_combine(data$OTCNSAIDS,data$RxNSAIDS)


# re-define Asprin use #
# combine 2 variables #
variable_combine=function(x,y){
  x[is.na(x)]= -9
  y[is.na(y)]= -9
  z=ifelse( (x+y) %in% c(-8,1,2), 1,ifelse( (x+y) == 0,0,NA))
  return(z)
}

## any aspirin ##
data$Aspirin_any = variable_combine(data$Aspirin, data$LowAspirin )
##regular aspirin ##
data$Aspirin_reg = ifelse(!data$Aspirin == 1 & data$Aspirin_any ==1, NA, data$Aspirin_any )

## total 22426 enrollment, Prostate_Cancer AgeDX , Cancer Registry 1058 case ##
## vital month 156 median, 143 mean, 6627 death ##

data$pr=ifelse(data$Prostate_Cancer==1,"PrCa",NA)
data$pr[is.na(data$pr)]="NoPrCa"
data$pr=as.numeric(as.factor(data$pr))-1

## drop 24 month to reflect the difference of case and control in 2016.case 2018.control ##
## Currently NDI deaths are complete through 12/31/2018 and 
# the 12 SCCS state cancer registries are reporting to be complete through at least 12/31/2016. 

# For the Cox analyses in the study, the length of follow-up time is not comparable between cases and non-cases. 
# The tumor diagnoses information we got from state registries usually fall behind the vital status linkage. 
# Most of the registries currently have complete data through the end of year 2015, 
# whereas the vital status data are through 2016. The end date of follow-up needs to be adjusted for non-cases.

data$t1=ifelse(is.na(data$Prostate_Cancer), as.numeric(data$VitalStatus_Months)-24,as.numeric(data$VitalStatus_Months))
data$time=ifelse(is.na(data$Incidence_Months),data$t1,data$Incidence_Months)
data$time=ifelse(data$time > 0, data$time, NA)

#write.csv(data,"aspirin_prostatecan_cl.csv") 

################################################################################################################
#### function to get stat from cox ####
library(survival)
cox=function(var, data){
	data=subset(data,time > 0 )
	cox.sum=summary(coxph(Surv(time,pr)~data[,var]+Enrollment_Age,data=data))
	p=round(cox.sum$coefficients[1,"Pr(>|z|)"],digits=3)
	hr=round(cox.sum$coefficients[1,"exp(coef)"],digits=2)
	ci=round(cox.sum$conf.int[1,c("lower .95", "upper .95")],digits=2)
	hr.ci.p=paste0(hr,"(",ci[1],"-",ci[2],")"," p=",p)
	return(hr.ci.p)
}
 #cox(var="Aspirin_re",data=data)

## cox.adj adjust for age,bmi,edu,hhincome,fa_history,smk,diabetes,bph,statin,stroke,psa.scr,dre.scr, acet, nsaids ##
cox.adj=function(var, data){
	data=subset(data,time > 0 )
	data$SmokingStatus = ifelse(data$SmokingStatus == '7777',NA,data$SmokingStatus)
	smk <- model.matrix.lm(~factor(SmokingStatus,levels=c(3,2,1)), data = data,na.action="na.pass")
	cox.sum=summary(coxph(Surv(time,pr)~ data[,var]+ Enrollment_Age + BMI + Education +
	                        HHIncome + family + smk[,2]+smk[,3] + Diabetes + BPH
	                        #+StatinUse 
	                        + Stroke_TIA + PSA_Ever+DRE_Ever + Acetaminophen+nsaids,data=data))
	
	p=round(cox.sum$coefficients[1,"Pr(>|z|)"],digits=3)
	hr=round(cox.sum$coefficients[1,"exp(coef)"],digits=2)
	ci=round(cox.sum$conf.int[1,c("lower .95", "upper .95")],digits=2)
	hr.ci.p=paste0(hr,"(",ci[1],"-",ci[2],")"," p=",p)
	return(hr.ci.p)
}
#cox.adj(var="Aspirin_re",data=data)

#######################################################################################################
### construct aspirin use table ###
#### start of function rr.table #####

rr.table=function(temp){
n=nrow(temp) #n=22426
asp.never=subset(temp,Aspirin_re==0)
asp.use=subset(temp,Aspirin_re==1)
#asp.use.1=subset(asp.use,Aspirin_PillsPerWk <= 7)
asp.use.1=subset(asp.use,Aspirin_PillsPerWk < 7)

#asp.use.gt1=subset(asp.use,Aspirin_PillsPerWk > 7)
asp.use.gt1=subset(asp.use,Aspirin_PillsPerWk >= 7)

asp.use.3= subset(asp.use,Aspirin_Yrs <= 3)
asp.use.gt3= subset(asp.use,Aspirin_Yrs > 3)

asp=rbind(table(asp.never$pr),
	table(asp.use$pr),
	table(asp.use.1$pr),
	table(asp.use.gt1$pr),
	table(asp.use.3$pr),
	table(asp.use.gt3$pr)
	)
colnames(asp)=c("Control","Case")
rownames(asp)=c("Never","Use","less than Daily","Daily and more","<=3 years",">3 years")
asp=data.frame(asp)
asp$PCT=round(100*asp$Case/(asp$Case+asp$Control),digits=2)
asp$py=c(sum(asp.never$time,na.rm = T),sum(asp.use$time,na.rm = T),sum(asp.use.1$time,na.rm = T),sum(asp.use.gt1$time,na.rm = T),sum(asp.use.3$time,na.rm = T),sum(asp.use.gt3$time,na.rm = T))/12

####################################################### 
#### RR table 2  #####################################
#cox(var="Aspirin_re",data=data)
#cox.adj(var="Aspirin_re",data=data)


temp$asp.use=ifelse(temp$Aspirin_re %in% c(0,1),1,0)## asp.use=1
## pills per week ##
temp$asp.use.1=ifelse(temp$Aspirin_re == 0 | (temp$Aspirin_re == 1 & temp$Aspirin_PillsPerWk < 7),1,0)
temp$asp.use.gt1=ifelse(temp$Aspirin_re == 0 | (temp$Aspirin_re == 1 & temp$Aspirin_PillsPerWk >= 7),1,0)
temp$asp.trend=ifelse(temp$Aspirin_re ==0,0,ifelse(temp$Aspirin_PillsPerWk < 7,1,ifelse(temp$Aspirin_PillsPerWk >= 7,2,NA)))
## aspirin yrs ##
temp$asp.use.3=ifelse(temp$Aspirin_re == 0 | (temp$Aspirin_re == 1 & temp$Aspirin_Yrs <= 3),1,0)
temp$asp.use.gt3=ifelse(temp$Aspirin_re == 0 | (temp$Aspirin_re == 1 & temp$Aspirin_Yrs > 3),1,0)
temp$asp.yr.trend=ifelse(temp$Aspirin_re ==0,0,ifelse(temp$Aspirin_Yrs <= 3,1,ifelse(temp$Aspirin_Yrs >3,2,NA)))

# asp.use #
test=subset(temp,asp.use==1)
rr=c(cox(var="Aspirin_re",data=test),cox.adj(var="Aspirin_re",data=test))

# by pill #
test=subset(temp,asp.use.1==1)
rr.1=c(cox(var="Aspirin_re",data=test),cox.adj(var="Aspirin_re",data=test))

test=subset(temp,asp.use.gt1==1)
rr.gt1=c(cox(var="Aspirin_re",data=test),cox.adj(var="Aspirin_re",data=test))

test=temp
rr.1.trend =c(cox(var="asp.trend",data=test),cox.adj(var="asp.trend",data=test))

# by yr #
test=subset(temp,asp.use.3==1)
rr.3=c(cox(var="Aspirin_re",data=test),cox.adj(var="Aspirin_re",data=test))

test=subset(temp,asp.use.gt3==1)
rr.gt3=c(cox(var="Aspirin_re",data=test),cox.adj(var="Aspirin_re",data=test))

test=temp
rr.3.trend =c(cox(var="asp.yr.trend",data=test),cox.adj(var="asp.yr.trend",data=test))


asp$RR = c("",rr[1],rr.1[1],rr.gt1[1],rr.3[1],rr.gt3[1])
asp$RR.adj= c("",rr[2],rr.1[2],rr.gt1[2],rr.3[2],rr.gt3[2])
asp$RR.trend=c("","",rr.1.trend[1],"","",rr.3.trend[1])
asp$RR.trend.adj=c("","",rr.1.trend[2],"","",rr.3.trend[2])
return(asp)
}
##### end of function #############################################################

data$Aspirin_re = data$Aspirin_any

data$Aspirin_Yrs = ifelse(data$Aspirin_Yrs =='7777', NA, data$Aspirin_Yrs)
data$LowAspirin_Yrs = ifelse(data$LowAspirin_Yrs =='7777', NA, data$LowAspirin_Yrs)

data$Aspirin_PillsPerWk = ifelse(data$Aspirin_PillsPerWk =='7777', NA, data$Aspirin_PillsPerWk)
data$LowAspirin_PillsPerWk = ifelse(data$LowAspirin_PillsPerWk =='7777', NA, data$LowAspirin_PillsPerWk)

data$Aspirin_Yrs = ifelse(is.na(data$Aspirin_Yrs), data$LowAspirin_Yrs, data$Aspirin_Yrs)
data$Aspirin_PillsPerWk = ifelse(is.na(data$Aspirin_PillsPerWk), data$LowAspirin_PillsPerWk, data$Aspirin_PillsPerWk)


temp=data
asp.rr=rr.table(temp)

# low stage #
temp=data
temp$pr=ifelse(temp$SCCS_SummAJCCStage %in% c("I","II"),1,0)
asp.rr.lowstage=rr.table(temp)

# high stage #
temp=data
temp$pr=ifelse(temp$SCCS_SummAJCCStage %in% c("III","IV"),1,0)
asp.rr.highstage=rr.table(temp)


# low Gleason  #
temp=data
temp$pr=ifelse(temp$Earliest_Gleason_Score %in% c("002","003","004","005","006","007"),1,0)
asp.rr.lowGleason=rr.table(temp)


# high Gleason #
temp=data
temp$pr=ifelse(temp$Earliest_Gleason_Score %in% c("008","009","010"),1,0)
asp.rr.highGleason=rr.table(temp)

# aggressive disease #
temp=data
highG = c("008","009","010")
temp$pr=ifelse(temp$Earliest_Gleason_Score %in% highG | temp$SCCS_SummAJCCStage == 'IV'  |
                 temp$TNM_Clin_M %in% c('1','1B','1C') |	temp$TNM_Path_M %in% c('1A','1B') |
                 temp$TNM_Clin_N == '1' | 	temp$TNM_Path_N == '1' |	
                 temp$TNM_Clin_T == '4' |  temp$TNM_Path_T == '4',1,0)
asp.rr.aggressive=rr.table(temp)



table2=qpcR:::rbind.na(asp.rr, "Low Stage", asp.rr.lowstage, "High Stage", asp.rr.highstage,
				"Low Gleason", asp.rr.lowGleason, "High Gleason", asp.rr.highGleason,"aggressive",asp.rr.aggressive )

write.csv(table2,"table2.any.csv")


###############################################
###############################################
###############################################

#### regular #################################


###################################################################################################
data=data.frame(dat)
# bmi convert #
bmi=cut(data$BMI, breaks=c(0,18.5,25,30,35,40,100),right=F)
data$BMI=as.numeric(bmi)

# combine fa.x and br.y prostate history, with Fa_Cancer.z #
history_combine=function(x,y,z){
  x=ifelse(z==0,0,x)
  #y=ifelse(z==0,0,y)
  
  x[x=='7777'] = 0
  y[y=='7777'] = 0
  
  e=ifelse(x==1|y==1,1,NA)
  e=ifelse(is.na(x) & is.na(y) ,-9,e)
  e[is.na(e)] = 0
  e[e == -9] = NA
  return(e)
}


data$family=history_combine(data$Fa_ProstateCancer,data$Br_ProstateCancer,data$Fa_Cancer)

# combine OTCNSAIDS.x and RxNSAIDS.y as nsaids #
nsaid_combine=function(x,y){
  z=ifelse(is.na(x),y,x)
  z[z=='7777']=NA
  return(z)
}
data$nsaids=nsaid_combine(data$OTCNSAIDS,data$RxNSAIDS)


# re-define Asprin use #
# combine 2 variables #
variable_combine=function(x,y){
  x[is.na(x)]= -9
  y[is.na(y)]= -9
  z=ifelse( (x+y) %in% c(-8,1,2), 1,ifelse( (x+y) == 0,0,NA))
  return(z)
}

## any aspirin ##
data$Aspirin_any = variable_combine(data$Aspirin, data$LowAspirin )
##regular aspirin ##
data$Aspirin_reg = ifelse(!data$Aspirin == 1 & data$Aspirin_any ==1, NA, data$Aspirin_any )

## total 22426 enrollment, Prostate_Cancer AgeDX , Cancer Registry 1058 case ##
## vital month 156 median, 143 mean, 6627 death ##

data$pr=ifelse(data$Prostate_Cancer==1,"PrCa",NA)
data$pr[is.na(data$pr)]="NoPrCa"
data$pr=as.numeric(as.factor(data$pr))-1

## drop 24 month to reflect the difference of case and control in 2016.case 2018.control ##
## Currently NDI deaths are complete through 12/31/2018 and 
# the 12 SCCS state cancer registries are reporting to be complete through at least 12/31/2016. 

# For the Cox analyses in the study, the length of follow-up time is not comparable between cases and non-cases. 
# The tumor diagnoses information we got from state registries usually fall behind the vital status linkage. 
# Most of the registries currently have complete data through the end of year 2015, 
# whereas the vital status data are through 2016. The end date of follow-up needs to be adjusted for non-cases.

data$t1=ifelse(is.na(data$Prostate_Cancer), as.numeric(data$VitalStatus_Months)-24,as.numeric(data$VitalStatus_Months))
data$time=ifelse(is.na(data$Incidence_Months),data$t1,data$Incidence_Months)
data$time=ifelse(data$time > 0, data$time, NA)

##############
data$Aspirin_re = data$Aspirin_reg

data$Aspirin_Yrs = ifelse(data$Aspirin_Yrs =='7777', NA, data$Aspirin_Yrs)
data$Aspirin_PillsPerWk = ifelse(data$Aspirin_PillsPerWk =='7777', NA, data$Aspirin_PillsPerWk)

temp=data
asp.rr=rr.table(temp)


# low stage #
temp=data
temp$pr=ifelse(temp$SCCS_SummAJCCStage %in% c("I","II"),1,0)
asp.rr.lowstage=rr.table(temp)

# high stage #
temp=data
temp$pr=ifelse(temp$SCCS_SummAJCCStage %in% c("III","IV"),1,0)
asp.rr.highstage=rr.table(temp)


# low Gleason  #
temp=data
temp$pr=ifelse(temp$Earliest_Gleason_Score %in% c("002","003","004","005","006","007"),1,0)
asp.rr.lowGleason=rr.table(temp)


# high Gleason #
temp=data
temp$pr=ifelse(temp$Earliest_Gleason_Score %in% c("008","009","010"),1,0)
asp.rr.highGleason=rr.table(temp)

# aggressive disease #
temp=data
highG = c("008","009","010")
temp$pr=ifelse(temp$Earliest_Gleason_Score %in% highG | temp$SCCS_SummAJCCStage == 'IV'  |
                 temp$TNM_Clin_M %in% c('1','1B','1C') |	temp$TNM_Path_M %in% c('1A','1B') |
                 temp$TNM_Clin_N == '1' | 	temp$TNM_Path_N == '1' |	
                 temp$TNM_Clin_T == '4' |  temp$TNM_Path_T == '4',1,0)
asp.rr.aggressive=rr.table(temp)



table2.reg=qpcR:::rbind.na(asp.rr, "Low Stage", asp.rr.lowstage, "High Stage", asp.rr.highstage,
                       "Low Gleason", asp.rr.lowGleason, "High Gleason", asp.rr.highGleason,"aggressive",asp.rr.aggressive )

write.csv(table2.reg,"table2.reg.csv")

