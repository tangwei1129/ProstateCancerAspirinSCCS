# source('code.R')

# source('code.R')


library(survival)
library(knitr)
library(haven)
library(cmprsk)

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
####
data$t1=ifelse(is.na(data$Prostate_Cancer), as.numeric(data$VitalStatus_Months)-24,as.numeric(data$VitalStatus_Months))
data$time=ifelse(is.na(data$Incidence_Months),data$t1,data$Incidence_Months)
data$time=ifelse(data$time > 0, data$time, NA)

###############################################################################################################
###

### overall survival  ###
#switch to higher version R for plotting#
library(survival)
library("survminer")
############################# function  #######################

### variable with OS ###
cox.OS.var=function(data,var){
  survival_data=data
  survival_data=subset(data,VitalStatus_Months > 0 & time > 0 )
  #survival_data$vitalstatus_new=ifelse(survival_data$Year_Death > 2015 | is.na(survival_data$Year_Death),"A",survival_data$VitalStatus)
  survival_data$death=as.numeric(as.factor(survival_data$VitalStatus))-1
  survival_data[,var] = ifelse(survival_data[,var] == '7777',NA,survival_data[,var])
  if(length(levels(factor(data[,var]))) ==2) {
    cox.sum=summary(coxph(Surv(VitalStatus_Months, death)~ factor(get(var)),data=survival_data))
  } else{cox.sum=summary(coxph(Surv(VitalStatus_Months, death)~ get(var),data=survival_data))}
  p=round(cox.sum$coefficients[1,"Pr(>|z|)"],digits=3)
  hr=round(cox.sum$coefficients[1,"exp(coef)"],digits=2)
  ci=round(cox.sum$conf.int[1,c("lower .95", "upper .95")],digits=2)
  hr.ci.p=paste0(hr,"(",ci[1],"-",ci[2],")"," p=",p)
  print(paste(var, hr.ci.p))
  #return(hr.ci.p)
  #return(cox.sum)
}


### event time ####

event.person.mo=function(data,var) {
  survival_data=subset(data, time > 0 & VitalStatus_Months > 0 & !(is.na(data[,var])))
	survival_data$VitalStatus=ifelse(survival_data$Year_Death > yrs | is.na(survival_data$Year_Death),"A",survival_data$VitalStatus)
	survival_data$death=as.numeric(as.factor(survival_data$VitalStatus))-1
      tl=table(survival_data$death,survival_data[,var])
	set0=sum(subset(survival_data,  get(var) == 0)$VitalStatus_Months	)
      set1=sum(subset(survival_data,  get(var) == 1)$VitalStatus_Months	)
	#sum.tl= qpcR:::rbind.na(c(tl[2,1],set0),c(tl[2,2],set1))
	sum.tl= rbind(c(tl[2,1],set0),c(tl[2,2],set1))
      colnames(sum.tl)=c("events","person-month")
	rownames(sum.tl)=c("No","Yes")
      return(sum.tl)
}

### OS use Aspirin_re ###
cox.OS=function(data){
	survival_data=subset(data,VitalStatus_Months > 0 & time > 0 & !is.na(Aspirin_re) )
	survival_data$VitalStatus=ifelse(survival_data$Year_Death > yrs | is.na(survival_data$Year_Death),"A",survival_data$VitalStatus)
	survival_data$death=as.numeric(as.factor(survival_data$VitalStatus))-1
	cox.sum=summary(coxph(Surv(VitalStatus_Months, death)~Aspirin_re + Enrollment_Age ,data=survival_data))
	p=round(cox.sum$coefficients[1,"Pr(>|z|)"],digits=3)
	hr=round(cox.sum$coefficients[1,"exp(coef)"],digits=2)
	ci=round(cox.sum$conf.int[1,c("lower .95", "upper .95")],digits=2)
	hr.ci.p=paste0(hr,"(",ci[1],"-",ci[2],")"," p=",p)
	return(hr.ci.p)
	#return(cox.sum)
}

cox.OS.age=function(data){
  survival_data=subset(data,VitalStatus_Months > 0 & time > 0 & !is.na(Aspirin_re) )
  survival_data$VitalStatus=ifelse(survival_data$Year_Death > yrs | is.na(survival_data$Year_Death),"A",survival_data$VitalStatus)
  survival_data$death=as.numeric(as.factor(survival_data$VitalStatus))-1
  cox.sum=summary(coxph(Surv(Enrollment_Age, VitalStatus_Age, death) ~ Aspirin_re,data=survival_data))
  p=round(cox.sum$coefficients[1,"Pr(>|z|)"],digits=3)
  hr=round(cox.sum$coefficients[1,"exp(coef)"],digits=2)
  ci=round(cox.sum$conf.int[1,c("lower .95", "upper .95")],digits=2)
  hr.ci.p=paste0(hr,"(",ci[1],"-",ci[2],")"," p=",p)
  return(hr.ci.p)
  #return(cox.sum)
}

cox.OS.adj=function(data){
	survival_data=subset(data,VitalStatus_Months > 0 & time > 0 &!is.na(Aspirin_re) )
	survival_data$VitalStatus=ifelse(survival_data$Year_Death > yrs | is.na(survival_data$Year_Death),"A",survival_data$VitalStatus)
	survival_data$death=as.numeric(as.factor(survival_data$VitalStatus))-1
	#stage=ifelse(survival_data$SCCS_SummAJCCStage %in% c("III","IV"),"High",ifelse(survival_data$SCCS_SummAJCCStage %in% c("I","II"),"Low",NA))
	#gl=ifelse(survival_data$Earliest_Gleason_Score %in% c("002","003","004","005","006","007"),"Low",ifelse(survival_data$Earliest_Gleason_Score %in% c("008","009","010"),"High",NA))
	
	survival_data$SmokingStatus = ifelse(survival_data$SmokingStatus == '7777',NA,survival_data$SmokingStatus)
	smk <- model.matrix.lm(~factor(SmokingStatus,levels=c(3,2,1)), data = survival_data,na.action="na.pass")
	
	aggressive = ifelse(survival_data$Earliest_Gleason_Score %in% c("008","009","010") | survival_data$SCCS_SummAJCCStage == 'IV'  |
	                      survival_data$TNM_Clin_M %in% c('1','1B','1C') |	survival_data$TNM_Path_M %in% c('1A','1B') |
	                      survival_data$TNM_Clin_N == '1' | 	survival_data$TNM_Path_N == '1' |	
	                      survival_data$TNM_Clin_T == '4' |  survival_data$TNM_Path_T == '4',1,0)
	
	cox.sum=summary(coxph(Surv(VitalStatus_Months, death)~ Aspirin_re 
		                        + Enrollment_Age 
	                          + BMI + Education 
		                        #+ HHIncome 
	                          + family 
	                          + smk[,2]+smk[,3] 
	                          + Diabetes + BPH
		                        #+ StatinUse 
	                          #+ Stroke_TIA 
		                        + factor(Enrollment_Year)
		                        + aggressive
		                        + PSA_Ever 
	                          + DRE_Ever + Acetaminophen + nsaids
	                          , data=survival_data))
		
	p=round(cox.sum$coefficients[1,"Pr(>|z|)"],digits=3)
	hr=round(cox.sum$coefficients[1,"exp(coef)"],digits=2)
	ci=round(cox.sum$conf.int[1,c("lower .95", "upper .95")],digits=2)
	hr.ci.p=paste0(hr,"(",ci[1],"-",ci[2],")"," p=",p)
	return(hr.ci.p)
	#return(cox.sum)
}
#data$Aspirin_re=data$Aspirin_any
#data$Aspirin_re=data$Aspirin_reg
#cox.OS.adj(data)
cox.OS.adj.age=function(data){
  survival_data=subset(data,VitalStatus_Months > 0 & time > 0 &!is.na(Aspirin_re) )
  survival_data$VitalStatus=ifelse(survival_data$Year_Death > yrs | is.na(survival_data$Year_Death),"A",survival_data$VitalStatus)
  survival_data$death=as.numeric(as.factor(survival_data$VitalStatus))-1
  #stage=ifelse(survival_data$SCCS_SummAJCCStage %in% c("III","IV"),"High",ifelse(survival_data$SCCS_SummAJCCStage %in% c("I","II"),"Low",NA))
  #gl=ifelse(survival_data$Earliest_Gleason_Score %in% c("002","003","004","005","006","007"),"Low",ifelse(survival_data$Earliest_Gleason_Score %in% c("008","009","010"),"High",NA))
  
  survival_data$SmokingStatus = ifelse(survival_data$SmokingStatus == '7777',NA,survival_data$SmokingStatus)
  smk <- model.matrix.lm(~factor(SmokingStatus,levels=c(3,2,1)), data = survival_data,na.action="na.pass")
  
  aggressive = ifelse(survival_data$Earliest_Gleason_Score %in% c("008","009","010") | survival_data$SCCS_SummAJCCStage == 'IV'  |
                        survival_data$TNM_Clin_M %in% c('1','1B','1C') |	survival_data$TNM_Path_M %in% c('1A','1B') |
                        survival_data$TNM_Clin_N == '1' | 	survival_data$TNM_Path_N == '1' |	
                        survival_data$TNM_Clin_T == '4' |  survival_data$TNM_Path_T == '4',1,0)
  
  cox.sum=summary(coxph(Surv(Enrollment_Age,VitalStatus_Age, death)~ Aspirin_re 
                        + BMI + Education 
                        #+ HHIncome 
                        + family 
                        + smk[,2]+smk[,3] 
                        + Diabetes + BPH
                        #+ StatinUse 
                        #+ Stroke_TIA 
                        + factor(Enrollment_Year)
                        + aggressive
                        + PSA_Ever 
                        + DRE_Ever + Acetaminophen + nsaids
                        , data=survival_data))
  
  p=round(cox.sum$coefficients[1,"Pr(>|z|)"],digits=3)
  hr=round(cox.sum$coefficients[1,"exp(coef)"],digits=2)
  ci=round(cox.sum$conf.int[1,c("lower .95", "upper .95")],digits=2)
  hr.ci.p=paste0(hr,"(",ci[1],"-",ci[2],")"," p=",p)
  return(hr.ci.p)
  #return(cox.sum)
}



cmp=function(cif){
  cif=subset(cif,VitalStatus_Months > 0 & time > 0 & !is.na(Aspirin_re) )
  cif$VitalStatus=ifelse(cif$Year_Death > yrs | is.na(cif$Year_Death),"A",cif$VitalStatus)
  cif$status= ifelse( cif$VitalStatus == 'A', 0, ifelse(cif$VitalStatus == 'D',1,NA))
  cif$status= ifelse(cif$VitalStatus == 'D' & cif$CauseOfDeath == '', 2, cif$status)
  
  ### Fine-Gray model ###
  cov1=model.matrix.lm(~Aspirin_re + Enrollment_Age, data = cif, na.action = 'na.pass')[,-1]
  #cov1= model.matrix.lm(~Aspirin_re +  Education + HHIncome+ Enrollment_Age+BMI+ PSA_Ever  + Fa_ProstateCancer, data= cif, na.action = 'na.pass')[,-1]
  fg_model= crr(cif$VitalStatus_Months,cif$status,cov1=cov1)
  fg_model.sum = summary(fg_model)
  p=round(fg_model.sum$coef[1,"p-value"],digits=3)
  hr=round(fg_model.sum$coef[1,"exp(coef)"],digits=2)
  ci=round(fg_model.sum$conf.int[1,c("2.5%", "97.5%")],digits=2)
  hr.ci.p=paste0(hr,"(",ci[1],"-",ci[2],")"," p=",p)
  return(hr.ci.p)
}
cmp.adj=function(cif){
  cif=subset(cif,VitalStatus_Months > 0 & time > 0 & !is.na(Aspirin_re) )
  cif$VitalStatus=ifelse(cif$Year_Death > yrs | is.na(cif$Year_Death),"A",cif$VitalStatus)
  cif$status= ifelse( cif$VitalStatus == 'A', 0, ifelse(cif$VitalStatus == 'D',1,NA))
  cif$status= ifelse(cif$VitalStatus == 'D' & cif$CauseOfDeath == '', 2, cif$status)
  
  ### Fine-Gray model ###
  #cov1=model.matrix.lm(~Aspirin_re + Enrollment_Age, data = cif, na.action = 'na.pass')[,-1]
  smk <- model.matrix.lm(~factor(SmokingStatus,levels=c(3,2,1)), data = cif,na.action="na.pass")
  
  aggressive = ifelse(cif$Earliest_Gleason_Score %in% c("008","009","010") | cif$SCCS_SummAJCCStage == 'IV'  |
                        cif$TNM_Clin_M %in% c('1','1B','1C') |	cif$TNM_Path_M %in% c('1A','1B') |
                        cif$TNM_Clin_N == '1' | 	cif$TNM_Path_N == '1' |	
                        cif$TNM_Clin_T == '4' |  cif$TNM_Path_T == '4',1,0)
  
  cov1= model.matrix.lm(~Aspirin_re 
                        + Enrollment_Age + BMI + Education 
                        #+ HHIncome 
                        + family + smk[,2]+smk[,3] + Diabetes + BPH
                        #+ StatinUse 
                        #+ Stroke_TIA 
                        + factor(Enrollment_Year)
                        + aggressive
                        + PSA_Ever
                        + DRE_Ever 
                        + Acetaminophen + nsaids
                        , data=cif,na.action='na.pass')[,-1]
  
  fg_model= crr(cif$VitalStatus_Months,cif$status,cov1=cov1)
  fg_model.sum = summary(fg_model)
  p=round(fg_model.sum$coef[1,"p-value"],digits=3)
  hr=round(fg_model.sum$coef[1,"exp(coef)"],digits=2)
  ci=round(fg_model.sum$conf.int[1,c("2.5%", "97.5%")],digits=2)
  hr.ci.p=paste0(hr,"(",ci[1],"-",ci[2],")"," p=",p)
  return(hr.ci.p)
}



######################## function end ############################
# prostate cancer specific #
temp=data
temp$VitalStatus = ifelse(temp$VitalStatus == "D" & temp$CauseOfDeath=="", "A",temp$VitalStatus )

# variables with mortality #
var.table=rbind(
  cox.OS.var(temp,"Aspirin_any"),
  cox.OS.var(temp,"Aspirin_reg"),
  cox.OS.var(temp,"family"),
  cox.OS.var(temp,"Diabetes"),
  cox.OS.var(temp,"BPH"),
  cox.OS.var(temp,"StatinUse"),
  cox.OS.var(temp,"Stroke_TIA"),
  cox.OS.var(temp,"PSA_Ever"),
  cox.OS.var(temp,"DRE_Ever"),
  cox.OS.var(temp,"Acetaminophen"),
  cox.OS.var(temp,"nsaids"),
  cox.OS.var(temp,"Enrollment_Age"),
  cox.OS.var(temp,"SmokingStatus"),
  cox.OS.var(temp,"BMI"),
  cox.OS.var(temp,"Education"),
  cox.OS.var(temp,"HHIncome"),
  cox.OS.var(temp, "Enrollment_Year")
)
colnames(var.table) = c('var HR p')
#write.csv(var.table,"temp.var.OS.csv",quote = F)


## table of mortality ##
pr.os.table= function(data){
  
temp$Aspirin_re=temp$Aspirin_any
any.t = event.person.mo(temp,"Aspirin_re")
any.os= cox.OS(temp)
any.os.adj= cox.OS.adj(temp)
any.cmp = cmp(temp)
any.cmp.adj=cmp.adj(temp)

temp$Aspirin_re=temp$Aspirin_reg
reg.t = event.person.mo(temp,"Aspirin_re")
reg.os= cox.OS(temp)
reg.os.adj= cox.OS.adj(temp)
reg.cmp = cmp(temp)
reg.cmp.adj=cmp.adj(temp)

os.table=qpcR:::rbind.na(cbind(any.t,c(NA,any.os),c(NA,any.os.adj),c(NA,any.cmp),c(NA,any.cmp.adj)),
                         cbind(reg.t,c(NA,reg.os),c(NA,reg.os.adj),c(NA,reg.cmp),c(NA,reg.cmp.adj)))
colnames(os.table)=c('Events',	'PY',	'HR (95% CI)	P value',	'adj HR (95% CI)	P value','CMP HR (95% CI)	P value',	'CMP adj HR (95% CI)	P value')
return(os.table)
}
pr.os.table.age = function(data){
  
  temp$Aspirin_re=temp$Aspirin_any
  any.t = event.person.mo(temp,"Aspirin_re")
  any.os= cox.OS.age(temp)
  any.os.adj= cox.OS.adj.age(temp)
  any.cmp = cmp(temp)
  any.cmp.adj=cmp.adj(temp)
  
  temp$Aspirin_re=temp$Aspirin_reg
  reg.t = event.person.mo(temp,"Aspirin_re")
  reg.os= cox.OS.age(temp)
  reg.os.adj= cox.OS.adj.age(temp)
  reg.cmp = cmp(temp)
  reg.cmp.adj=cmp.adj(temp)
  
  os.table=qpcR:::rbind.na(cbind(any.t,c(NA,any.os),c(NA,any.os.adj),c(NA,any.cmp),c(NA,any.cmp.adj)),
                           cbind(reg.t,c(NA,reg.os),c(NA,reg.os.adj),c(NA,reg.cmp),c(NA,reg.cmp.adj)))
  colnames(os.table)=c('Events',	'PY',	'HR (95% CI)	P value',	'adj HR (95% CI)	P value','CMP HR (95% CI)	P value',	'CMP adj HR (95% CI)	P value')
  return(os.table)
}





yrs=2015
table2015=pr.os.table(temp)
table2015.age= pr.os.table.age(temp)
yrs=2016
table2016=pr.os.table(temp)
table2016.age= pr.os.table.age(temp)
yrs=2017
table2017=pr.os.table(temp)
table2017.age= pr.os.table.age(temp)
yrs=2018
table2018=pr.os.table(temp)
table2018.age= pr.os.table.age(temp)

all.table=rbind(table2015,table2016,table2017,table2018,
                NA,
                table2015.age,table2016.age,table2017.age,table2018.age)

write.csv(all.table,'table3.csv')
















# high stage #
temp=subset(data,pr==0 | SCCS_SummAJCCStage %in% c("III","IV"))
km.OS(temp,"High stage")

# high Gleason #
temp=subset(data,pr==0 | Earliest_Gleason_Score %in% c("008","009","010"))
km.OS(temp,"High Gleason")




km.OS=function(survival_data,title){
  survival_data=subset(survival_data,time > 0 )
  survival_data$vitalstatus_new=ifelse(survival_data$Year_Death > 2015 | is.na(survival_data$Year_Death),"A",survival_data$VitalStatus)
  survival_data$death=as.numeric(as.factor(survival_data$vitalstatus_new))-1
  fit<- survfit(Surv(VitalStatus_Months, death) ~ Aspirin_re, data = survival_data)
  # Drawing survival curves
  ggsurv=ggsurvplot(fit, data=survival_data,pval = F, conf.int = T,risk.table = T,title = title,ylim=c(0,0.01),xlab="Time(months)",fun="event")
  ggsurv$plot <- ggsurv$plot + 
    ggplot2::annotate("text", 
                      x = 50, y = 0.98, # x and y coordinates of the text
                      label = paste0("HR=",cox.OS(survival_data)), size = 5)
  ggsurv
  return(ggsurv)
}
temp=data
km.OS(temp,"")

Fig1=km.OS(temp,"Prostate Specfic Death")
pdf("Fig1.pdf")
print(Fig1, newpage = FALSE)
dev.off()

################################################################################################################
vars=c("Aspirin_re","Enrollment_Age","Education","HHIncome","BMI","SmokingStatus",
       "family","Diabetes","BPH","StatinUse","Stroke_TIA","Enrollment_Year","PSA_Ever","DRE_Ever","Acetaminophen","nsaids")

+ Enrollment_Age 
+ BMI + Education 
#+ HHIncome 
+ family 
+ smk[,2]+smk[,3] 
+ Diabetes + BPH
#+ StatinUse 
#+ Stroke_TIA 
+ factor(Enrollment_Year)
+ aggressive
+ PSA_Ever 
+ DRE_Ever + Acetaminophen + nsaids



for (i in vars) {
vars.sum = cox.OS.var(temp,i)
print(c(vars.sum,i))
}
 
for (i in vars) {
vars.sum = cox.OS.var(data,i)
print(c(vars.sum,i))
}
 
event.person.mo(temp,"Aspirin_re")
cox.OS(temp)
cox.OS.adj(temp)

non_scr=subset(temp, PSA_Ever ==0)
event.person.mo(non_scr,"Aspirin_re")
cox.OS(non_scr)
cox.OS.adj(non_scr)

scr=subset(temp, PSA_Ever == 1)
event.person.mo(scr,"Aspirin_re")
cox.OS(scr)
cox.OS.adj(scr)


### competing risk ####

### CIF ###

cmp.adj(data)
