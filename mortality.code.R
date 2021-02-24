# source('code.R')


library(survival)
library(knitr)
library(haven)
library(forestplot)

## data import ##
SASData <- read_sas("aspirin_prostatecan_apr2019.sas7bdat")
write.csv(SASData,"aspirin_prostatecan.csv") 

d=SASData
#d=read.csv("aspirin_prostatecan.csv",header=T)

dat=d
dat[dat == "9999"]=NA
dat[dat == "7777"]=NA
dat[dat == "8888"]=NA
 

###################################################################################################
data=data.frame(dat)
 
## total 22427 enrollment, Prostate_Cancer AgeDX , Cancer Registry 929 ##

data$pr=ifelse(data$Prostate_Cancer==1,"PrCa",NA)
data$pr[is.na(data$pr)]="NoPrCa"
data$pr=as.numeric(as.factor(data$pr))-1

## drop 12 month to reflect the difference of case and control in 2016 ##
data$t1=ifelse(is.na(data$Year_Death) | data$Year_Death < 2016, as.numeric(data$VitalStatus_Months)-12,as.numeric(data$VitalStatus_Months))
data$t2=12*(as.numeric(data$AgeDX)-as.numeric(data$Enrollment_Age))
data$time=ifelse(!is.na(data$t2),data$t2,data$t1)

#write.csv(data,"aspirin_prostatecan_cl1.csv") 
### construct case and control characteristics ###
n=nrow(data) #n=22427

# combine 2 variables #
variable_combine=function(x,y){
x[is.na(x)]= -9
y[is.na(y)]= -9
z=ifelse( (x+y) %in% c(-8,1,2), 1,ifelse( (x+y) %in% c(-9,0),0,NA))
return(z)
}

# combine fa and br prostate history #
history_combine=function(x,y){
x[is.na(x)]=0
y[is.na(y)]=0
z=ifelse(x+y >0, 1,0)
return(z)
}
data$Fa_ProstateCancer=history_combine(data$Br_ProstateCancer,data$Fa_ProstateCancer)
data$Fa_ProstateCancerF1=history_combine(data$Br_ProstateCancerF1,data$Fa_ProstateCancerF1)

# combine OTCNSAIDS and RxNSAIDS as nsaids #
data$nsaids=history_combine(data$RxNSAIDS,data$OTCNSAIDS)


# re-define Asprin use #
data$Aspirin_re = ifelse(data$Aspirin== 0 & data$LowAspirin ==1, NA ,data$Aspirin )
data$AspirinF1_re = ifelse(data$AspirinF1 == 0 & data$LowAspirinF1 ==1, NA ,data$AspirinF1 )

cases=subset(data,pr==1)
controls=subset(data,pr==0)

###################################################################################################
### table1 ###
# case - control #
cases.n=nrow(cases)
controls.n=nrow(controls)
cc=t(c(cases.n,100*cases.n/n, controls.n, 100*controls.n/n))
rownames(cc)="Group"
# Enrollment_Age #
age=t(c(median(cases$Enrollment_Age),IQR(cases$Enrollment_Age),median(controls$Enrollment_Age),IQR(controls$Enrollment_Age)))
rownames(age)="Age median IQR"
# BMI #
bmi=t(c(mean(cases$BMI,na.rm=T),sd(cases$BMI,na.rm=T),mean(controls$BMI,na.rm=T),sd(controls$BMI,na.rm=T)))
rownames(bmi)="BMI mean SD"
# Education ,Less High, High school, college, afterCollege #
edu=cbind(table(transform(cases, bin = cut(Education, breaks=c(1,3,5,6,8),right=F))$bin),
		100*table(transform(cases, bin = cut(Education, breaks=c(1,3,5,6,8),right=F))$bin)/cases.n,
	    table(transform(controls, bin = cut(Education, breaks=c(1,3,5,6,8),right=F))$bin),
		100*table(transform(controls, bin = cut(Education, breaks=c(1,3,5,6,8),right=F))$bin)/controls.n
	)
rownames(edu)=c("Education Less High", "High school", "College", "AfterCollege")
# HHIncome  <$15000, $15000-25000,$25000-50000,$50000-100000, $100,000 #
income=cbind(table(cases$HHIncome),
		100*table(cases$HHIncome)/cases.n,
		table(controls$HHIncome),
		100*table(controls$HHIncome)/controls.n
		)
rownames(income)=c("<$15,000","$15,000-25,000","$25,000-50,000","$50,000-100,000",">$100,000")

## Fa_ProstateCancer ##
fa_history=cbind(table(cases$Fa_ProstateCancer)["1"],
		100*table(cases$Fa_ProstateCancer)["1"]/cases.n,
		table(controls$Fa_ProstateCancer)["1"],
		100*table(controls$Fa_ProstateCancer)["1"]/controls.n
		)
rownames(fa_history)="Family History of Prostate Cancer Yes"

##aspirin use ##
asp=cbind(table(cases$Aspirin_re)[c("0","1")],
		100*table(cases$Aspirin_re)[c("0","1")]/cases.n,
		table(controls$Aspirin_re)[c("0","1")],
		100*table(controls$Aspirin_re)[c("0","1")]/controls.n 
		)
rownames(asp)=c("Aspirin No","Yes")

# SmokingStatus current former never #
#1=Current
#2=Former
#3=Never
smk=cbind(table(cases$SmokingStatus)[c("1","2","3")],
		100*table(cases$SmokingStatus)[c("1","2","3")]/cases.n,
		table(controls$SmokingStatus)[c("1","2","3")],
		100*table(controls$SmokingStatus)[c("1","2","3")]/controls.n
		)
rownames(smk)=c("Smoke Current","Former","Never")
# Diabetes #
diabetes=cbind(table(cases$Diabetes)[c("0","1")],
		100*table(cases$Diabetes)[c("0","1")]/cases.n,
		table(controls$Diabetes)[c("0","1")],
		100*table(controls$Diabetes)[c("0","1")]/controls.n 
		)
rownames(diabetes)=c("Diabetes No","Yes")
# BPH #
bph=cbind(table(cases$BPH)[c("0","1")],
		100*table(cases$BPH)[c("0","1")]/cases.n,
		table(controls$BPH)[c("0","1")],
		100*table(controls$BPH)[c("0","1")]/controls.n
		)
rownames(bph)=c("BPH No","Yes")
# SCCS_SummAJCCStage #
stage=cbind(table(cases$SCCS_SummAJCCStage)[c("I","II","III","IV")],
		100*table(cases$SCCS_SummAJCCStage)[c("I","II","III","IV")]/cases.n
		)
rownames(stage)=c("Stage I","II","III","IV")
# Earliest_Gleason_Score, 000 test not done #
gl=   rbind(cbind(sum(table(cases$Earliest_Gleason_Score)[c("002","003","004","005","006","007")]),
		100*sum(table(cases$Earliest_Gleason_Score)[c("002","003","004","005","006","007")])/cases.n),
		cbind(sum(table(cases$Earliest_Gleason_Score)[c("008","009","010")]),
		100*sum(table(cases$Earliest_Gleason_Score)[c("008","009","010")])/cases.n)
		)
rownames(gl)=c("Gleason <=7",">7")
# diseaes aggressiveness gleason >7 and T3T4 #
ag=subset(cases,(Earliest_Gleason_Score %in% c("008","009","010")) & (SCCS_SummAJCCStage %in% c("III","IV")))
non_ag=subset(cases,(Earliest_Gleason_Score %in% c("002","003","004","005","006","007")) & (SCCS_SummAJCCStage %in% c("I","II")))

aggressive=rbind(cbind(nrow(ag),100*nrow(ag)/cases.n),
			cbind(nrow(non_ag),100*nrow(non_ag)/cases.n))
rownames(aggressive)=c("Aggressive Disease","Nonaggressive Disease")


# freePSA_ngml #
freepsa=c("free PSA median IQR",median(cases$freePSA_ngml,na.rm=T),IQR(cases$freePSA_ngml,na.rm=T))
# intactPSA_ngml #
intactpsa=c("intact PSA median IQR",median(cases$intactPSA_ngml,na.rm=T),IQR(cases$intactPSA_ngml,na.rm=T))
# totPSA_ngml #
psa=t(c(median(cases$totPSA_ngml,na.rm=T),IQR(cases$totPSA_ngml,na.rm=T)))
rownames(psa)="total PSA median IQR"


# PSA_Ever #
psa.ever=cbind(table(cases$PSA_Ever)[c("0","1")],
		100*table(cases$PSA_Ever)[c("0","1")]/cases.n,
		table(controls$PSA_Ever)[c("0","1")],
		100*table(controls$PSA_Ever)[c("0","1")]/controls.n 
		)
rownames(psa.ever)=c("PSA_Ever No","Yes")

# DRE_Ever #
dre.ever=cbind(table(cases$DRE_Ever)[c("0","1")],
		100*table(cases$DRE_Ever)[c("0","1")]/cases.n,
		table(controls$DRE_Ever)[c("0","1")],
		100*table(controls$DRE_Ever)[c("0","1")]/controls.n 
		)
rownames(dre.ever)=c("DRE_Ever No","Yes")






###############################################################################################################
###

### overall survival  ###
#switch to higher version R for plotting#
library(survival)
library("survminer")
############################# function  #######################
event.person.mo=function(survival_data,var) {
	
	survival_data$vitalstatus_new=ifelse(survival_data$Year_Death > 2015 | is.na(survival_data$Year_Death),"A",survival_data$VitalStatus)
	survival_data$death=as.numeric(as.factor(survival_data$vitalstatus_new))-1
      tl=table(survival_data$death,survival_data[,var])
	set0=sum(subset(survival_data,   get(var) ==0)$VitalStatus_Months	)
      set1=sum(subset(survival_data,   get(var) ==1)$VitalStatus_Months	)
	sum.tl= rbind(c(tl[2,1],set0),c(tl[2,2],set1))
      colnames(sum.tl)=c("events","person-month")
	rownames(sum.tl)=c("No","Yes")
      return(sum.tl)
			}
event.person.mo(data,"Aspirin_re")


cox.OS.var=function(data,var){
	survival_data=subset(data,time > 0 )
	survival_data$vitalstatus_new=ifelse(survival_data$Year_Death > 2015 | is.na(survival_data$Year_Death),"A",survival_data$VitalStatus)
	survival_data$death=as.numeric(as.factor(survival_data$vitalstatus_new))-1
	cox.sum=summary(coxph(Surv(VitalStatus_Months, death)~ get(var) + Enrollment_Age,data=survival_data))
	p=round(cox.sum$coefficients[1,"Pr(>|z|)"],digits=3)
	hr=round(cox.sum$coefficients[1,"exp(coef)"],digits=2)
	ci=round(cox.sum$conf.int[1,c("lower .95", "upper .95")],digits=2)
	hr.ci.p=paste0(hr,"(",ci[1],"-",ci[2],")"," p=",p)
	return(hr.ci.p)
	#return(cox.sum)
}
cox.OS.var(data,"Aspirin_re")




cox.OS=function(data){
	survival_data=subset(data,time > 0 )
	survival_data$vitalstatus_new=ifelse(survival_data$Year_Death > 2015 | is.na(survival_data$Year_Death),"A",survival_data$VitalStatus)
	survival_data$death=as.numeric(as.factor(survival_data$vitalstatus_new))-1
	cox.sum=summary(coxph(Surv(VitalStatus_Months, death)~Aspirin_re + Enrollment_Age,data=survival_data))
	p=round(cox.sum$coefficients[1,"Pr(>|z|)"],digits=3)
	hr=round(cox.sum$coefficients[1,"exp(coef)"],digits=2)
	ci=round(cox.sum$conf.int[1,c("lower .95", "upper .95")],digits=2)
	hr.ci.p=paste0(hr,"(",ci[1],"-",ci[2],")"," p=",p)
	return(hr.ci.p)
	#return(cox.sum)
}
#cox.OS(data)
cox.OS.adj=function(data){
	survival_data=subset(data,time > 0 )
	survival_data$vitalstatus_new=ifelse(survival_data$Year_Death > 2015 | is.na(survival_data$Year_Death),"A",survival_data$VitalStatus)
	survival_data$death=as.numeric(as.factor(survival_data$vitalstatus_new))-1
	stage=ifelse(survival_data$SCCS_SummAJCCStage %in% c("III","IV"),"High",ifelse(survival_data$SCCS_SummAJCCStage %in% c("I","II"),"Low",NA))
	gl=ifelse(survival_data$Earliest_Gleason_Score %in% c("002","003","004","005","006","007"),"Low",ifelse(survival_data$Earliest_Gleason_Score %in% c("008","009","010"),"High",NA))
	smk <- model.matrix.lm(~factor(SmokingStatus,levels=c(3,2,1)), data = survival_data,na.action="na.pass")
	cox.sum=summary(coxph(Surv(VitalStatus_Months, death)~ Aspirin_re 
	+ Education + HHIncome+ Enrollment_Age+BMI +as.factor(SmokingStatus)
	+ PSA_Ever  + Fa_ProstateCancer
	,data=survival_data))
	p=round(cox.sum$coefficients[1,"Pr(>|z|)"],digits=3)
	hr=round(cox.sum$coefficients[1,"exp(coef)"],digits=2)
	ci=round(cox.sum$conf.int[1,c("lower .95", "upper .95")],digits=2)
	hr.ci.p=paste0(hr,"(",ci[1],"-",ci[2],")"," p=",p)
	return(hr.ci.p)
	#return(cox.sum)
}
cox.OS.adj(temp)




km.OS=function(survival_data,title){
	survival_data=subset(survival_data,time > 0 )
	survival_data$vitalstatus_new=ifelse(survival_data$Year_Death > 2015 | is.na(survival_data$Year_Death),"A",survival_data$VitalStatus)
	survival_data$death=as.numeric(as.factor(survival_data$vitalstatus_new))-1
	fit<- survfit(Surv(VitalStatus_Months, death) ~ Aspirin_re, data = survival_data)
      # Drawing survival curves
	ggsurv=ggsurvplot(fit, data=survival_data,pval = F, conf.int = T,risk.table = T,title = title,ylim=c(0,0.01),xlab="Time(months)",fun="event")
	ggsurv$plot <- ggsurv$plot + 
              ggplot2::annotate("text", 
                                x = 50, y = 0.0098, # x and y coordinates of the text
                                label = paste0("HR=",cox.OS(survival_data)), size = 5)
	ggsurv
	return(ggsurv)
}
temp=data
km.OS(temp,"")
######################## function end ############################

# low stage #
temp=subset(data,pr==0 | SCCS_SummAJCCStage %in% c("I","II"))
km.OS(temp,"Low stage")

# high stage #
temp=subset(data,pr==0 | SCCS_SummAJCCStage %in% c("III","IV"))
km.OS(temp,"High stage")

# low Gleason  #
temp=subset(data,pr==0 | Earliest_Gleason_Score %in% c("002","003","004","005","006","007"))
km.OS(temp,"Low Gleason")

# high Gleason #
temp=subset(data,pr==0 | Earliest_Gleason_Score %in% c("008","009","010"))
km.OS(temp,"High Gleason")

# prostate cancer specific #
temp=data
temp$VitalStatus = ifelse(temp$VitalStatus == "D" & temp$CauseOfDeath=="", "A",temp$VitalStatus )
Fig1=km.OS(temp,"Prostate Specfic Death")
pdf("Fig1.pdf")
print(Fig1, newpage = FALSE)
dev.off()

################################################################################################################
vars=c("Aspirin_re","Education","HHIncome","BMI","SmokingStatus","Fa_ProstateCancer","Diabetes","BPH","PSA_Ever","DRE_Ever")
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

















