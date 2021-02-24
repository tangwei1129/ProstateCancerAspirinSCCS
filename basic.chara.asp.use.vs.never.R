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
# calculate missing % #
library(skimr)
write.csv(skim(dat),'missing.csv')

###################################################################################################
data=data.frame(dat)

n=nrow(data) #n=22426


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


###################################################################################################################

### set aspirin_re as any or reg ###
data$Aspirin_re= data$Aspirin_any
#data$Aspirin_re= data$Aspirin_reg

###################################################################################################################
# aspirin use  #
data=subset(data,!is.na(Aspirin_re))
n=nrow(data)#21851
asp.use=subset(data,Aspirin_re == 1)
asp.never=subset(data,Aspirin_re == 0)
asp.use.n=nrow(asp.use)##5486
asp.never.n=nrow(asp.never)#16365
cc=t(c(nrow(asp.use),100*nrow(asp.use)/n, nrow(asp.never), 100*nrow(asp.never)/n))
rownames(cc)="Group"

# Enrollment_Age #
age=t(c(median(asp.use$Enrollment_Age),IQR(asp.use$Enrollment_Age),median(asp.never$Enrollment_Age),IQR(asp.never$Enrollment_Age),
        t.test(asp.use$Enrollment_Age,asp.never$Enrollment_Age)$p.value))
rownames(age)="Age median IQR"
# BMI #
bmi=t(c(median(asp.use$BMI,na.rm=T),IQR(asp.use$BMI,na.rm=T),median(asp.never$BMI,na.rm=T),IQR(asp.never$BMI,na.rm=T),
        t.test(asp.use$BMI,asp.never$BMI)$p.value))
rownames(bmi)="BMI mean IQR"
# Education ,Less High, High school, college, afterCollege #

edu=cbind(table(transform(asp.use, bin = cut(Education, breaks=c(1,3,5,6,8),right=F))$bin),
          100*table(transform(asp.use, bin = cut(Education, breaks=c(1,3,5,6,8),right=F))$bin)/asp.use.n,
          table(transform(asp.never, bin = cut(Education, breaks=c(1,3,5,6,8),right=F))$bin),
          100*table(transform(asp.never, bin = cut(Education, breaks=c(1,3,5,6,8),right=F))$bin)/asp.never.n
)
rownames(edu)=c("Education Less High", " High school", " College", " AfterCollege")
edu=cbind(edu,c(NA,NA,NA,chisq.test(edu[,c(1,3)])$p.value))

# HHIcome "Less than 15k","15K-25K","25K-50","50K-100K","More than 100k" ##
hhincome=cbind(table(asp.use$HHIncome)[c("1","2","3","4","5")],
          100*table(asp.use$HHIncome)[c("1","2","3","4","5")]/asp.use.n,
          table(asp.never$HHIncome)[c("1","2","3","4","5")],
          100*table(asp.never$HHIncome)[c("1","2","3","4","5")]/asp.never.n
)
rownames(hhincome)=c("Less than 15k","15K-25K","25K-50","50K-100K","More than 100k")
hhincome=cbind(hhincome,c(NA,NA,NA,NA,chisq.test(hhincome[,c(1,3)])$p.value))

## Fa_ProstateCancer ##
fa_history=cbind(table(asp.use$family)[c("0","1")],
                 100*table(asp.use$family)[c("0","1")]/asp.use.n,
                 table(asp.never$family)[c("0","1")],
                 100*table(asp.never$family)[c("0","1")]/asp.never.n 
)
rownames(fa_history)=c("Family history No"," Yes")
fa_history=cbind(fa_history,c(NA,chisq.test(fa_history[,c(1,3)])$p.value))


# SmokingStatus current former never #
#1=Current
#2=Former
#3=Never
smk=cbind(table(asp.use$SmokingStatus)[c("1","2","3")],
          100*table(asp.use$SmokingStatus)[c("1","2","3")]/asp.use.n,
          table(asp.never$SmokingStatus)[c("1","2","3")],
          100*table(asp.never$SmokingStatus)[c("1","2","3")]/asp.never.n
)
rownames(smk)=c("Smoke Current"," Former"," Never")
smk=cbind(smk,c(NA,NA,chisq.test(smk[,c(1,3)])$p.value))



# Diabetes #
diabetes=cbind(table(asp.use$Diabetes)[c("0","1")],
               100*table(asp.use$Diabetes)[c("0","1")]/asp.use.n,
               table(asp.never$Diabetes)[c("0","1")],
               100*table(asp.never$Diabetes)[c("0","1")]/asp.never.n 
)
rownames(diabetes)=c("Diabetes No"," Yes")
diabetes=cbind(diabetes,c(NA,chisq.test(diabetes[,c(1,3)])$p.value))

# BPH #
bph=cbind(table(asp.use$BPH)[c("0","1")],
          100*table(asp.use$BPH)[c("0","1")]/asp.use.n,
          table(asp.never$BPH)[c("0","1")],
          100*table(asp.never$BPH)[c("0","1")]/asp.never.n
)
rownames(bph)=c("BPH No"," Yes")
bph=cbind(bph,c(NA,chisq.test(bph[,c(1,3)])$p.value))

# StatinUse #
statin=cbind(table(asp.use$StatinUse)[c("0","1")],
          100*table(asp.use$StatinUse)[c("0","1")]/asp.use.n,
          table(asp.never$StatinUse)[c("0","1")],
          100*table(asp.never$StatinUse)[c("0","1")]/asp.never.n
          )
rownames(statin)=c("StatinUse No"," Yes")
statin=cbind(statin,c(NA,chisq.test(statin[,c(1,3)])$p.value))

# Stroke_TIA #
stroke=cbind(table(asp.use$Stroke_TIA)[c("0","1")],
          100*table(asp.use$Stroke_TIA)[c("0","1")]/asp.use.n,
          table(asp.never$Stroke_TIA)[c("0","1")],
          100*table(asp.never$Stroke_TIA)[c("0","1")]/asp.never.n
)
rownames(stroke)=c("Stroke No"," Yes")
stroke=cbind(stroke,c(NA,chisq.test(stroke[,c(1,3)])$p.value))

# PSA screen #
psa.scr =cbind(table(asp.use$PSA_Ever)[c("0","1")],
          100*table(asp.use$PSA_Ever)[c("0","1")]/asp.use.n,
          table(asp.never$PSA_Ever)[c("0","1")],
          100*table(asp.never$PSA_Ever)[c("0","1")]/asp.never.n
)
rownames(psa.scr)=c("PSA No"," Yes")
psa.scr=cbind(psa.scr,c(NA,chisq.test(psa.scr[,c(1,3)])$p.value))

# DRE screen #
dre.scr =cbind(table(asp.use$DRE_Ever)[c("0","1")],
               100*table(asp.use$DRE_Ever)[c("0","1")]/asp.use.n,
               table(asp.never$DRE_Ever)[c("0","1")],
               100*table(asp.never$DRE_Ever)[c("0","1")]/asp.never.n
)
rownames(dre.scr)=c("DRE No"," Yes")
dre.scr=cbind(dre.scr,c(NA,chisq.test(dre.scr[,c(1,3)])$p.value))

# SCCS_SummAJCCStage #
asp.use.n.case=table(data$Aspirin_re,data$Prostate_Cancer)[2,1]
asp.never.n.case = table(data$Aspirin_re,data$Prostate_Cancer)[1,1]

stage=cbind(table(asp.use$SCCS_SummAJCCStage)[c("I","II","III","IV")],
            100*table(asp.use$SCCS_SummAJCCStage)[c("I","II","III","IV")]/asp.use.n.case,
            table(asp.never$SCCS_SummAJCCStage)[c("I","II","III","IV")],
            100*table(asp.never$SCCS_SummAJCCStage)[c("I","II","III","IV")]/asp.never.n.case
)
rownames(stage)=c("Stage I"," II"," III"," IV")
stage=cbind(stage,c(NA,NA,NA,chisq.test(stage[,c(1,3)])$p.value))

# Earliest_Gleason_Score, 000 test not done #
lowG = c("002","003","004","005","006","007")
highG = c("008","009","010")
gl=   cbind(rbind(cbind(table(asp.use$Earliest_Gleason_Score %in% lowG)[2],
                        100*table(asp.use$Earliest_Gleason_Score %in% lowG)[2]/asp.use.n.case),
                  cbind(table(asp.use$Earliest_Gleason_Score %in% highG)[2],
                        100*table(asp.use$Earliest_Gleason_Score %in% highG)[2]/asp.use.n.case)
),
rbind(cbind(table(asp.never$Earliest_Gleason_Score %in% lowG)[2],
            100*table(asp.never$Earliest_Gleason_Score %in% lowG)[2]/asp.never.n.case),
      cbind(table(asp.never$Earliest_Gleason_Score %in% highG)[2],
            100*table(asp.never$Earliest_Gleason_Score %in% highG)[2]/asp.never.n.case)
)
)

rownames(gl)=c("Gleason <=7"," >=8")
gl=cbind(gl,c(NA,chisq.test(gl[,c(1,3)])$p.value))

# diseaes aggressiveness gleason >7 and T3T4 as stage T4 or N1 or M1 or Gleason score >=8 ##
ag.filter = "Earliest_Gleason_Score %in% highG | SCCS_SummAJCCStage == 'IV'  |
              TNM_Clin_M %in% c('1','1B','1C') |	TNM_Path_M %in% c('1A','1B') |
              TNM_Clin_N == '1' | 	TNM_Path_N == '1' |	
              TNM_Clin_T == '4' |   TNM_Path_T == '4' "
            
ag=subset(asp.use, Earliest_Gleason_Score %in% highG | SCCS_SummAJCCStage == "IV"  |
            TNM_Clin_M %in% c('1','1B','1C') |	TNM_Path_M %in% c('1A','1B') |
            TNM_Clin_N == '1' | 	TNM_Path_N == '1' |	
            TNM_Clin_T == '4' |   TNM_Path_T == '4')

ag.asp.never=subset(asp.never,Earliest_Gleason_Score %in% highG | SCCS_SummAJCCStage == "IV"  |
                      TNM_Clin_M %in% c('1','1B','1C') |	TNM_Path_M %in% c('1A','1B') |
                      TNM_Clin_N == '1' | 	TNM_Path_N == '1' |	
                      TNM_Clin_T == '4' |   TNM_Path_T == '4')

aggressive=cbind(rbind(cbind(nrow(ag),100*nrow(ag)/asp.use.n.case),
                       cbind(asp.use.n.case - nrow(ag),100-100*nrow(ag)/asp.use.n.case)),
                 rbind(cbind(nrow(ag.asp.never),100*nrow(ag.asp.never)/asp.never.n.case),
                       cbind(asp.never.n.case - nrow(ag.asp.never),100-100*nrow(ag)/asp.never.n.case)
                  ))


rownames(aggressive)=c("Aggressive Disease","Nonaggressive Disease")
aggressive=cbind(aggressive,c(NA,chisq.test(aggressive[,c(1,3)])$p.value))


# BMI record #

bmi_re=cbind(table(transform(asp.use, bin = cut(BMI, breaks=c(0,18.5,25,30,35,40,100),right=F))$bin),
          100*table(transform(asp.use, bin = cut(BMI, breaks=c(0,18.5,25,30,35,40,100),right=F))$bin)/asp.use.n,
          table(transform(asp.never, bin = cut(BMI, breaks=c(0,18.5,25,30,35,40,100),right=F))$bin),
          100*table(transform(asp.never, bin = cut(BMI, breaks=c(0,18.5,25,30,35,40,100),right=F))$bin)/asp.never.n
)
rownames(bmi_re)=c("BMI <18.5", " [18.5-25)", " [25-30)"," [30-35)", " [35-40)", " >40")
bmi_re=cbind(bmi_re,c(NA,NA,NA,NA,NA,chisq.test(bmi_re[,c(1,3)])$p.value))


## sum up ##
table.asp.use.asp.never= qpcR:::rbind.na(cc,age,bmi,edu,hhincome,fa_history,smk,diabetes,bph,statin,stroke,psa.scr,dre.scr,stage,gl,aggressive,bmi_re)
colnames(table.asp.use.asp.never)=c("asp.use, N","asp.use, %","asp.never, N","asp.never, %","P value")
write.csv(table.asp.use.asp.never,"table1.any.asp.useVSasp.never.csv")





#####################################################################################################
########################################################################################################


#data$Aspirin_re= data$Aspirin_reg


data=data.frame(dat)

n=nrow(data) #n=22427
 
data$family=history_combine(data$Fa_ProstateCancer,data$Br_ProstateCancer,data$Fa_Cancer)
data$nsaids=nsaid_combine(data$OTCNSAIDS,data$RxNSAIDS)

## any aspirin ##
data$Aspirin_any = variable_combine(data$Aspirin, data$LowAspirin )
##regular aspirin ##
data$Aspirin_reg = ifelse(!data$Aspirin == 1 & data$Aspirin_any ==1, NA, data$Aspirin_any )


###################################################################################################################

### set aspirin_re as any or reg ###
#data$Aspirin_re= data$Aspirin_any
data$Aspirin_re= data$Aspirin_reg

###################################################################################################################
# aspirin use  #
data=subset(data,!is.na(Aspirin_re))
n=nrow(data)
asp.use=subset(data,Aspirin_re == 1)
asp.never=subset(data,Aspirin_re == 0)
asp.use.n=nrow(asp.use)
asp.never.n=nrow(asp.never)
cc=t(c(nrow(asp.use),100*nrow(asp.use)/21851, nrow(asp.never), 100*nrow(asp.never)/21851))
rownames(cc)="Group"


# Enrollment_Age #
age=t(c(median(asp.use$Enrollment_Age),IQR(asp.use$Enrollment_Age),median(asp.never$Enrollment_Age),IQR(asp.never$Enrollment_Age),
        t.test(asp.use$Enrollment_Age,asp.never$Enrollment_Age)$p.value))
rownames(age)="Age median IQR"
# BMI #
bmi=t(c(median(asp.use$BMI,na.rm=T),IQR(asp.use$BMI,na.rm=T),median(asp.never$BMI,na.rm=T),IQR(asp.never$BMI,na.rm=T),
        t.test(asp.use$BMI,asp.never$BMI)$p.value))
rownames(bmi)="BMI median IQR"
# Education ,Less High, High school, college, afterCollege #

edu=cbind(table(transform(asp.use, bin = cut(Education, breaks=c(1,3,5,6,8),right=F))$bin),
          100*table(transform(asp.use, bin = cut(Education, breaks=c(1,3,5,6,8),right=F))$bin)/asp.use.n,
          table(transform(asp.never, bin = cut(Education, breaks=c(1,3,5,6,8),right=F))$bin),
          100*table(transform(asp.never, bin = cut(Education, breaks=c(1,3,5,6,8),right=F))$bin)/asp.never.n
)
rownames(edu)=c("Education Less High", " High school", " College", " AfterCollege")
edu=cbind(edu,c(NA,NA,NA,chisq.test(edu[,c(1,3)])$p.value))

# HHIcome "Less than 15k","15K-25K","25K-50","50K-100K","More than 100k" ##
hhincome=cbind(table(asp.use$HHIncome)[c("1","2","3","4","5")],
               100*table(asp.use$HHIncome)[c("1","2","3","4","5")]/asp.use.n,
               table(asp.never$HHIncome)[c("1","2","3","4","5")],
               100*table(asp.never$HHIncome)[c("1","2","3","4","5")]/asp.never.n
)
rownames(hhincome)=c("Less than 15k","15K-25K","25K-50","50K-100K","More than 100k")
hhincome=cbind(hhincome,c(NA,NA,NA,NA,chisq.test(hhincome[,c(1,3)])$p.value))

## Fa_ProstateCancer ##
fa_history=cbind(table(asp.use$family)[c("0","1")],
               100*table(asp.use$family)[c("0","1")]/asp.use.n,
               table(asp.never$family)[c("0","1")],
               100*table(asp.never$family)[c("0","1")]/asp.never.n 
)
rownames(fa_history)=c("Family history No"," Yes")
fa_history=cbind(fa_history,c(NA,chisq.test(fa_history[,c(1,3)])$p.value))


# SmokingStatus current former never #
#1=Current
#2=Former
#3=Never
smk=cbind(table(asp.use$SmokingStatus)[c("1","2","3")],
          100*table(asp.use$SmokingStatus)[c("1","2","3")]/asp.use.n,
          table(asp.never$SmokingStatus)[c("1","2","3")],
          100*table(asp.never$SmokingStatus)[c("1","2","3")]/asp.never.n
)
rownames(smk)=c("Smoke Current"," Former"," Never")
smk=cbind(smk,c(NA,NA,chisq.test(smk[,c(1,3)])$p.value))



# Diabetes #
diabetes=cbind(table(asp.use$Diabetes)[c("0","1")],
               100*table(asp.use$Diabetes)[c("0","1")]/asp.use.n,
               table(asp.never$Diabetes)[c("0","1")],
               100*table(asp.never$Diabetes)[c("0","1")]/asp.never.n 
)
rownames(diabetes)=c("Diabetes No"," Yes")
diabetes=cbind(diabetes,c(NA,chisq.test(diabetes[,c(1,3)])$p.value))

# BPH #
bph=cbind(table(asp.use$BPH)[c("0","1")],
          100*table(asp.use$BPH)[c("0","1")]/asp.use.n,
          table(asp.never$BPH)[c("0","1")],
          100*table(asp.never$BPH)[c("0","1")]/asp.never.n
)
rownames(bph)=c("BPH No"," Yes")
bph=cbind(bph,c(NA,chisq.test(bph[,c(1,3)])$p.value))

# StatinUse #
statin=cbind(table(asp.use$StatinUse)[c("0","1")],
             100*table(asp.use$StatinUse)[c("0","1")]/asp.use.n,
             table(asp.never$StatinUse)[c("0","1")],
             100*table(asp.never$StatinUse)[c("0","1")]/asp.never.n
)
rownames(statin)=c("StatinUse No"," Yes")
statin=cbind(statin,c(NA,chisq.test(statin[,c(1,3)])$p.value))

# Stroke_TIA #
stroke=cbind(table(asp.use$Stroke_TIA)[c("0","1")],
             100*table(asp.use$Stroke_TIA)[c("0","1")]/asp.use.n,
             table(asp.never$Stroke_TIA)[c("0","1")],
             100*table(asp.never$Stroke_TIA)[c("0","1")]/asp.never.n
)
rownames(stroke)=c("Stroke No"," Yes")
stroke=cbind(stroke,c(NA,chisq.test(stroke[,c(1,3)])$p.value))

# PSA screen #
psa.scr =cbind(table(asp.use$PSA_Ever)[c("0","1")],
               100*table(asp.use$PSA_Ever)[c("0","1")]/asp.use.n,
               table(asp.never$PSA_Ever)[c("0","1")],
               100*table(asp.never$PSA_Ever)[c("0","1")]/asp.never.n
)
rownames(psa.scr)=c("PSA No"," Yes")
psa.scr=cbind(psa.scr,c(NA,chisq.test(psa.scr[,c(1,3)])$p.value))

# DRE screen #
dre.scr =cbind(table(asp.use$DRE_Ever)[c("0","1")],
               100*table(asp.use$DRE_Ever)[c("0","1")]/asp.use.n,
               table(asp.never$DRE_Ever)[c("0","1")],
               100*table(asp.never$DRE_Ever)[c("0","1")]/asp.never.n
)
rownames(dre.scr)=c("DRE No"," Yes")
dre.scr=cbind(dre.scr,c(NA,chisq.test(dre.scr[,c(1,3)])$p.value))



# SCCS_SummAJCCStage #
asp.use.n.case=table(data$Aspirin_re,data$Prostate_Cancer)[2,1]
asp.never.n.case = table(data$Aspirin_re,data$Prostate_Cancer)[1,1]

stage=cbind(table(asp.use$SCCS_SummAJCCStage)[c("I","II","III","IV")],
            100*table(asp.use$SCCS_SummAJCCStage)[c("I","II","III","IV")]/asp.use.n.case,
            table(asp.never$SCCS_SummAJCCStage)[c("I","II","III","IV")],
            100*table(asp.never$SCCS_SummAJCCStage)[c("I","II","III","IV")]/asp.never.n.case
)
rownames(stage)=c("Stage I"," II"," III"," IV")
stage=cbind(stage,c(NA,NA,NA,chisq.test(stage[,c(1,3)])$p.value))

# Earliest_Gleason_Score, 000 test not done #
lowG = c("002","003","004","005","006","007")
highG = c("008","009","010")
gl=   cbind(rbind(cbind(table(asp.use$Earliest_Gleason_Score %in% lowG)[2],
                        100*table(asp.use$Earliest_Gleason_Score %in% lowG)[2]/asp.use.n.case),
                  cbind(table(asp.use$Earliest_Gleason_Score %in% highG)[2],
                        100*table(asp.use$Earliest_Gleason_Score %in% highG)[2]/asp.use.n.case)
),
rbind(cbind(table(asp.never$Earliest_Gleason_Score %in% lowG)[2],
            100*table(asp.never$Earliest_Gleason_Score %in% lowG)[2]/asp.never.n.case),
      cbind(table(asp.never$Earliest_Gleason_Score %in% highG)[2],
            100*table(asp.never$Earliest_Gleason_Score %in% highG)[2]/asp.never.n.case)
)
)

rownames(gl)=c("Gleason <=7"," >=8")
gl=cbind(gl,c(NA,chisq.test(gl[,c(1,3)])$p.value))

# diseaes aggressiveness gleason >7 and T3T4 as stage T4 or N1 or M1 or Gleason score >=8 ##
ag.filter = "Earliest_Gleason_Score %in% highG | SCCS_SummAJCCStage == 'IV'  |
              TNM_Clin_M %in% c('1','1B','1C') |	TNM_Path_M %in% c('1A','1B') |
              TNM_Clin_N == '1' | 	TNM_Path_N == '1' |	
              TNM_Clin_T == '4' |   TNM_Path_T == '4' "

ag=subset(asp.use, Earliest_Gleason_Score %in% highG | SCCS_SummAJCCStage == "IV"  |
            TNM_Clin_M %in% c('1','1B','1C') |	TNM_Path_M %in% c('1A','1B') |
            TNM_Clin_N == '1' | 	TNM_Path_N == '1' |	
            TNM_Clin_T == '4' |   TNM_Path_T == '4')

ag.asp.never=subset(asp.never,Earliest_Gleason_Score %in% highG | SCCS_SummAJCCStage == "IV"  |
                      TNM_Clin_M %in% c('1','1B','1C') |	TNM_Path_M %in% c('1A','1B') |
                      TNM_Clin_N == '1' | 	TNM_Path_N == '1' |	
                      TNM_Clin_T == '4' |   TNM_Path_T == '4')

aggressive=cbind(rbind(cbind(nrow(ag),100*nrow(ag)/asp.use.n.case),
                       cbind(asp.use.n.case - nrow(ag),100-100*nrow(ag)/asp.use.n.case)),
                 rbind(cbind(nrow(ag.asp.never),100*nrow(ag.asp.never)/asp.never.n.case),
                       cbind(asp.never.n.case - nrow(ag.asp.never),100-100*nrow(ag)/asp.never.n.case)
                 ))


rownames(aggressive)=c("Aggressive Disease","Nonaggressive Disease")
aggressive=cbind(aggressive,c(NA,chisq.test(aggressive[,c(1,3)])$p.value))

# BMI record #

bmi_re=cbind(table(transform(asp.use, bin = cut(BMI, breaks=c(0,18.5,25,30,35,40,100),right=F))$bin),
             100*table(transform(asp.use, bin = cut(BMI, breaks=c(0,18.5,25,30,35,40,100),right=F))$bin)/asp.use.n,
             table(transform(asp.never, bin = cut(BMI, breaks=c(0,18.5,25,30,35,40,100),right=F))$bin),
             100*table(transform(asp.never, bin = cut(BMI, breaks=c(0,18.5,25,30,35,40,100),right=F))$bin)/asp.never.n
)
rownames(bmi_re)=c("BMI <18.5", " [18.5-25)", " [25-30)"," [30-35)", " [35-40)", " >40")
bmi_re=cbind(bmi_re,c(NA,NA,NA,NA,NA,chisq.test(bmi_re[,c(1,3)])$p.value))

## sum up ##
table.asp.use.asp.never= qpcR:::rbind.na(cc,age,bmi,edu,hhincome,fa_history,smk,diabetes,bph,statin,stroke,psa.scr,dre.scr,stage,gl,aggressive,bmi_re)
colnames(table.asp.use.asp.never)=c("asp.use, N","asp.use, %","asp.never, N","asp.never, %","P value")
write.csv(table.asp.use.asp.never,"table1.reg.asp.useVSasp.never.csv")



