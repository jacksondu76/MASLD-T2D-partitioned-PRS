
# MASLD genetic subtyping and GCKR: selected analysis code.
# This file contains part of the analysis workflow, rather than all study scripts.
# Update local data paths and provide your own LDlink token where needed.
###NAFLD AND T2D####
setwd('F:\\nafld-t2d')
library(data.table)
library(tidyverse)
library(cowplot)
library(survival)
library(ggsci)


# 1. Import SNP genotypes and align allele counts to the PDFF-increasing direction.
snps <- fread('nafld_snps_query_data.csv')
###convert all SNP allele to PDFF-increasing mode###
snps <- snps %>% mutate(
  rs1260326 = case_match(rs1260326,2~0,0~2,1~1),
  rs429358  = case_match(rs429358,2~0,0~2,1~1),
  rs28601761= case_match(rs28601761,2~0,0~2,1~1),
  rs2792751 = case_match(rs2792751,2~0,0~2,1~1),
  rs13389219= case_match(rs13389219,2~0,0~2,1~1),
  rs7029757 = case_match(rs7029757,2~0,0~2,1~1),
  rs11944752= case_match(rs11944752,2~0,0~2,1~1),
  rs2862954 = case_match(rs2862954,2~0,0~2,1~1)
)
pdff <- fread('F:\\E\\UkbData\\PDFF_40314.csv')
pdff <- pdff %>% select(eid.63726,pdff) %>% 
  left_join(snps %>% select(-c('ID_1','missing','sex')),
            by = c('eid.63726'='ID_2'))


# 2. Prepare UKBB phenotypes, covariates and the discovery dataset.
# Rename UKBB fields, handle missing values and merge genotype/PDFF data.
ukb <- fread('F:\\ukb_alco_mafld.csv')
#500938 184
ukb2 <- ukb %>% mutate_if(
  is.character,
  function(x) ifelse(x == '',NA,x)) %>% 
  mutate(
    # DiaInterval = time_length(ymd(f.41280.0.0)-ymd(f.53.0.0),'day'),
    DeathInterval = time_length(ymd(f.40000.0.0)-ymd(f.53.0.0),'day'),
    CancerInterval = time_length(ymd(f.40005.0.0)-ymd(f.53.0.0),'day'),
    survive = ifelse(is.na(f.40000.0.0),'survive','die')
  )

ukb3 <- ukb2 %>% rename(
  eid = f.eid,sex = f.31.0.0,YearOfBirth = f.34.0.0,
  attendingCentre = f.54.0.0,
  MonthOfBirth = f.52.0.0,
  WC = f.30000.0.0,RC = f.30010.0.0,HB = f.30020.0.0,
  HaePer = f.30030.0.0,MCV = f.30040.0.0,MCH = f.30050.0.0,
  MCHC = f.30060.0.0,RBCDW = f.30070.0.0,PLT = f.30080.0.0,
  PLTC = f.30090.0.0,MPV = f.30100.0.0,PLTDW = f.30100.0.0,
  LYM = f.30120.0.0,MON = f.30130.0.0,NEU = f.30140.0.0,
  EOS = f.30150.0.0,BAS = f.30160.0.0,NRBC = f.30170.0.0,
  LYMPer = f.30180.0.0,MonPer = f.30190.0.0,NeuPer = f.30200.0.0,
  EsoPer = f.30210.0.0,BasPer = f.30220.0.0,NRBCPer = f.30230.0.0,
  RetPer = f.30240.0.0,RET = f.30250.0.0,MRV = f.30260.0.0,
  MSCV = f.30270.0.0,IRF = f.30280.0.0,HLSRPer = f.30290.0.0,
  HLSR = f.30300.0.0,MALBUR = f.30500.0.0,CREUR = f.30510.0.0,
  KaUR = f.30520.0.0,NaUR = f.30530.0.0,ALB = f.30600.0.0,
  ALP = f.30610.0.0,ALT = f.30620.0.0,ApoA = f.30630.0.0,
  ApoB = f.30640.0.0,AST = f.30650.0.0,DBi = f.30660.0.0,
  UREA = f.30670.0.0,Ca = f.30680.0.0,CHL = f.30690.0.0,
  CRE = f.30700.0.0,CRP = f.30710.0.0,CysC = f.30720.0.0,
  GGT = f.30730.0.0,Glucose = f.30740.0.0,HbA1c = f.30750.0.0,
  HDL = f.30760.0.0,IGF1 = f.30770.0.0,LDL = f.30780.0.0,
  LipA = f.30790.0.0,OES = f.30800.0.0,PHOS = f.30810.0.0,
  RHE = f.30820.0.0,SHBG = f.30830.0.0,TBi = f.30840.0.0,
  TES = f.30850.0.0,TPro = f.30860.0.0,TRG = f.30870.0.0,
  Urate = f.30880.0.0,VD = f.30890.0.0,
  DateDeath = f.40000.0.0,DeathCause = f.40001.0.0,
  DateCaDia = f.40005.0.0,TypeofCa = f.40006.0.0,
  AgeDeath = f.40007.0.0,AgeCaDia = f.40008.0.0,
  HistoCa = f.40011.0.0
)
ukb3 <- ukb3 %>% rename(handGripStrength = f.47.0.0,waistCir = f.48.0.0,hipCir = f.49.0.0,
                        standingHeight = f.50.0.0,seatedHeight = f.51.0.0,dateAttendingCentre = f.53.0.0,
                        cancerAgeOccurred = f.84.0.0,noncancerAgeOccurred = f.87.0.0,
                        BirthWeightType = f.120.0.0,no.selfreportCA = f.134.0.0,
                        no.selfreportNonCA = f.135.0.0,no.selfreportTreat = f.137.0.0,
                        dateLossFollowup = f.191.0.0, averageTotalIncome = f.738.0.0,
                        ageCompleteFulltimeEdu = f.845.0.0,
                        no.daysWalk10 = f.864.0.0,
                        no.daysModeAct = f.884.0.0,
                        no.daysVigorAct = f.904.0.0,
                        sleepDuration = f.1160.0.0,insomnia = f.1200.0.0,
                        currentSmoke = f.1239.0.0,pastSmoke = f.1249.0.0,
                        aolIntakeFreq = f.1558.0.0,redwine1 = f.1568.0.0,
                        whitewine1 = f.1578.0.0,beer1 = f.1588.0.0,
                        spirits1 = f.1598.0.0,fortifiedwine1 = f.1608.0.0,
                        otherwine1 = f.5364.0.0,
                        redwine2 = f.4407.0.0,
                        whitewine2 = f.4418.0.0,beer2 = f.4429.0.0,
                        spirits2 = f.4440.0.0,fortifiedwine2 = f.4451.0.0,
                        otherwine2 = f.4462.0.0,
                        wineChange = f.1628.0.0,
                        CompBodysizeAt10 = f.1687.0.0,CompHeightsizeAt10= f.1697.0.0,
                        overallHealth = f.2178.0.0,longstandingIll = f.2188.0.0,
                        weightChange = f.2306.0.0,hadMajorOperation = f.2415.0.0,
                        diabetesBydoctor = f.2443.0.0,caBydoctor = f.2453.0.0,
                        otherSeriousIllBydoc = f.2473.0.0,
                        hadPrescription = f.2492.0.0,ageMenarche = f.2714.0.0,
                        hadOralContraceptivePill = f.2784.0.0,ageContraceptivePill=f.2794.0.0,
                        ageLastUseCP = f.2804.0.0,ageSmoking = f.2867.0.0,
                        forcedVitalCapacity = f.3062.0.0,FEV1 = f.3063.0.0,
                        PeakExpiratoryFlow = f.3064.0.0,coffeeDrink = f.3089.0.0,
                        DBP = f.4079.0.0,SBP = f.4080.0.0,
                        Nevereateggsdairywheatsugar = f.6144.0.0,
                        allowance = f.6146.0.0,dentalProbelm = f.6149.0.0,
                        vascularDisBydoc = f.6150.0.0,otherDisBydoc = f.6152.0.0,
                        medCholBPDiaF = f.6153.0.0,medPain = f.6154.0.0,
                        vitaMineralSupp = f.6155.0.0,medCholBPDiaM = f.6177.0.0,
                        otherDietarySupp = f.6179.0.0,birthweight = f.20022.0.0,
                        illnessMom = f.20110.0.0,smokingStatus = f.20116.0.0,
                        alcoholStatus = f.20117.0.0, NeuroticismScore = f.20127.0.0,
                        ethnic = f.21000.0.0,BMI = f.21001.0.0,weight = f.21002.0.0,
                        ageAssessment = f.21003.0.0,ageRecruitment = f.21022.0.0,
                        IBS = f.21024.0.0,
                        bodyFatPerc = f.23099.0.0,wholeBodyFatMass = f.23100.0.0,
                        wholeBodyFatFreeMass = f.23101.0.0, BMI2 = f.23104.0.0,
                        multipleDeprivationIndex = f.26410.0.0,
                        IncomeScore = f.26411.0.0,
                        healthScore = f.26413.0.0,
                        educationScore = f.26414.0.0,
                        housingScore = f.26415.0.0,
                        crimeScore = f.26416.0.0,
                        MPV = PLTDW, PLTDW = f.30110.0.0)
ukb3 <- ukb3 %>% mutate(
  age = year(dateAttendingCentre)-YearOfBirth,
  education = ifelse(!is.na(educationScore),educationScore,
                     ifelse(!is.na(f.26421.0.0),f.26421.0.0,f.26431.0.0))
)

ukb4 <- ukb3 %>% left_join(pdff,by = c('eid'='eid.63726'))
ukb4 <- ukb4 %>%
  mutate(
    education = ifelse(is.na(education),median(education,na.rm = T),education),
    averageTotalIncome = ifelse(is.na(averageTotalIncome)|averageTotalIncome<0,
                                median(averageTotalIncome,na.rm = T),averageTotalIncome),
    aolIntakeFreq = ifelse(!aolIntakeFreq %in% 1:6,3,aolIntakeFreq),
    BMI = ifelse(is.na(BMI),median(BMI,na.rm = T),BMI),
    smokingStatus = ifelse(!smokingStatus%in% 0:2,0,smokingStatus),
    no.daysModeAct = ifelse(no.daysModeAct<0|is.na(no.daysModeAct),3,no.daysModeAct),
    TRG = ifelse(is.na(TRG),median(TRG,na.rm = T),TRG),
    HbA1c = ifelse(is.na(HbA1c),median(HbA1c,na.rm = T),HbA1c),
    Glucose = ifelse(is.na(Glucose),median(Glucose,na.rm = T),Glucose),
    pdff2 = scale(pdff)
  ) %>% mutate(
    no.daysModeAct2 = ifelse(no.daysModeAct<2,'0-1',ifelse(no.daysModeAct<5,'2-4','5-7')),
    aolIntakeFreq2 = ifelse(aolIntakeFreq>=5,'Never',ifelse(aolIntakeFreq<=2,'Excess','Moderate')),
    aolIntakeFreq2 = factor(aolIntakeFreq2,levels = c('Never','Moderate','Excess')),
    TRG2= scale(TRG),
    HbA1c2 = scale(HbA1c),
    Glucose2 = scale(Glucose)
  )
ukb5 <- ukb4 %>% filter(ethnic %in% c(1001,1002,1003),
                        !is.na(rs738409))
#458407

# 3. Link NMR metabolites to the UKBB data using the project-ID mapping.
# Select the triglyceride-related metabolite fields for clustering.
metabolites <- fread('F:\\UKB-metabolites-updated.csv')
ukb_eid_match <- read.csv('F:\\ukb_match_63726_92718.csv',sep = '\t')
metabolites <- metabolites %>% left_join(ukb_eid_match,by = c('eid'='eid.92718'))
metabolites <- metabolites %>% filter(!is.na(eid.63726))
ukb5 <- ukb5 %>% left_join(metabolites,by = c('eid'='eid.63726'))

ukb_mets <- readxl::read_excel('F:\\ukb_mets.xlsx')
ukb_mets <- ukb_mets %>% mutate(
  metID= paste0('f.',eid,'.0.0')
) %>% rename(ID = eid)

tg <- ukb_mets$metID[grep('Triglyceride',ukb_mets$metabolites)]
rsids <- names(snps)[5:22]


# 4. Estimate covariate-adjusted SNP associations with TG-related metabolites.
# Save the association estimates in snps_tg.csv.
df2 <- data.frame()
for(i in 1:18){
  dat = ukb5 %>% select(age,sex,BMI,no.daysModeAct2,smokingStatus,aolIntakeFreq,
                        education,averageTotalIncome,rsids,tg[i])
  for(j in 1:18){
    snpi = scale(get(rsids[j],dat))
    met = scale(get(tg[i],dat))
    fit = lm(met~age+sex+BMI+no.daysModeAct2+smokingStatus+
               education+aolIntakeFreq+averageTotalIncome+snpi,data=dat)
    x <-  broom::tidy(fit)
    df <- data.frame(
      Met = c(tg[i]),
      snp = rsids[j],
      beta = c(x$estimate[12]),
      se= c(x$std.error[12]),
      p = c(x$p.value[12]))
    df2 <- rbind(df2,df)
  }
}

df2 <- df2 %>% left_join(ukb_mets %>% select(1,3),by = c('Met'='metID'))
write.csv(df2,'snps_tg.csv',row.names = F)


# 5. Construct SNP metabolic profiles and perform k-means clustering.
# Display the elbow plot, PCA, t-SNE and cluster-centre heatmap.
df2 <- read.csv('snps_tg.csv')
snp_mets <- df2 %>% mutate(
  fdr = p.adjust(p,'fdr'),
  beta2 = ifelse(fdr<0.05,beta,0)
)



library(factoextra)  
library(cluster)    
library(Rtsne)      
library(patchwork)  
library(pheatmap) 

data <- snp_mets

df <- data %>%
  select(snp, metabolites, beta) %>%
  spread(key = metabolites, value = beta, fill = 0)

rownames(df) <- df$snp
df <- df[, -1]
df_scaled <- scale(df)

elbow_plot <- fviz_nbclust(df_scaled, kmeans, method = "wss",)

set.seed(42)
k <- 4
kmeans_result <- kmeans(df_scaled, centers = k, nstart = 25)
pca_plot <- fviz_cluster(kmeans_result, data = df_scaled, 
                         geom = "point", ellipse.type = "convex", 
                         palette = c('#ED0000','#0099b4',"#42b540",'#00468b'),
                         ggtheme = theme_minimal())

set.seed(518)
tsne <- Rtsne(df_scaled, dims = 2, perplexity = 3)
tsne_df <- data.frame(tsne$Y, cluster = as.factor(kmeans_result$cluster))
tsne_plot <- ggplot(tsne_df, aes(x = X1, y = X2, color = cluster)) +
  scale_color_manual(values = c('#ED0000','#0099b4',"#42b540",'#00468b'))+
  geom_point(size = 3) +
  theme_minimal()

# sil <- silhouette(kmeans_result$cluster, dist(df_scaled))
# sil_plot <- fviz_silhouette(sil) +
#   ggtitle("D. Silhouette Plot of Clusters")

centers <- as.data.frame(kmeans_result$centers)
color_palette <- colorRampPalette(c("#0099b4", "white", "#ED0000"))(50)  

pheatmap(
  centers,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "row",  
  color = color_palette,  
  border_color = NA,  
  #main = "Cluster Centers Heatmap (Optimized Colors)",
  fontsize_row = 10,
  fontsize_col = 10,
  cellwidth = 25,
  cellheight = 15,
  treeheight_row = 50,
  treeheight_col = 50,
  angle_col = 45
)

(elbow_plot | pca_plot|tsne_plot )




# 6. Estimate individual-SNP associations with PDFF, TG, HbA1c and glucose.
# Save estimates and plot hepatic versus circulating TG effects.
df2 <- data.frame()
for(i in 1:18){
  dat = ukb5 %>% select(pdff2,rsids[i],age,sex,BMI,smokingStatus,aolIntakeFreq,
                        averageTotalIncome,no.daysModeAct2,education)
  fit1 <- lm(pdff2 ~ .,data = dat)
  x <-  broom::tidy(fit1)
  df <- data.frame(SNP = rsids[i],
                   Biomarker = 'PDFF',
                   beta = x$estimate[2],
                   se=x$std.error[2],
                   p = x$p.value[2])
  df2 <- bind_rows(df2,df)
}

for(i in 1:18){
  dat = ukb5 %>% select(TRG2,rsids[i],age,sex,BMI,smokingStatus,aolIntakeFreq,
                         averageTotalIncome,no.daysModeAct2,education)
  fit1 <- lm(TRG2 ~ .,data = dat)
  x <-  broom::tidy(fit1)
  df <- data.frame(SNP = rsids[i],
                   Biomarker = 'TRG',
                   beta = x$estimate[2],
                   se=x$std.error[2],
                   p = x$p.value[2])
  df2 <- bind_rows(df2,df)
}
for(i in 1:18){
  dat = ukb5 %>% select(HbA1c2,rsids[i],age,sex,BMI,smokingStatus,aolIntakeFreq,
                        averageTotalIncome,no.daysModeAct2,education)
  fit1 <- lm(HbA1c2 ~ .,data = dat)
  x <-  broom::tidy(fit1)
  df <- data.frame(SNP = rsids[i],
                   Biomarker = 'HbA1c',
                   beta = x$estimate[2],
                   se=x$std.error[2],
                   p = x$p.value[2])
  df2 <- bind_rows(df2,df)
}
for(i in 1:18){
  dat = ukb5 %>% select(Glucose2,rsids[i],age,sex,BMI,smokingStatus,aolIntakeFreq,
                        averageTotalIncome,no.daysModeAct2,education)
  fit1 <- lm(Glucose2 ~ .,data = dat)
  x <-  broom::tidy(fit1)
  df <- data.frame(SNP = rsids[i],
                   Biomarker = 'FBG',
                   beta = x$estimate[2],
                   se=x$std.error[2],
                   p = x$p.value[2])
  df2 <- bind_rows(df2,df)
}

write.csv(df2,'snp_bio.csv',row.names = F)

df2 <- read.csv('snp_bio.csv')
# Load the data

df3 <- df2 %>% filter(Biomarker %in% c('PDFF','TRG')) %>%
  mutate(
    beta = ifelse(p>0.05,0,beta)
  ) %>% pivot_wider(names_from = 'Biomarker',values_from = beta:p)
df3 <- df3 %>% mutate(
  group = ifelse(SNP %in% c('rs1260326','rs28601761'),'b',
                 ifelse(SNP %in% c('rs58542926','rs187429064'),'a',
                        ifelse(SNP %in% c('rs13389219','rs11944752',
                                          'rs140201358'),'d','c')))
)
b <- ggplot(df3,aes(beta_PDFF,beta_TRG,color=group))+
  geom_hline(yintercept = 0,linetype=2)+
  geom_point(size=2,show.legend=F)+
  scale_color_lancet()+
  scale_x_continuous(limits = c(0,0.35),expand = expansion(add=0))+
  ggrepel::geom_text_repel(aes(label=SNP),max.overlaps = 20,
                           show.legend = F)+
  theme_bw()+
  theme(
    axis.text = element_text(size = 11,color='black'),
    axis.title = element_text(size=12))+
  labs(x = 'Effect size on hepatic TG',
       y = 'Effect size on circulating TG')


ukb5 <- ukb5 %>%
  mutate(

# 7. Prepare additional biomarkers and calculate the four weighted PRSs.
    CHL = ifelse(is.na(CHL),median(CHL,na.rm = T),CHL),
    LDL = ifelse(is.na(LDL),median(LDL,na.rm = T),LDL),
    HDL = ifelse(is.na(HDL),median(HDL,na.rm = T),HDL),
    ALT = ifelse(is.na(ALT),median(ALT,na.rm = T),ALT),
    AST = ifelse(is.na(AST),median(AST,na.rm = T),AST),
    SBP = ifelse(is.na(SBP),median(SBP,na.rm = T),SBP),
    DBP = ifelse(is.na(DBP),median(DBP,na.rm = T),DBP),
    GGT = ifelse(is.na(GGT),median(GGT,na.rm = T),GGT),
    TRG = ifelse(is.na(TRG),median(TRG,na.rm = T),TRG),
    PLT = ifelse(is.na(PLT),median(PLT,na.rm = T),PLT),
    ALB = ifelse(is.na(ALB),median(ALB,na.rm = T),ALB),
    waistCir = ifelse(is.na(waistCir),median(waistCir,na.rm = T),waistCir)
  )

ukb5 <- ukb5 %>% mutate(
  prs1 = rs1260326*0.067+rs28601761*0.072,
  prs4 = rs58542926*0.34+rs187429064*0.39,
  prs2 = rs13389219*0.054+rs140201358*0.22+rs11944752*0.047,
  prs3 = rs429358*0.110+rs738409*0.28+rs2792751*0.066+
    rs641738*0.052+rs2642438*0.053+rs7029757*0.078+rs1801689*0.13+
    rs1229984*0.14+
    +rs6955582*0.038+rs1800562*0.071+rs2862954*0.031
)



# 8. Define the incident cohort and follow-up for T2D and liver outcomes.
# Exclude baseline disease and merge the incident diagnosis records.
###association between PRS and T2D###
diagData <- fread('F:\\E\\UkbData\\diagData_icd10_20221031.csv')
diagData <- diagData %>% dplyr::rename(eid = eid.63726)
baseline_ca <- diagData %>%
  dplyr:: filter(followUp_days<=0,diseases  %in% c('Cancer')) %>%
  dplyr::select(eid)

baseline_cvd <- diagData %>% dplyr::filter(grepl('I20|I21|I22|I23|I24|I25|I60|I61|I62|I63|I64|I69',Diagnosis)) %>%
  dplyr::filter(followUp_days<=0) %>%
  dplyr::select(eid)

baseline_t2d <- diagData %>%dplyr::filter(grepl('E11',Diagnosis)) %>%
  dplyr::filter(followUp_days<=0) %>%
  dplyr::select(eid)

baseline_cld <- diagData %>%
  filter(followUp_days<=0,diseases  %in% c('MLD','MSLD')) %>%
  select(eid)

baseline_diseases <- union(baseline_ca$eid,baseline_cvd$eid) %>%
  union(baseline_t2d$eid) %>% union(baseline_cld$eid)

ukb6 <- ukb5 %>% filter(!eid %in% baseline_diseases)
#n = 415167

T2D <- diagData %>% filter(grepl('E11',Diagnosis)) %>% 
  select(eid,DiagDate,followUp_days) %>%
  rename(followUp_days_T2D = followUp_days) %>%
  arrange(eid,DiagDate) %>%
  distinct(eid,.keep_all = TRUE) %>%
  filter(followUp_days_T2D>0) %>%
  mutate(T2D=1)

HCC <- diagData %>% filter(HCC==1) %>%
  select(eid,DiagDate,followUp_days,HCC) %>%
  rename(followUp_days_HCC = followUp_days) %>%
  arrange(eid,DiagDate) %>%
  distinct(eid,.keep_all = TRUE) %>%
  filter(followUp_days_HCC>0)

Cirrhosis <- diagData %>% filter(Cirrhosis==1) %>%
  select(eid,DiagDate,followUp_days,Cirrhosis) %>%
  rename(followUp_days_CIRR = followUp_days) %>%
  arrange(eid,DiagDate) %>%
  distinct(eid,.keep_all = TRUE) %>%
  filter(followUp_days_CIRR>0)

MASLD <- diagData %>% filter(grepl('K760|K758',Diagnosis)) %>%
  select(eid,DiagDate,followUp_days) %>%
  rename(followUp_days_MASLD = followUp_days) %>%
  arrange(eid,DiagDate) %>%
  distinct(eid,.keep_all = TRUE) %>%
  filter(followUp_days_MASLD>0) %>%
  mutate(MASLD=1)

ukb6 <- ukb6 %>% 
  left_join(T2D,by='eid') %>%
  #left_join(CVD,by='eid') %>%
  #left_join(CKD,by='eid') %>%
  left_join(MASLD,by='eid') %>%
  left_join(HCC,by='eid') %>%
  left_join(Cirrhosis,by='eid') %>%
  mutate(followUp_days_T2D = ifelse(is.na(T2D),time_length(ymd('2022-10-31')-ymd(dateAttendingCentre),'day'),followUp_days_T2D),
         T2D = ifelse(is.na(T2D),0,1)) %>%
  mutate(followUp_days_MASLD = ifelse(is.na(MASLD),time_length(ymd('2022-10-31')-ymd(dateAttendingCentre),'day'),followUp_days_MASLD),
         MASLD = ifelse(is.na(MASLD),0,1)) %>%
  # mutate(followUp_days_CKD = ifelse(is.na(CKD),time_length(ymd('2022-10-31')-ymd(dateAttendingCentre),'day'),followUp_days_CKD),
  #        CKD = ifelse(is.na(CKD),0,1)) %>%
  mutate(followUp_days_HCC = ifelse(is.na(HCC),time_length(ymd('2022-10-31')-ymd(dateAttendingCentre),'day'),followUp_days_HCC),
         HCC = ifelse(is.na(HCC),0,1)) %>%
  # mutate(followUp_days_stroke = ifelse(is.na(stroke),time_length(ymd('2022-10-31')-ymd(dateAttendingCentre),'day'),followUp_days_stroke),
  #        stroke = ifelse(is.na(stroke),0,1))
  mutate(followUp_days_CIRR = ifelse(is.na(Cirrhosis),time_length(ymd('2022-10-31')-ymd(dateAttendingCentre),'day'),followUp_days_CIRR),
       Cirrhosis = ifelse(is.na(Cirrhosis),0,1))


ukb6 <- ukb6 %>% mutate(
  loweringINS = ifelse(((sex==1 & medCholBPDiaM==3) |(sex==0 & medCholBPDiaF==3)),1,0)
  
)
fwrite(ukb6,'analyzingData.csv')


# 9. Estimate covariate-adjusted associations between PRSs and biomarkers.
ukb6 <- fread('analyzingData.csv')
ukb6 <- ukb6 %>% mutate(
  IR = TRG/HDL
)
###associations between PRS and biomarkers@###
biomarker_names <- c('Alanine aminotransferase',
                     'Aspartate aminotransferase','Gamma glutamyltransferase',
                     'Direct bilirubin',
                     'Total bilirubin','Fasting blood glucose',
                     'Glycosylated hemoglobin')
biomarkers <- c('ALT',
                'AST','GGT','DBi','TBi',
                'Glucose','HbA1c')
df2 <- data.frame()
for(i in 1:7){
  y <- get(biomarkers[i],pos = ukb6)
  fit <- lm(scale(y)~scale(prs1)+sex+age+smokingStatus+averageTotalIncome+
              aolIntakeFreq2+no.daysModeAct2+education+BMI,data = ukb6)
  x <-  broom::tidy(fit)
  df <- data.frame(SNP = 'PRS1',
                   Biomarker = biomarkers[i],
                   BioName   = biomarker_names[i],
                   beta = x$estimate[2],
                   se=x$std.error[2],
                   p = x$p.value[2])
  df2 <- bind_rows(df2,df)
}
for(i in 1:7){
  y <- get(biomarkers[i],pos = ukb6)
  fit <- lm(scale(y)~scale(prs2)+sex+age+smokingStatus+averageTotalIncome+
              aolIntakeFreq2+no.daysModeAct2+education+BMI,data = ukb6)
  x <-  broom::tidy(fit)
  df <- data.frame(SNP = 'PRS2',
                   Biomarker = biomarkers[i],
                   BioName   = biomarker_names[i],
                   beta = x$estimate[2],
                   se=x$std.error[2],
                   p = x$p.value[2])
  df2 <- bind_rows(df2,df)
}
for(i in 1:7){
  y <- get(biomarkers[i],pos = ukb6)
  fit <- lm(scale(y)~scale(prs3)+sex+age+smokingStatus+averageTotalIncome+
              aolIntakeFreq2+no.daysModeAct2+education+BMI,data = ukb6)
  x <-  broom::tidy(fit)
  df <- data.frame(SNP = 'PRS3',
                   Biomarker = biomarkers[i],
                   BioName   = biomarker_names[i],
                   beta = x$estimate[2],
                   se=x$std.error[2],
                   p = x$p.value[2])
  df2 <- bind_rows(df2,df)
}
for(i in 1:7){
  y <- get(biomarkers[i],pos = ukb6)
  fit <- lm(scale(y)~scale(prs4)+sex+age+smokingStatus+averageTotalIncome+
              aolIntakeFreq2+no.daysModeAct2+education+BMI,data = ukb6)
  x <-  broom::tidy(fit)
  df <- data.frame(SNP = 'PRS4',
                   Biomarker = biomarkers[i],
                   BioName   = biomarker_names[i],
                   beta = x$estimate[2],
                   se=x$std.error[2],
                   p = x$p.value[2])
  df2 <- bind_rows(df2,df)
}


write.csv(df2,'PRS_biomarkers.csv',row.names = F)

df2 <- read.csv('PRS_biomarkers.csv')
df2 <- df2 %>% mutate(
  BioName = factor(BioName,
                   levels = biomarker_names),
  SNP = factor(SNP,levels = c('PRS4','PRS3','PRS2','PRS1'))
)
a <- ggplot(df2,aes(beta,BioName,color=SNP))+
  annotate('rect',xmin = rep(-0.035,4),xmax = rep(0.1,4),
           ymin = seq(0.5,6.5,by=2),
           ymax = seq(1.5,7.5,by=2),fill = '#0081cf',
           alpha=.1,color=NA)+
  geom_vline(xintercept = 0,linetype=2)+
  geom_errorbarh(aes(xmin=lwr,xmax=upr),height=0,
                 position = position_dodge(width = 0.75),
                 linewidth=.5)+
  geom_point(size=1,position = position_dodge(width=.75))+
  theme_test()+
  scale_color_manual(name = 'Partitioned PRS',
                     values = c('#ED0000','#0099b4',"#42b540",'#00468b'),
                     breaks = c('PRS1','PRS2','PRS3','PRS4'))+
  scale_x_continuous(limits = c(-0.035,0.1),
                     expand = expansion(add=0),
                     breaks = seq(-0.03,0.1,by=0.03))+
  scale_y_discrete(limits = rev,
                   expand = expansion(add=c(0,0.5)))+
  labs(x = 'Effect size (95% CI)',
       y = NULL)+
  theme(axis.text = element_text(size = 11,color='black'),
        axis.title = element_text(size=12))



# 10. Fit Cox models for incident disease, including all four PRSs together.
# Summarise hazard ratios and confidence intervals and draw forest plots.
fit_cox <- coxph(Surv(followUp_days_T2D,T2D)~sex+age+averageTotalIncome+aolIntakeFreq2+
                   BMI+education+smokingStatus+no.daysModeAct2+loweringINS+
                   prs1+prs2+prs3+prs4,
                 data = ukb6)
x <- broom::tidy(fit_cox)

df <- data.frame(
  prs     = c('PRS1','PRS2','PRS3','PRS4'),
  outcome = 'T2D',
  hr      = exp(x$estimate[12:15]),
  lwr     = exp(x$estimate[12:15]-1.96*x$std.error[12:15]),
  upr     = exp(x$estimate[12:15]+1.96*x$std.error[12:15])
)


fit_cox <- coxph(Surv(followUp_days_HCC,HCC)~sex+age+averageTotalIncome+aolIntakeFreq2+
                   BMI+education+smokingStatus+no.daysModeAct2+loweringINS+
                   prs1+prs2+prs3+prs4,
                 data = ukb6)
x <- broom::tidy(fit_cox)

df <- df %>% bind_rows(
  data.frame(
  prs     = c('PRS1','PRS2','PRS3','PRS4'),
  outcome = 'HCC',
  hr      = exp(x$estimate[12:15]),
  lwr     = exp(x$estimate[12:15]-1.96*x$std.error[12:15]),
  upr     = exp(x$estimate[12:15]+1.96*x$std.error[12:15])
))


fit_cox <- coxph(Surv(followUp_days_CIRR,Cirrhosis)~sex+age+averageTotalIncome+aolIntakeFreq2+
                   BMI+education+smokingStatus+no.daysModeAct2+loweringINS+
                   prs1+prs2+prs3+prs4,
                 data = ukb6)
x <- broom::tidy(fit_cox)

df <- df %>% bind_rows(
  data.frame(
    prs     = c('PRS1','PRS2','PRS3','PRS4'),
    outcome = 'Cirrhosis',
    hr      = exp(x$estimate[12:15]),
    lwr     = exp(x$estimate[12:15]-1.96*x$std.error[12:15]),
    upr     = exp(x$estimate[12:15]+1.96*x$std.error[12:15])
  ))

fit_cox <- coxph(Surv(followUp_days_MASLD,MASLD)~sex+age+averageTotalIncome+aolIntakeFreq2+
                   BMI+education+smokingStatus+no.daysModeAct2+loweringINS+
                   prs1+prs2+prs3+prs4,
                 data = ukb6)
x <- broom::tidy(fit_cox)

df <- df %>% bind_rows(
  data.frame(
    prs     = c('PRS1','PRS2','PRS3','PRS4'),
    outcome = 'MASLD',
    hr      = exp(x$estimate[12:15]),
    lwr     = exp(x$estimate[12:15]-1.96*x$std.error[12:15]),
    upr     = exp(x$estimate[12:15]+1.96*x$std.error[12:15])
  ))


write.csv(df,'prs_hr.csv',row.names = F)
df <- read.csv('prs_hr.csv')
df <- df %>% mutate(
  prs = factor(prs,levels = c('PRS4','PRS3','PRS2','PRS1')),
  outcome=factor(outcome,levels = c('MASLD','HCC','Cirrhosis','T2D'))
)

b <- ggplot(df,aes(hr,outcome,color=prs))+
  annotate('rect',xmin = rep(0.2,2),xmax = rep(30,2),
           ymin = seq(0.5,3.5,by=2),
           ymax = seq(1.5,4.5,by=2),fill = '#0081cf',
           alpha=.1,color=NA)+
  geom_vline(xintercept = 1,linetype=2)+
  geom_errorbarh(aes(xmin=lwr,xmax=upr),height=0,
                 position = position_dodge(width = 0.75),
                 linewidth=.5)+
  geom_point(size=1,position = position_dodge(width=.75))+
  theme_test()+
  scale_color_manual(name = 'Partitioned PRS',
                     values = c('#ED0000','#0099b4',"#42b540",'#00468b'),
                     breaks = c('PRS1','PRS2','PRS3','PRS4'))+
  scale_x_continuous(limits = c(0.2,30),
                     expand = expansion(add=0),
                     breaks = c(1,2,5,10,20),
                     trans = "log10")+
  scale_y_discrete(limits = rev,
                   expand = expansion(add=c(0,0.5)))+
  labs(x = 'Hazard ratio (95% CI)',
       y = NULL)+
  theme(axis.text = element_text(size = 11,color='black'),
        axis.title = element_text(size=12))
legend <- get_legend(a)
plot_grid(plot_grid(a+theme(legend.position = 'none'),
          b+theme(legend.position = 'none'),
          ncol = 2,rel_widths = c(1,0.73),
          labels = 'AUTO'),
          legend,ncol = 2,rel_widths = c(1,.15))



# 11. Examine combined and alternative genetic-score groupings.
###integrating prs2-4 together
ukb6 <- ukb6 %>% mutate(
  prs234 = rs58542926*0.34+rs187429064*0.39+rs13389219*0.054+
    rs140201358*0.22+rs11944752*0.047+rs429358*0.110+rs738409*0.28+
    rs2792751*0.066+rs641738*0.052+rs2642438*0.053+rs7029757*0.078+
    rs1801689*0.13+rs1229984*0.14+
    +rs6955582*0.038+rs1800562*0.071+rs2862954*0.031
)

fit_cox <- coxph(Surv(followUp_days_T2D,T2D)~sex+age+averageTotalIncome+aolIntakeFreq2+
                   BMI+education+smokingStatus+no.daysModeAct2+loweringINS+
                   prs1+prs234,
                 data = ukb6)
summary(fit_cox)

###cutting prs by direction of pdff and tg
ukb6 <- ukb6 %>% mutate(
  prs13 = rs1260326*0.067+rs28601761*0.072+rs13389219*0.054+
    rs140201358*0.22+rs11944752*0.047,
  prs24 = rs58542926*0.34+rs187429064*0.39+rs429358*0.110+rs738409*0.28+
    rs2792751*0.066+rs641738*0.052+rs2642438*0.053+
    rs1801689*0.13+rs1229984*0.14+
    +rs6955582*0.038+rs2862954*0.031
)

fit_cox <- coxph(Surv(followUp_days_T2D,T2D)~sex+age+averageTotalIncome+aolIntakeFreq2+
                   BMI+education+smokingStatus+no.daysModeAct2+loweringINS+
                   prs13+prs24,
                 data = ukb6)
summary(fit_cox)

# prs13                    1.3431       1.1716    1.5397
# prs24                    1.2433       1.1820    1.3099


# 12. Plot T2D cumulative risk across PRS categories.
# Predictions below use the specified age and sex values.
###cumulative incidence of T2D
ukb6 <- ukb6 %>% mutate(
  prs1G = ifelse(prs1>0.211,'> 90%',ifelse(prs1<0.067,'0-10%','10-90%'))
)

model1 <- coxph(Surv(followUp_days_T2D,T2D)~prs1G+sex+age,
                ukb6)
new_df <- with(ukb6,
               data.frame(prs1G = c('0-10%','10-90%','> 90%'),
                          sex = 1,
                          age = 60
               ))

fit <- survfit(model1, newdata = new_df)
dat = data.frame(surv = fit$surv,lower= fit$lower,
                 upper = fit$upper,time= fit$time)
data = melt(dat, id="time")
data_wide <- data %>%
  separate(col = variable, into = c("type", "pPRS1"), sep = "\\.") %>%
  spread(key = type, value = value)

a <- ggplot(data = data_wide) +
  geom_line(aes(x = time, y = 1-surv, group = pPRS1,
                colour = pPRS1),
            size = 1.2,
            linetype = 1) +
  geom_ribbon(aes(x = time, ymin = 1-lower, ymax = 1-upper,
                  group = pPRS1, fill = pPRS1),
              alpha = 0.5) +
  ggsci::scale_colour_npg(name ='FBG-lowering PRS',
                          labels = c('0-10%',
                                     '10-90%',
                                     '>90%')) +
  #scale_fill_manual(values = c('#112d4e','#3f72af')) +
  scale_x_continuous(limits = c(0, 3650),
                     breaks = seq(0, 3650, 365),
                     labels = seq(0,10,by=1)) +
  theme_bw()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = 'inside',
        legend.position.inside = c(.15,.8)) +
  ylab("Cumulative probability")+
  xlab('Years of follow-up')+
  scale_y_continuous(labels = scales::percent,
                     limits = c(0,0.075))


ukb6 <- ukb6 %>% mutate(
  prs1G = ifelse(prs1>0.211,'> 90%',ifelse(prs1<0.067,'0-10%','10-90%'))
)

model1 <- coxph(Surv(followUp_days_T2D,T2D)~prs1G+sex+age,
                ukb6)
new_df <- with(ukb6,
               data.frame(prs1G = c('0-10%','10-90%','> 90%'),
                          sex = 1,
                          age = 60
               ))

fit <- survfit(model1, newdata = new_df)
dat = data.frame(surv = fit$surv,lower= fit$lower,
                 upper = fit$upper,time= fit$time)
data = melt(dat, id="time")
data_wide <- data %>%
  separate(col = variable, into = c("type", "pPRS1"), sep = "\\.") %>%
  spread(key = type, value = value)

a <- ggplot(data = data_wide) +
  geom_line(aes(x = time, y = 1-surv, group = pPRS1,
                colour = pPRS1),
            size = 1.2,
            linetype = 1) +
  # geom_ribbon(aes(x = time, ymin = 1-lower, ymax = 1-upper,
  #                 group = pPRS1, fill = pPRS1),
  #             alpha = 0.5) +
  ggsci::scale_colour_npg(name ='FBG-lowering PRS',
                          labels = c('0-10%',
                                     '10-90%',
                                     '>90%')) +
  #scale_fill_manual(values = c('#112d4e','#3f72af')) +
  scale_x_continuous(limits = c(0, 5475),
                     breaks = seq(0, 5475, by=365*3),
                     labels = seq(0,15,by=3)) +
  theme_bw()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = 'inside',
        legend.position.inside = c(.3,.8)) +
  ylab("Cumulative probability")+
  xlab('Years of follow-up')+
  scale_y_continuous(labels = scales::percent,
                     limits = c(0,0.125))
data_wide %>% filter(time == 365*15-1)
1-c(0.9031676,0.9077012,0.9116620)
#9.68% vs 8.83%

inc_top10 <- 0.0883   # Top 10% 发病率

# Additional calculation using the entered 15-year risks and group sizes.
inc_bottom10 <- 0.0968  # Bottom 10% 发病率
n_top10 <- 21873     # Top 10% 样本量
n_bottom10 <- 26912   # Bottom 10% 样本量
#P值: 0.03177275
# 计算率差 (Risk Difference)
RD <- inc_top10 - inc_bottom10
cat("率差 (RD):", round(RD * 100, 3), "%\n")

# 计算标准误 (SE)
SE_RD <- sqrt((inc_top10 * (1 - inc_top10) / n_top10) + (inc_bottom10 * (1 - inc_bottom10) / n_bottom10))
cat("标准误 (SE):", SE_RD, "\n")

# 计算Z值
Z <- RD / SE_RD
cat("Z值:", Z, "\n")

# 计算P值（双侧）
P_value <- 2 * (1 - pnorm(abs(Z)))
cat("P值:", P_value, "\n")


ukb6 <- ukb6 %>% mutate(
  prs234 = rs58542926*0.34+rs187429064*0.39+rs13389219*0.054+
    rs140201358*0.22+rs11944752*0.047+rs429358*0.110+rs738409*0.28+
    rs2792751*0.066+rs641738*0.052+rs2642438*0.053+rs7029757*0.078+
    rs1801689*0.13+rs1229984*0.14+
    +rs6955582*0.038+rs1800562*0.071+rs2862954*0.031
) %>% mutate(
  prs234G = ifelse(prs234>1.51,'> 90%',ifelse(prs234<0.884,'0-10%','10-90%'))
)

model1 <- coxph(Surv(followUp_days_T2D,T2D)~prs234G+sex+age,
                ukb6)
new_df <- with(ukb6,
               data.frame(prs234G = c('0-10%','10-90%','> 90%'),
                          sex = 1,
                          age = 60
               ))

fit <- survfit(model1, newdata = new_df)
dat = data.frame(surv = fit$surv,lower= fit$lower,
                 upper = fit$upper,time= fit$time)
data = melt(dat, id="time")
data_wide2 <- data %>%
  separate(col = variable, into = c("type", "pPRS1"), sep = "\\.") %>%
  spread(key = type, value = value)

b <- ggplot(data = data_wide2) +
  geom_line(aes(x = time, y = 1-surv, group = pPRS1,
                colour = pPRS1),
            size = 1.2,
            linetype = 1) +
  # geom_ribbon(aes(x = time, ymin = 1-lower, ymax = 1-upper,
  #                 group = pPRS1, fill = pPRS1),
  #             alpha = 0.5) +
  ggsci::scale_colour_npg(name ='FBG-increasing PRS',
                          labels = c('0-10%',
                                     '10-90%',
                                     '>90%')) +
  #scale_fill_manual(values = c('#112d4e','#3f72af')) +
  scale_x_continuous(limits = c(0, 5475),
                     breaks = seq(0, 5475, by=365*3),
                     labels = seq(0,15,by=3)) +
  theme_bw()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = 'inside',
        legend.position.inside = c(.25,.8)) +
  ylab("Cumulative probability")+
  xlab('Years of follow-up')+
  scale_y_continuous(labels = scales::percent,
                     limits = c(0,0.125))
data_wide2 %>% filter(time == 365*15-1)
1-c(0.9161051,0.8956013)
#8.38% vs 10.44%
plot_grid(a,b,ncol = 2,labels = 'auto')


# Plot TZL prevalence estimates supplied from the separate validation analyses.
df <- data.frame(
  PRS = factor(c('0-10%','10-90%','>90%'),levels = c('0-10%','10-90%','>90%')),
  prevalence=c(8.38,6.22,5.83),
  lwr = c(7.69,5.99,5.27),
  upr = c(9.01,6.56,6.40)
)
df2 <- data.frame(
  PRS = factor(c('0-10%','10-90%','>90%'),levels = c('0-10%','10-90%','>90%')),
  prevalence=c(4.38,7.22,9.11),
  lwr = c(2.89,5.99,8.27),
  upr = c(5.86,8.42,9.95)
)
a <- ggplot(df,aes(PRS,weight=prevalence,fill=PRS))+
  geom_bar(width=0.6)+
  geom_errorbar(aes(ymin = lwr,ymax=upr),width=0.1)+
  theme_classic()+
  scale_fill_npg()+
  scale_y_continuous(expand = expansion(add = c(0,0)),
                     limits = c(0,10),
                     breaks = c(0,2.5,5,7.5,10))+
  theme(axis.text = element_text(color='black'))+
  labs(x = 'FBG-lowering PRS category',
       y = 'T2D prevalence (95% CI)')
b <- ggplot(df2,aes(PRS,weight=prevalence,fill=PRS))+
  geom_bar(width=0.6)+
  geom_errorbar(aes(ymin = lwr,ymax=upr),width=0.1)+
  theme_classic()+
  scale_fill_npg()+
  scale_y_continuous(expand = expansion(add = c(0,0)),
                     limits = c(0,10),
                     breaks = c(0,2.5,5,7.5,10))+
  theme(axis.text = element_text(color='black'))+
  labs(x = 'FBG-increasing PRS category',
       y = 'T2D prevalence (95% CI)')
plot_grid(a,b,ncol = 2)



# 13. Perform two-sample MR using liver-fat SNPs and cluster-specific instruments.
# Prepare summary statistics, harmonise alleles and estimate MR associations.
####MR analysis
library(TwoSampleMR)
nafld <- readxl::read_excel('nafld-snps.xlsx',sheet = 2)
nafld <- nafld %>% rename(rsid=SNP) %>%
  mutate(
  eaf = eaf/100,
  SNP = paste0(chr,":",pos)
)
t2d <- fread('Mahajan.NG2018b.T2D.European.txt.gz')
t2d2 <- t2d %>% filter(SNP %in% nafld$SNP) %>% 
  left_join(nafld %>% select(rsid,SNP),by='SNP') %>%
  as.data.frame()

nafld_mr <- format_data(nafld,
                      type = 'exposure',
                      snp_col =  'rsid',
                      beta_col = 'beta',
                      se_col = 'se',
                      eaf_col = 'eaf',
                      effect_allele_col = 'ea',
                      other_allele_col = 'oa',
                      pval_col = 'p',
                      chr_col = 'chr',
                      pos_col = 'pos')
t2d_mr <- format_data(t2d2,
                        type = 'outcome',
                        snp_col =  'rsid',
                        beta_col = 'Beta',
                        se_col = 'SE',
                        eaf_col = 'EAF',
                        effect_allele_col = 'EA',
                        other_allele_col = 'NEA',
                        pval_col = 'Pvalue',
                        chr_col = 'Chr',
                        pos_col = 'Pos')

r2 = 2*nafld_mr$beta.exposure^2*nafld_mr$eaf.exposure*(1-nafld_mr$eaf.exposure)
mean(r2*(33592-2)/(1-r2))
#98.1

h_dat <- harmonise_data(nafld_mr,t2d_mr)
mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
                                    'mr_egger_regression',
                                    'mr_weighted_median',
                                    'mr_weighted_mode'))
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- data.frame(
  exposure = 'NAFLD',
  outcome  = 'T2D',
  method   = mr_res$method,
  beta     = mr_res$b,
  se       = mr_res$se,
  p        = mr_res$pval,
  p_Qval   = c(h_test$Q_pval[2],NA,NA,NA),
  p_pleio  = c(p_test$pval,NA,NA,NA)
)

prs1 <- nafld_mr %>% filter(SNP %in% c('rs1260326','rs28601761'))
prs4 <- nafld_mr %>% filter(SNP %in% c('rs58542926','rs187429064'))
prs3 <- nafld_mr %>% filter(SNP %in% c('rs429358','rs738409','rs2792751',
                                       'rs641738','rs2642438','rs7029757',
                                       'rs1801689','rs1229984','rs6955582',
                                       'rs1800562','rs2862954'))
prs2 <- nafld_mr %>% filter(SNP %in% c('rs13389219','rs140201358','rs11944752'))

h_dat <- harmonise_data(prs1,t2d_mr)
mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
                                  'mr_egger_regression',
                                  'mr_weighted_median',
                                  'mr_weighted_mode'))
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- df %>% bind_rows(
  data.frame(
  exposure = 'PRS1',
  outcome  = 'T2D',
  method   = mr_res$method,
  beta     = mr_res$b,
  se       = mr_res$se,
  p        = mr_res$pval,
  p_Qval   = c(h_test$Q_pval[1]),
  p_pleio  = c(p_test$pval)
))

h_dat <- harmonise_data(prs2,t2d_mr)
mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
                                  'mr_egger_regression',
                                  'mr_weighted_median',
                                  'mr_weighted_mode'))
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- df %>% bind_rows(
  data.frame(
    exposure = 'PRS2',
    outcome  = 'T2D',
    method   = mr_res$method,
    beta     = mr_res$b,
    se       = mr_res$se,
    p        = mr_res$pval,
    p_Qval   = c(h_test$Q_pval[1]),
    p_pleio  = c(p_test$pval)
  ))
h_dat <- harmonise_data(prs3,t2d_mr)
mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
                                  'mr_egger_regression',
                                  'mr_weighted_median',
                                  'mr_weighted_mode'))
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- df %>% bind_rows(
  data.frame(
    exposure = 'PRS3',
    outcome  = 'T2D',
    method   = mr_res$method,
    beta     = mr_res$b,
    se       = mr_res$se,
    p        = mr_res$pval,
    p_Qval   = c(h_test$Q_pval[1]),
    p_pleio  = c(p_test$pval)
  ))
h_dat <- harmonise_data(prs4,t2d_mr)
mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
                                  'mr_egger_regression',
                                  'mr_weighted_median',
                                  'mr_weighted_mode'))
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- df %>% bind_rows(
  data.frame(
    exposure = 'PRS4',
    outcome  = 'T2D',
    method   = mr_res$method,
    beta     = mr_res$b,
    se       = mr_res$se,
    p        = mr_res$pval,
    p_Qval   = c(h_test$Q_pval[1]),
    p_pleio  = c(p_test$pval)
  ))


# MR outcome: fasting blood glucose.
glucose <- fread('MAGIC1000G_FG_EUR.tsv.gz')
glucose2 <- glucose %>% filter(variant %in% nafld$rsid)

glu_mr <- format_data(as.data.frame(glucose2),
                      type = 'outcome',
                      snp_col =  'variant',
                      beta_col = 'beta',
                      se_col = 'standard_error',
                      eaf_col = 'effect_allele_frequency',
                      effect_allele_col = 'effect_allele',
                      other_allele_col = 'other_allele',
                      pval_col = 'p_value',
                      chr_col = 'chromosome',
                      pos_col = 'base_pair_location')

h_dat <- harmonise_data(nafld_mr,glu_mr)
mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
                                  'mr_egger_regression',
                                  'mr_weighted_median',
                                  'mr_weighted_mode'))
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- df %>% bind_rows(
  data.frame(
  exposure = 'NAFLD',
  outcome  = 'FBG',
  method   = mr_res$method,
  beta     = mr_res$b,
  se       = mr_res$se,
  p        = mr_res$pval,
  p_Qval   = c(h_test$Q_pval[2],NA,NA,NA),
  p_pleio  = c(p_test$pval,NA,NA,NA)
))
h_dat <- harmonise_data(prs1,glu_mr)
mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
                                  'mr_egger_regression',
                                  'mr_weighted_median',
                                  'mr_weighted_mode'))
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- df %>% bind_rows(
  data.frame(
    exposure = 'PRS1',
    outcome  = 'FBG',
    method   = mr_res$method,
    beta     = mr_res$b,
    se       = mr_res$se,
    p        = mr_res$pval,
    p_Qval   = c(h_test$Q_pval[1]),
    p_pleio  = c(p_test$pval)
  ))

h_dat <- harmonise_data(prs2,glu_mr)
mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
                                  'mr_egger_regression',
                                  'mr_weighted_median',
                                  'mr_weighted_mode'))
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- df %>% bind_rows(
  data.frame(
    exposure = 'PRS2',
    outcome  = 'FBG',
    method   = mr_res$method,
    beta     = mr_res$b,
    se       = mr_res$se,
    p        = mr_res$pval,
    p_Qval   = c(h_test$Q_pval[1]),
    p_pleio  = c(p_test$pval)
  ))
h_dat <- harmonise_data(prs3,glu_mr)
mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
                                  'mr_egger_regression',
                                  'mr_weighted_median',
                                  'mr_weighted_mode'))
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- df %>% bind_rows(
  data.frame(
    exposure = 'PRS3',
    outcome  = 'FBG',
    method   = mr_res$method,
    beta     = mr_res$b,
    se       = mr_res$se,
    p        = mr_res$pval,
    p_Qval   = c(h_test$Q_pval[1]),
    p_pleio  = c(p_test$pval)
  ))
h_dat <- harmonise_data(prs4,glu_mr)
mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
                                  'mr_egger_regression',
                                  'mr_weighted_median',
                                  'mr_weighted_mode'))
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- df %>% bind_rows(
  data.frame(
    exposure = 'PRS4',
    outcome  = 'FBG',
    method   = mr_res$method,
    beta     = mr_res$b,
    se       = mr_res$se,
    p        = mr_res$pval,
    p_Qval   = c(h_test$Q_pval[1]),
    p_pleio  = c(p_test$pval)
  ))


# MR outcome: HbA1c.
hba1c <- fread('MAGIC1000G_HbA1c_EUR.tsv.gz')
hba1c2 <- hba1c %>% filter(variant %in% nafld$rsid)

hba1c_mr <- format_data(as.data.frame(hba1c2),
                      type = 'outcome',
                      snp_col =  'variant',
                      beta_col = 'beta',
                      se_col = 'standard_error',
                      eaf_col = 'effect_allele_frequency',
                      effect_allele_col = 'effect_allele',
                      other_allele_col = 'other_allele',
                      pval_col = 'p_value',
                      chr_col = 'chromosome',
                      pos_col = 'base_pair_location')

h_dat <- harmonise_data(nafld_mr,hba1c_mr)
mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
                                  'mr_egger_regression',
                                  'mr_weighted_median',
                                  'mr_weighted_mode'))
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- df %>% bind_rows(
  data.frame(
    exposure = 'NAFLD',
    outcome  = 'hba1c',
    method   = mr_res$method,
    beta     = mr_res$b,
    se       = mr_res$se,
    p        = mr_res$pval,
    p_Qval   = c(h_test$Q_pval[2],NA,NA,NA),
    p_pleio  = c(p_test$pval,NA,NA,NA)
  ))
h_dat <- harmonise_data(prs1,hba1c_mr)
mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
                                  'mr_egger_regression',
                                  'mr_weighted_median',
                                  'mr_weighted_mode'))
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- df %>% bind_rows(
  data.frame(
    exposure = 'PRS1',
    outcome  = 'hba1c',
    method   = mr_res$method,
    beta     = mr_res$b,
    se       = mr_res$se,
    p        = mr_res$pval,
    p_Qval   = c(h_test$Q_pval[1]),
    p_pleio  = c(p_test$pval)
  ))

h_dat <- harmonise_data(prs2,hba1c_mr)
mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
                                  'mr_egger_regression',
                                  'mr_weighted_median',
                                  'mr_weighted_mode'))
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- df %>% bind_rows(
  data.frame(
    exposure = 'PRS2',
    outcome  = 'hba1c',
    method   = mr_res$method,
    beta     = mr_res$b,
    se       = mr_res$se,
    p        = mr_res$pval,
    p_Qval   = c(h_test$Q_pval[1]),
    p_pleio  = c(p_test$pval)
  ))
h_dat <- harmonise_data(prs3,hba1c_mr)
mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
                                  'mr_egger_regression',
                                  'mr_weighted_median',
                                  'mr_weighted_mode'))
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- df %>% bind_rows(
  data.frame(
    exposure = 'PRS3',
    outcome  = 'hba1c',
    method   = mr_res$method,
    beta     = mr_res$b,
    se       = mr_res$se,
    p        = mr_res$pval,
    p_Qval   = c(h_test$Q_pval[1]),
    p_pleio  = c(p_test$pval)
  ))
h_dat <- harmonise_data(prs4,hba1c_mr)
mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
                                  'mr_egger_regression',
                                  'mr_weighted_median',
                                  'mr_weighted_mode'))
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- df %>% bind_rows(
  data.frame(
    exposure = 'PRS4',
    outcome  = 'hba1c',
    method   = mr_res$method,
    beta     = mr_res$b,
    se       = mr_res$se,
    p        = mr_res$pval,
    p_Qval   = c(h_test$Q_pval[1]),
    p_pleio  = c(p_test$pval)
  ))


# MR outcome: fasting insulin.
FI <- fread('MAGIC1000G_FI_EUR.tsv.gz')
FI2 <- FI %>% filter(variant %in% nafld$rsid)

FI_mr <- format_data(as.data.frame(FI2),
                        type = 'outcome',
                        snp_col =  'variant',
                        beta_col = 'beta',
                        se_col = 'standard_error',
                        eaf_col = 'effect_allele_frequency',
                        effect_allele_col = 'effect_allele',
                        other_allele_col = 'other_allele',
                        pval_col = 'p_value',
                        chr_col = 'chromosome',
                        pos_col = 'base_pair_location')

h_dat <- harmonise_data(nafld_mr,FI_mr)
mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
                                  'mr_egger_regression',
                                  'mr_weighted_median',
                                  'mr_weighted_mode'))
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- df %>% bind_rows(
  data.frame(
    exposure = 'NAFLD',
    outcome  = 'FI',
    method   = mr_res$method,
    beta     = mr_res$b,
    se       = mr_res$se,
    p        = mr_res$pval,
    p_Qval   = c(h_test$Q_pval[2],NA,NA,NA),
    p_pleio  = c(p_test$pval,NA,NA,NA)
  ))
h_dat <- harmonise_data(prs1,FI_mr)
mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
                                  'mr_egger_regression',
                                  'mr_weighted_median',
                                  'mr_weighted_mode'))
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- df %>% bind_rows(
  data.frame(
    exposure = 'PRS1',
    outcome  = 'FI',
    method   = mr_res$method,
    beta     = mr_res$b,
    se       = mr_res$se,
    p        = mr_res$pval,
    p_Qval   = c(h_test$Q_pval[1]),
    p_pleio  = c(p_test$pval)
  ))

h_dat <- harmonise_data(prs2,FI_mr)
mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
                                  'mr_egger_regression',
                                  'mr_weighted_median',
                                  'mr_weighted_mode'))
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- df %>% bind_rows(
  data.frame(
    exposure = 'PRS2',
    outcome  = 'FI',
    method   = mr_res$method,
    beta     = mr_res$b,
    se       = mr_res$se,
    p        = mr_res$pval,
    p_Qval   = c(h_test$Q_pval[1]),
    p_pleio  = c(p_test$pval)
  ))
h_dat <- harmonise_data(prs3,FI_mr)
mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
                                  'mr_egger_regression',
                                  'mr_weighted_median',
                                  'mr_weighted_mode'))
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- df %>% bind_rows(
  data.frame(
    exposure = 'PRS3',
    outcome  = 'FI',
    method   = mr_res$method,
    beta     = mr_res$b,
    se       = mr_res$se,
    p        = mr_res$pval,
    p_Qval   = c(h_test$Q_pval[1]),
    p_pleio  = c(p_test$pval)
  ))
h_dat <- harmonise_data(prs4,FI_mr)
mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
                                  'mr_egger_regression',
                                  'mr_weighted_median',
                                  'mr_weighted_mode'))
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- df %>% bind_rows(
  data.frame(
    exposure = 'PRS4',
    outcome  = 'FI',
    method   = mr_res$method,
    beta     = mr_res$b,
    se       = mr_res$se,
    p        = mr_res$pval,
    p_Qval   = c(h_test$Q_pval[1]),
    p_pleio  = c(p_test$pval)
  ))



# MR outcome: hepatocellular carcinoma.
hcc <- fread('finngen_R11_C3_HEPATOCELLU_CARC_EXALLC.gz')
hcc2 <- hcc %>% filter(rsids %in% nafld$rsid)
hcc_mr <- format_data(as.data.frame(hcc2),
                   type = 'outcome',
                   snp_col =  'rsids',
                   beta_col = 'beta',
                   se_col = 'sebeta',
                   eaf_col = 'af_alt',
                   effect_allele_col = 'alt',
                   other_allele_col = 'ref',
                   pval_col = 'pval',
                   chr_col = '#chrom',
                   pos_col = 'pos')
h_dat <- harmonise_data(nafld_mr,hcc_mr)
mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
                                  'mr_egger_regression',
                                  'mr_weighted_median',
                                  'mr_weighted_mode'))
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- df %>% bind_rows(
  data.frame(
    exposure = 'NAFLD',
    outcome  = 'HCC',
    method   = mr_res$method,
    beta     = mr_res$b,
    se       = mr_res$se,
    p        = mr_res$pval,
    p_Qval   = c(h_test$Q_pval[2],NA,NA,NA),
    p_pleio  = c(p_test$pval,NA,NA,NA)
  ))
h_dat <- harmonise_data(prs1,hcc_mr)
mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
                                  'mr_egger_regression',
                                  'mr_weighted_median',
                                  'mr_weighted_mode'))
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- df %>% bind_rows(
  data.frame(
    exposure = 'PRS1',
    outcome  = 'HCC',
    method   = mr_res$method,
    beta     = mr_res$b,
    se       = mr_res$se,
    p        = mr_res$pval,
    p_Qval   = c(h_test$Q_pval[1]),
    p_pleio  = c(p_test$pval)
  ))

h_dat <- harmonise_data(prs2,hcc_mr)
mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
                                  'mr_egger_regression',
                                  'mr_weighted_median',
                                  'mr_weighted_mode'))
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- df %>% bind_rows(
  data.frame(
    exposure = 'PRS2',
    outcome  = 'HCC',
    method   = mr_res$method,
    beta     = mr_res$b,
    se       = mr_res$se,
    p        = mr_res$pval,
    p_Qval   = c(h_test$Q_pval[1]),
    p_pleio  = c(p_test$pval)
  ))
h_dat <- harmonise_data(prs3,hcc_mr)
mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
                                  'mr_egger_regression',
                                  'mr_weighted_median',
                                  'mr_weighted_mode'))
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- df %>% bind_rows(
  data.frame(
    exposure = 'PRS3',
    outcome  = 'HCC',
    method   = mr_res$method,
    beta     = mr_res$b,
    se       = mr_res$se,
    p        = mr_res$pval,
    p_Qval   = c(h_test$Q_pval[1]),
    p_pleio  = c(p_test$pval)
  ))
h_dat <- harmonise_data(prs4,hcc_mr)
mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
                                  'mr_egger_regression',
                                  'mr_weighted_median',
                                  'mr_weighted_mode'))
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- df %>% bind_rows(
  data.frame(
    exposure = 'PRS4',
    outcome  = 'HCC',
    method   = mr_res$method,
    beta     = mr_res$b,
    se       = mr_res$se,
    p        = mr_res$pval,
    p_Qval   = c(h_test$Q_pval[1]),
    p_pleio  = c(p_test$pval)
  ))


# MR outcome: cirrhosis.
cirr <- fread('finngen_R11_CHIRHEP_NAS.gz')
cirr2 <- cirr %>% filter(rsids %in% nafld$rsid)
cirr_mr <- format_data(as.data.frame(cirr2),
                       type = 'outcome',
                       snp_col =  'rsids',
                       beta_col = 'beta',
                       se_col = 'sebeta',
                       eaf_col = 'af_alt',
                       effect_allele_col = 'alt',
                       other_allele_col = 'ref',
                       pval_col = 'pval',
                       chr_col = '#chrom',
                       pos_col = 'pos')
h_dat <- harmonise_data(nafld_mr,cirr_mr)
mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
                                  'mr_egger_regression',
                                  'mr_weighted_median',
                                  'mr_weighted_mode'))
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- df %>% bind_rows(
  data.frame(
    exposure = 'NAFLD',
    outcome  = 'Cirrhosis',
    method   = mr_res$method,
    beta     = mr_res$b,
    se       = mr_res$se,
    p        = mr_res$pval,
    p_Qval   = c(h_test$Q_pval[2],NA,NA,NA),
    p_pleio  = c(p_test$pval,NA,NA,NA)
  ))
h_dat <- harmonise_data(prs1,cirr_mr)
mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
                                  'mr_egger_regression',
                                  'mr_weighted_median',
                                  'mr_weighted_mode'))
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- df %>% bind_rows(
  data.frame(
    exposure = 'PRS1',
    outcome  = 'Cirrhosis',
    method   = mr_res$method,
    beta     = mr_res$b,
    se       = mr_res$se,
    p        = mr_res$pval,
    p_Qval   = c(h_test$Q_pval[1]),
    p_pleio  = c(p_test$pval)
  ))

h_dat <- harmonise_data(prs2,cirr_mr)
mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
                                  'mr_egger_regression',
                                  'mr_weighted_median',
                                  'mr_weighted_mode'))
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- df %>% bind_rows(
  data.frame(
    exposure = 'PRS2',
    outcome  = 'Cirrhosis',
    method   = mr_res$method,
    beta     = mr_res$b,
    se       = mr_res$se,
    p        = mr_res$pval,
    p_Qval   = c(h_test$Q_pval[1]),
    p_pleio  = c(p_test$pval)
  ))
h_dat <- harmonise_data(prs3,cirr_mr)
mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
                                  'mr_egger_regression',
                                  'mr_weighted_median',
                                  'mr_weighted_mode'))
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- df %>% bind_rows(
  data.frame(
    exposure = 'PRS3',
    outcome  = 'Cirrhosis',
    method   = mr_res$method,
    beta     = mr_res$b,
    se       = mr_res$se,
    p        = mr_res$pval,
    p_Qval   = c(h_test$Q_pval[1]),
    p_pleio  = c(p_test$pval)
  ))
h_dat <- harmonise_data(prs4,cirr_mr)
mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
                                  'mr_egger_regression',
                                  'mr_weighted_median',
                                  'mr_weighted_mode'))
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- df %>% bind_rows(
  data.frame(
    exposure = 'PRS4',
    outcome  = 'Cirrhosis',
    method   = mr_res$method,
    beta     = mr_res$b,
    se       = mr_res$se,
    p        = mr_res$pval,
    p_Qval   = c(h_test$Q_pval[1]),
    p_pleio  = c(p_test$pval)
  ))


# MR outcome: insulin resistance (HOMA-IR).
IR <- fread('MAGIC_Manning_et_al_HOMAIR_MainEffect.txt.gz')
IR2 <- IR %>% filter(Snp %in% nafld$rsid)
IR_mr <- format_data(as.data.frame(IR2),
                       type = 'outcome',
                       snp_col =  'Snp',
                       beta_col = 'MainEffects',
                       se_col = 'MainSE',
                       eaf_col = 'maf',
                       effect_allele_col = 'effect_allele',
                       other_allele_col = 'other_allele',
                       pval_col = 'MainP')
h_dat <- harmonise_data(nafld_mr,IR_mr)
mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
                                  'mr_egger_regression',
                                  'mr_weighted_median',
                                  'mr_weighted_mode'))
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- df %>% bind_rows(
  data.frame(
    exposure = 'NAFLD',
    outcome  = 'Insulin resistance',
    method   = mr_res$method,
    beta     = mr_res$b,
    se       = mr_res$se,
    p        = mr_res$pval,
    p_Qval   = c(h_test$Q_pval[2],NA,NA,NA),
    p_pleio  = c(p_test$pval,NA,NA,NA)
  ))
h_dat <- harmonise_data(prs1,IR_mr)
mr_res = mr_wald_ratio(b_exp=h_dat$beta.exposure,
                       b_out=h_dat$beta.outcome,
                       se_exp=h_dat$se.exposure,
                       se_out=h_dat$se.outcome)
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- df %>% bind_rows(
  data.frame(
    exposure = 'PRS1',
    outcome  = 'Insulin resistance',
    method   = mr_res$method,
    beta     = mr_res$b,
    se       = mr_res$se,
    p        = mr_res$pval,
    p_Qval   = c(h_test$Q_pval[1]),
    p_pleio  = c(p_test$pval)
  ))

h_dat <- harmonise_data(prs2,IR_mr)
mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
                                  'mr_egger_regression',
                                  'mr_weighted_median',
                                  'mr_weighted_mode'))
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- df %>% bind_rows(
  data.frame(
    exposure = 'PRS2',
    outcome  = 'Insulin resistance',
    method   = mr_res$method,
    beta     = mr_res$b,
    se       = mr_res$se,
    p        = mr_res$pval,
    p_Qval   = c(h_test$Q_pval[1]),
    p_pleio  = c(p_test$pval)
  ))
h_dat <- harmonise_data(prs3,IR_mr)
mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
                                  'mr_egger_regression',
                                  'mr_weighted_median',
                                  'mr_weighted_mode'))
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- df %>% bind_rows(
  data.frame(
    exposure = 'PRS3',
    outcome  = 'Insulin resistance',
    method   = mr_res$method,
    beta     = mr_res$b,
    se       = mr_res$se,
    p        = mr_res$pval,
    p_Qval   = c(h_test$Q_pval[1]),
    p_pleio  = c(p_test$pval)
  ))
h_dat <- harmonise_data(prs4,IR_mr)
mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
                                  'mr_egger_regression',
                                  'mr_weighted_median',
                                  'mr_weighted_mode'))
h_test = mr_heterogeneity(h_dat)
p_test = mr_pleiotropy_test(h_dat)
df <- df %>% bind_rows(
  data.frame(
    exposure = 'PRS4',
    outcome  = 'Insulin resistance',
    method   = mr_res$method,
    beta     = mr_res$b,
    se       = mr_res$se,
    p        = mr_res$pval,
    p_Qval   = c(h_test$Q_pval[1]),
    p_pleio  = c(p_test$pval)
  ))




# Summarise MR results and plot continuous-trait effects and disease odds ratios.
df$exposure[df$exposure=='NAFLD'] <- 'Overall'
df$outcome[df$outcome=='hba1c'] <- 'HbA1c'
write.csv(df,'MR.csv',row.names = F)

df <- read.csv('MR.csv')
df2 <- df %>% filter(method == 'Inverse variance weighted (fixed effects)',
                     outcome %in% c('HCC','Cirrhosis','T2D'))%>%
  filter(exposure != 'Overall') %>%
  mutate(
  OR = exp(beta),
  lwr= exp(beta-1.96*se),
  upr= exp(beta+1.96*se),
  exposure = factor(exposure,levels = c('PRS4','PRS3','PRS2','PRS1')),
  outcome  = factor(outcome,levels = c('HCC','Cirrhosis','T2D'))
  )

b <- ggplot(df2,aes(OR,outcome,color=exposure))+
  annotate('rect',xmin = rep(0.4,2),xmax = rep(100,2),
           ymin = seq(0.5,3.5,by=2),
           ymax = seq(1.5,4.5,by=2),fill = '#0081cf',
           alpha=.1,color=NA)+
  geom_vline(xintercept = 1,linetype=2)+
  geom_errorbarh(aes(xmin=lwr,xmax=upr),height=0,
                 position = position_dodge(width = 0.75),
                 linewidth=.5)+
  geom_point(size=1,position = position_dodge(width=.75))+
  theme_test()+
  scale_color_manual(name = 'MASLD PRS',
                     values = c('PRS1'='#ED0000','PRS4'='#00468b','PRS2'='#0099b4','PRS3'="#42b540"),
                     breaks = c('PRS1','PRS2','PRS3','PRS4'))+
  scale_x_continuous(limits = c(0.4,100),
                     expand = expansion(add=0),
                     breaks = c(1,2,5,10,20,50),
                     trans = "log10")+
  scale_y_discrete(limits = rev,
                   expand = expansion(add=c(0,0.5)))+
  labs(x = 'Odds ratio (95% CI)',
       y = NULL)+
  theme(axis.text = element_text(size = 11,color='black'),
        axis.title = element_text(size=12))

df3 <- df %>% filter(method == 'Inverse variance weighted (fixed effects)',
                     !outcome %in% c('HCC','Cirrhosis','T2D'))%>%
  filter(exposure != 'Overall') %>%
  mutate(
    lwr= beta-1.96*se,
    upr= beta+1.96*se,
    exposure = factor(exposure,levels = c('PRS4','PRS3','PRS2','PRS1')))
a <- ggplot(df3,aes(beta,outcome,color=exposure))+
  annotate('rect',xmin = rep(-0.3,2),xmax = rep(0.4,2),
           ymin = c(0.5,2.5),
           ymax = c(1.5,3.5),fill = '#0081cf',
           alpha=.1,color=NA)+
  geom_vline(xintercept = 0,linetype=2)+
  geom_errorbarh(aes(xmin=lwr,xmax=upr),height=0,
                 position = position_dodge(width = 0.75),
                 linewidth=.5)+
  geom_point(size=1,position = position_dodge(width=.75))+
  theme_test()+
  scale_color_manual(name = 'MASLD PRS',
                     values = c('PRS1'='#ED0000','PRS4'='#00468b','PRS2'='#0099b4','PRS3'="#42b540"),
                     breaks = c('PRS1','PRS2','PRS3','PRS4'))+
  scale_x_continuous(limits = c(-0.3,0.4),
                     expand = expansion(add=0),
                     breaks = seq(-0.2,0.3,by=0.1))+
  scale_y_discrete(limits = rev,
                   expand = expansion(add=c(0,0.5)))+
  labs(x = 'Effect size (95% CI)',
       y = NULL)+
  theme(axis.text = element_text(size = 11,color='black'),
        axis.title = element_text(size=12))

legend <- get_legend(a)
plot_grid(plot_grid(a+theme(legend.position = 'none'),
                    b+theme(legend.position = 'none'),
                    ncol = 1,align = 'v'),
          legend,ncol = 2,rel_widths = c(1,.2))




# 14. Plot regional GWAS associations around GCKR and rs28601761.
# Retrieve LD information for the locus plots using your own LDlink token.
####colocalization for rs1260326
library(EnsDb.Hsapiens.v75)
library(locuszoomr)
pdff <- fread('F:\\English Papers for writing\\PDFF-GWAS-MR\\pdff_ukbb_31377.gz')
chr2 <- pdff %>% dplyr::filter(chr==2,pos > 27730940-1024*1024,
                        pos<27730940+1024*1024)
chr2 <- chr2 %>% dplyr::select(chr,pos,snp,beta,se,pval) %>% as.data.frame() 
loc <- locus(data = chr2, index_snp =  'rs1260326', flank = 1e5,
             LD = 'r2',
             ens_db = "EnsDb.Hsapiens.v75")
loc <- link_LD(loc, token = "YOUR_LDLINK_API_TOKEN")
summary(loc)
locus_plot(loc, labels = c("index"),
           label_x = c(4, -5))

t2d <- fread('Mahajan.NG2018b.T2D.European.txt.gz')
chr2 <- t2d %>% dplyr::filter(Chr==2,Pos > 27730940-1024*1024,
                               Pos<27730940+1024*1024)
chr2 <- chr2 %>% dplyr::select(Chr,Pos,name,Beta,SE,Pvalue) %>%
  dplyr::rename(SNP=name) %>% as.data.frame() 
chr2$Pvalue <- as.numeric(chr2$Pvalue) 
loc <- locus(data = chr2, index_snp =  'rs1260326', flank = 1e5,
             LD = 'r2',
             ens_db = "EnsDb.Hsapiens.v75")
loc <- link_LD(loc, token = "YOUR_LDLINK_API_TOKEN")
summary(loc)
locus_plot(loc, labels = c("index"),
           label_x = c(4, -5))


chr2 <- pdff %>% dplyr::filter(chr==8,pos > 126500031-1024*1024,
                               pos<126500031+1024*1024)
chr2 <- chr2 %>% dplyr::select(chr,pos,snp,beta,se,pval) %>% as.data.frame() 
loc <- locus(data = chr2, index_snp =  'rs28601761', flank = 1e5,
             LD = 'r2',
             ens_db = "EnsDb.Hsapiens.v75")
loc <- link_LD(loc, token = "YOUR_LDLINK_API_TOKEN")
summary(loc)
locus_plot(loc, labels = c("index"),
           label_x = c(4, -5))

chr2 <- t2d %>% dplyr::filter(Chr==8,Pos > 126500031-1024*1024,
                              Pos<126500031+1024*1024)
chr2 <- chr2 %>% dplyr::select(Chr,Pos,name,Beta,SE,Pvalue) %>%
  dplyr::rename(SNP=name) %>% as.data.frame() 
chr2$Pvalue <- as.numeric(chr2$Pvalue) 
loc <- locus(data = chr2, index_snp =  'rs28601761', flank = 1e5,
             LD = 'r2',
             ens_db = "EnsDb.Hsapiens.v75")
loc <- link_LD(loc, token = "YOUR_LDLINK_API_TOKEN")
summary(loc)
locus_plot(loc, labels = c("index"),
           label_x = c(4, -5))




# 15. Annotate and plot rs1260326 associations with circulating proteins.
###volcano plot for proteins(Pietzner*, Wheeler* et al., Science, 2021)
pros <- read.csv('rs1260326-Protein-omicsicence.csv')
somalogic <- readxl::read_excel('SomalogicProteinIDs.xlsx')
pros <- pros %>% mutate(
  SeqID = gsub('-','_',som_seqid),
  logP = -log10(Pvalue)
)
pros <- pros %>% dplyr::left_join(somalogic %>% 
                                    dplyr::select(SeqId,`Protein (short name)`,
                                                  Gene,Ensembl.Gene.ID),
                                  by= c('SeqID'='SeqId'))

ggplot(pros,aes(Effect,logP,color=strata(Pvalue<0.05/1590)))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05/1590),linetype=2)+
  theme_classic()+
  scale_color_jama(name = 'Pvalue',labels = c('>0.05/1590','<0.05/1590'))+
  labs(x = 'Effect size',y = '-log10(P)')+
  #scale_y_continuous(limits = c(0,40))+
  ggrepel::geom_text_repel(data=pros %>% filter(Gene %in% c('GCKR','GCK')),
                           aes(Effect,logP,label=`Gene`),
                           color='black',
                           max.overlaps = 30)



# 16. Inspect protein QTLs, select/clump candidate instruments and plot effects.
###MR GCKR-->T2D
gckr_decode <- fread('5223_59_GCKR_GCKR.txt.gz')
gckr_decode %>% filter(rsids == 'rs1260326')

gckr_decode %>% filter(Chrom == 'chr2') %>%
  dplyr::slice(which.min(Pval))
#rs1260326 C>T beta=-0.0866,se=0.0084, p = 5.745e-25
#lead cis-pQTL, PIP = ?

pqtlGCKR <- gckr_decode %>% filter(Pval<5e-8,
                                   Pos < 27496839+1024*1024,
                                   Pos > 27496839-1024*1024)
library(TwoSampleMR)
library(ieugwasr)
x = ld_clump(data.frame(rsid=pqtlGCKR$rsids,pval=pqtlGCKR$Pval),
             clump_r2=0.01,
             plink_bin='F:\\E\\plink\\plink',
             bfile='F:\\E\\EUR-LD\\EUR')
pro_mr <- pqtlGCKR %>% filter(rsids %in% x$rsid)



gck <- fread('12960_9_GCK_HXK4.txt.gz')
gck %>% filter(rsids == 'rs1260326')
#rs1260326 C>T beta=0.0427,se=0.0085, p = 4.892e-7
snx17 <- fread('12845_18_SNX17_SNX17.txt.gz')
snx17 %>% filter(rsids == 'rs1260326')
#rs1260326 C>T beta=0.0288,se=0.0086, p = 7.75e-4


df <- data.frame(
  study = rep(c('Decode','Omicscience'),each=2),
  protein=rep(c('GCKR','GCK'),2),
  beta= c(-0.0866,0.0427,-0.322,0.050),
  se  = c(0.0084,0.0085,0.013,0.013),
  p   = c(5.75e-25,4.89e-7,7.31e-128,1.95e-4)
)
library(scales)
# 绘制误差条图
a <- ggplot(df, aes(x = beta, y = protein, color = study, group = study)) +
  geom_vline(xintercept = 0,linetype=2)+
  geom_errorbarh(aes(xmin = beta - se, xmax = beta + se), height = 0.05, 
                 size = 0.8, position = position_dodge(width = 0.3)) +
  geom_point(size = 3, position = position_dodge(width = 0.3)) +
  geom_text(aes(label = paste0("p=", formatC(p, format = "e", digits = 2))), 
            hjust = 0.5, vjust = -1.2, size = 4, 
            position = position_dodge(width = 0.3),show.legend = F) +
  scale_color_brewer(palette = "Set1") +
   labs(
    x = "Effect size (95% CI)",
    y = ""
  ) +
  theme_test() +
  theme(legend.position = 'inside',
        legend.position.inside = c(.2,.2),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )


# 17. Import T2D fine-mapping results and examine rs1260326-associated T2D risk.
###fine mapping of rs1260326 on T2D
fm_t2d <- fread('fineMappingT2D.tsv')
fm_t2d <- fm_t2d %>% dplyr::filter(grepl('Type 2 diabetes',x = reportedTrait))
write.csv(fm_t2d,'rs1260326FineMappingT2D.csv',row.names = F)

###effects of rs1260326 on T2D risk
ukb6 <- ukb6 %>% mutate(
  rs1260326G = factor(rs1260326,labels = c('p.P446P','p.P446L','p.L446L'))
)
model1 <- coxph(Surv(followUp_days_T2D,T2D)~rs1260326G+sex+age,
                ukb6)
summary(model1)
# rs1260326Gp.P446L    0.9408     1.0629    0.9158    0.9665
# rs1260326Gp.L446L    0.8913     1.1220    0.8580    0.9258
new_df <- with(ukb6,
               data.frame(rs1260326G = c('p.P446P','p.P446L','p.L446L'),
                          sex = 1,
                          age = 60
               ))

fit <- survfit(model1, newdata = new_df)
dat = data.frame(surv = fit$surv,lower= fit$lower,
                 upper = fit$upper,time= fit$time)
data = melt(dat, id="time")
data_wide <- data %>%
  separate(col = variable, into = c("type", "pPRS1"), sep = "\\.") %>%
  spread(key = type, value = value)

b <- ggplot(data = data_wide) +
  geom_line(aes(x = time, y = 1-surv, group = pPRS1,
                colour = pPRS1),
            size = 1.2,
            linetype = 1) +
  geom_ribbon(aes(x = time, ymin = 1-lower, ymax = 1-upper,
                  group = pPRS1, fill = pPRS1),
              alpha = 0.5,show.legend = F) +
  ggsci::scale_colour_npg(name ='GCKR rs1260326',
                          labels = c('p.P446P','p.P446L','p.L446L')) +
  ggsci::scale_fill_npg()+
  #scale_fill_manual(values = c('#112d4e','#3f72af')) +
  scale_x_continuous(limits = c(0, 5475),
                     breaks = seq(0, 5475, 365*3),
                     labels = seq(0,15,by=3)) +
  theme_bw()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = 'inside',
        legend.position.inside = c(.25,.8),
        legend.text = element_text(size=11)) +
  ylab("Cumulative probability")+
  xlab('Years of follow-up')+
  scale_y_continuous(labels = scales::percent,
                     limits = c(0,0.1))
data_wide %>% filter(time == 365*15-1)
1-c(0.9034488,0.9088933,0.9134797)
plot_grid(a,b,ncol = 2)

model1 <- coxph(Surv(followUp_days_T2D,T2D)~rs1260326+sex+age,
                ukb6)
summary(model1)
#6% decrease in T2D risk per allele (0.5% absolute decrease)


# 18. Compare imported fine-mapping posterior probabilities across resources.
###Fine mapping


fm_t2d <- fread('fineMappingT2D.tsv')
fm_nafld <- fm_t2d %>% dplyr::filter(grepl('(fatty liver)|(Liver proton density fat fraction)',x = reportedTrait))
fm_nafld <- fm_nafld %>% mutate(
  study = factor(study,levels = study[order(posteriorProbability,decreasing = T)])
)

fm_t2d <- read.csv('rs1260326FineMappingT2D.csv')
fm_t2d <- fm_t2d %>% filter(leadVariant == '2_27508073_T_C')
fm_t2d <- fm_t2d %>% mutate(
  study = factor(study,levels = study[order(posteriorProbability,decreasing = T)])
)

df <- data.frame(
  traits = c(rep(c('GCKR','MASLD','T2D'),c(2,4,12))),
  study  = c('Decode','Omicscience',as.character(fm_nafld$study)
             ,as.character(fm_t2d$study)),
  PIP    = c(0.871,1,fm_nafld$posteriorProbability,fm_t2d$posteriorProbability)
)
df <- df %>% mutate(
  traits = factor(traits,levels = rev(c('GCKR','MASLD','T2D'))),
  study = factor(study,levels = rev(study[order(traits,PIP,decreasing = F)]))
)

ggplot(df,aes(PIP,study,color=traits))+
  geom_linerange(aes(xmin=0,xmax=PIP))+
  geom_point(size=2)+
  scale_color_aaas(name='Trait',
                  breaks = c('GCKR','MASLD','T2D'))+
  theme_minimal()+
  scale_y_discrete(limits=rev)+
  scale_x_continuous(expand = expansion(add=c(0,0.02)),
                     limits = c(0,1))+
  theme(axis.text = element_text(size=11,color='black'),
        legend.position = 'top')+
  labs(x = 'Posterior inclusion probability',
       y = NULL)


# 19. Extract reference genomic sequences for the GCKR region and export FASTA.
# 安装人类基因组
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)
genome <- BSgenome.Hsapiens.UCSC.hg38

# 示例：提取 chr2 上某个区域（如 GCKR 第9个外显子）
seq1 <- getSeq(genome, names = "chr2", start = 27507976, end = 27508250)
seq1
seq2 <- getSeq(genome, names = "chr2", start = 27507680, end = 27508250)
seq2
library(Biostrings)

# 包装成 DNAStringSet
seq_set <- DNAStringSet(seq1)

# 添加名称
names(seq_set) <- "GCKR_exon15_intron15_exon16"

# 导出为 FASTA
writeXStringSet(seq_set, filepath = "GCKR_exon15-16_region.fasta")

seq_set <- DNAStringSet(seq2)

# 添加名称
names(seq_set) <- "GCKR_exon14_exon16"

# 导出为 FASTA
writeXStringSet(seq_set, filepath = "GCKR_exon14-16_region.fasta")
