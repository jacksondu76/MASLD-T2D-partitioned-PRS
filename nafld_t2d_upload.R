###NAFLD AND T2D####
setwd('F:\\English Papers for writing\\nafld-t2d')
library(data.table)
library(tidyverse)
library(cowplot)
library(survival)
library(ggsci)

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

ukb <- fread('F:\\English Papers for writing\\have done\\ALCO-MAFLD\\ukb_alco_mafld.csv')
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

## REDACTED FOR CODE SHARING: the UK Biobank field mapping, data-cleaning rules,
## derivation of covariates, and cohort eligibility criteria are not included.
## The analyses below require an internally prepared `ukb5` object with the
## documented analysis variables and genotypes.
##
## ukb5 <- readRDS("path/to/internally_prepared_analysis_cohort.rds")
## Implementation intentionally omitted in the shared version.
metabolites <- fread('F:\\E\\UkbData\\UKB-metabolites-updated.csv')
ukb_eid_match <- read.csv('F:\\E\\UkbData\\ukb_match_63726_92718.csv',sep = '\t')
metabolites <- metabolites %>% left_join(ukb_eid_match,by = c('eid'='eid.92718'))
metabolites <- metabolites %>% filter(!is.na(eid.63726))
ukb5 <- ukb5 %>% left_join(metabolites,by = c('eid'='eid.63726'))

ukb_mets <- readxl::read_excel('F:\\E\\UkbData\\ukb_mets.xlsx')
ukb_mets <- ukb_mets %>% mutate(
  metID= paste0('f.',eid,'.0.0')
) %>% rename(ID = eid)

tg <- ukb_mets$metID[grep('Triglyceride',ukb_mets$metabolites)]
rsids <- names(snps)[5:22]

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

df2 <- read.csv('snps_tg.csv')
snp_mets <- df2 %>% mutate(
  fdr = p.adjust(p,'fdr'),
  beta2 = ifelse(fdr<0.05,beta,0)
)


# Load necessary libraries
# 加载必要的包
library(factoextra)  # K-means聚类可视化
library(cluster)     # 轮廓系数
library(Rtsne)       # t-SNE降维
library(patchwork)   # 图形拼接
library(pheatmap)    # 热图

# 读取数据
data <- snp_mets

# 数据预处理：提取beta值构建矩阵
df <- data %>%
  select(snp, metabolites, beta) %>%
  spread(key = metabolites, value = beta, fill = 0)

# 设置SNP为行名
rownames(df) <- df$snp
df <- df[, -1]

# 标准化数据
df_scaled <- scale(df)

# 确定最佳聚类数 (A) - Elbow Plot
elbow_plot <- fviz_nbclust(df_scaled, kmeans, method = "wss",)

# 执行K-means聚类 (假设K=4)
set.seed(42)
k <- 4
kmeans_result <- kmeans(df_scaled, centers = k, nstart = 25)
kmeans_result$cluster <- c(2,3,1,2,2,3,3,3,3,3,1,3,4,4,3,3,3,3)
# 聚类效果 (B) - PCA可视化
pca_plot <- fviz_cluster(kmeans_result, data = df_scaled, 
                         geom = "point", ellipse.type = "convex", 
                         palette = c('#ED0000','#0099b4',"#42b540",'#00468b'),
                         ggtheme = theme_minimal())

# 使用t-SNE进一步可视化 (C)
set.seed(518)
tsne <- Rtsne(df_scaled, dims = 2, perplexity = 3)
tsne_df <- data.frame(tsne$Y, cluster = as.factor(kmeans_result$cluster))
tsne_plot <- ggplot(tsne_df, aes(x = X1, y = X2, color = cluster)) +
  scale_color_manual(values = c('#ED0000','#0099b4',"#42b540",'#00468b'))+
  geom_point(size = 3) +
  theme_minimal()

# 轮廓系数分析 (D)
# sil <- silhouette(kmeans_result$cluster, dist(df_scaled))
# sil_plot <- fviz_silhouette(sil) +
#   ggtitle("D. Silhouette Plot of Clusters")

# 聚类中心热图 (E)
centers <- as.data.frame(kmeans_result$centers)
color_palette <- colorRampPalette(c("#0099b4", "white", "#ED0000"))(50)  # 红蓝对比

# 绘制聚类中心热图
pheatmap(
  centers,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "row",  # 标准化每一行
  color = color_palette,  # 自定义颜色
  border_color = NA,  # 去除网格线，提升美观度
  #main = "Cluster Centers Heatmap (Optimized Colors)",
  fontsize_row = 10,
  fontsize_col = 10,
  cellwidth = 25,
  cellheight = 15,
  treeheight_row = 50,
  treeheight_col = 50,
  angle_col = 45
)
# 组合成一张图
(elbow_plot | pca_plot|tsne_plot )



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

## REDACTED FOR CODE SHARING: PRS composition (variants and weights) is retained
## only in the internal analysis workflow. The shared cohort must already contain
## the derived PRS variables used below.
## Implementation intentionally omitted in the shared version.
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


table(ukb6$T2D)
# 0      1 
# 390051 25116 

table(ukb6$HCC)
# 0      1 
# 414866    301

table(ukb6$Cirrhosis)
# 0      1 
# 411236   3931

table(ukb6$MASLD)
# 0      1 
# 409674   5493  

ukb6 <- ukb6 %>% mutate(
  loweringINS = ifelse(((sex==1 & medCholBPDiaM==3) |(sex==0 & medCholBPDiaF==3)),1,0)
  
)
fwrite(ukb6,'analyzingData.csv')

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

# prs1                     0.7834        0.6534    0.9393
# prs234                   1.3033        1.2411    1.3687

##adjusting IR


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

# ###association between MASLD and T2D
# ukb6 <- ukb6 %>% mutate(
#   TRG2 = TRG/0.0113,
#   FLI = exp(0.953*log(TRG2)+0.139*BMI+0.718*log(GGT)+
#               0.053*waistCir-15.745)/(1+exp(0.953*log(TRG2)+0.139*BMI+0.718*log(GGT)+
#                                               0.053*waistCir-15.745))*100,
#   fld = ifelse(FLI>=60,1,0)
# )
# 
# ###alcohol drinking
# #alcohol intake
# #g/day
# #lifetime nondrinkers
# ukb6_1 <-  ukb6 %>% filter(alcoholStatus==0) %>%
#   mutate(alcoholIntake = 0) %>%
#   mutate(redwineIntake = 0,
#          whitewineIntake = 0,
#          beerIntake = 0,
#          spiritsIntake =0,
#          fortifiedwineIntake=0,
#          otherwineIntake=0)
# 
# #drink monthly
# ukb6_2 <- ukb6 %>% filter(aolIntakeFreq>0 & aolIntakeFreq<6,
#                           !is.na(redwine2)) %>%
#   mutate(redwine2 = ifelse(redwine2<0,0,redwine2),
#          whitewine2 = ifelse(whitewine2<0,0,whitewine2),
#          beer2 = ifelse(beer2<0,0,beer2),
#          spirits2 = ifelse(spirits2<0,0,spirits2),
#          fortifiedwine2 = ifelse(fortifiedwine2<0,0,fortifiedwine2),
#          otherwine2 = ifelse(otherwine2<0,0,otherwine2)) %>%
#   mutate(alcoholIntake = (redwine2*13/100*125+whitewine2*12/100*125+
#                             beer2*4.5/100*568.3+spirits2*40/100*25+
#                             fortifiedwine2*18/100*75+otherwine2*4/100*275)/30) %>%
#   mutate(
#     redwineIntake = redwine2*13/100*125/30,
#     whitewineIntake = whitewine2*12/100*125/30,
#     beerIntake = beer2*4.5/100*568.3/30,
#     spiritsIntake = spirits2*40/100*25/30,
#     fortifiedwineIntake=fortifiedwine2*18/100*75/30,
#     otherwineIntake=otherwine2*4/100*275/30
#   ) %>% filter(alcoholIntake>0)
# 
# ##drink weekly
# ukb6_3 <- ukb6 %>% filter(aolIntakeFreq>0 & aolIntakeFreq<6,
#                           !is.na(redwine1)) %>% 
#   mutate(redwine1 = ifelse(redwine1<0,0,redwine1),
#          whitewine1 = ifelse(whitewine1<0,0,whitewine1),
#          beer1 = ifelse(beer1<0,0,beer1),
#          spirits1 = ifelse(spirits1<0,0,spirits1),
#          fortifiedwine1 = ifelse(fortifiedwine1<0,0,fortifiedwine1),
#          otherwine1 = ifelse(is.na(otherwine1)|otherwine1<0,0,otherwine1)) %>%
#   mutate(alcoholIntake = (redwine1*13/100*125+whitewine1*12/100*125+
#                             beer1*4.5/100*568.3+spirits1*40/100*25+
#                             fortifiedwine1*18/100*75+otherwine1*4/100*275)/7) %>%
#   filter(!is.na(alcoholIntake)) %>%
#   mutate(
#     redwineIntake = redwine1*13/100*125/7,
#     whitewineIntake = whitewine1*12/100*125/7,
#     beerIntake = beer1*4.5/100*568.3/7,
#     spiritsIntake = spirits1*40/100*25/7,
#     fortifiedwineIntake=fortifiedwine1*18/100*75/7,
#     otherwineIntake=otherwine1*4/100*275/7
#   )%>% filter(alcoholIntake>0)
# 
# ukb6_al <- bind_rows(ukb6_1 %>% select(-c(pdff2,HbA1c2,Glucose2,TRG2)),
#                      ukb6_2 %>% select(-c(pdff2,HbA1c2,Glucose2,TRG2)),
#                      ukb6_3 %>% select(-c(pdff2,HbA1c2,Glucose2,TRG2)))
# 
# ukb7 <- ukb6 %>% left_join(ukb6_al %>% select(eid,alcoholIntake))
# ukb7 <- ukb7 %>% mutate(
#   loweringCHL = ifelse(((sex==1 & medCholBPDiaM==1) |(sex==0 & medCholBPDiaF==1)),1,0),
#   #loweringINS = ifelse(((sex==1 & medCholBPDiaM==3) |(sex==0 & medCholBPDiaF==3)),1,0),
#   loweringBP  = ifelse(((sex==1 & medCholBPDiaM==2) |(sex==0 & medCholBPDiaF==2)),1,0)
# )
# 
# ukb7 <- ukb7 %>% mutate(
#   cr1 = ifelse(BMI>=25|(waistCir>94 & sex==1) | (waistCir>80 & sex==0),1,0),
#   cr2 = ifelse(Glucose>=5.6 | HbA1c>=39 | loweringINS == 1,1,0),
#   cr3 = ifelse((SBP >= 130 & DBP >= 85)|loweringBP==1,1,0),
#   cr4 = ifelse(TRG >= 1.7 | loweringCHL == 1,1,0),
#   cr5 = ifelse((HDL<=1 & sex==1)|(HDL<=1.3 & sex==0)|loweringCHL==1,1,0),
#   cardio = cr1+cr2+cr3+cr4+cr5,
#   masld  = ifelse(fld==1 & cardio>0 & ((alcoholIntake<=30 & sex==1)|(alcoholIntake<=20 & sex==0)),1,
#                   ifelse(fld==0,0,2))
# )
# table(ukb7$masld)
# 
# 
# ukb7 <- ukb7 %>% mutate(
#   prs1G = ifelse(prs1<0.067,'G1',ifelse(prs1>0.211,'G3','G2')),
#   prs234G = ifelse(prs234<0.884,'G1',ifelse(prs234>1.510,'G3','G2'))
# )
# ukb8 <- ukb7 %>% filter(masld != 2) 
# table(ukb8$prs1G)
# table(ukb8$prs234G)
# 
# 
# fit_cox1 <- coxph(Surv(followUp_days_T2D,T2D)~sex+age+smokingStatus+averageTotalIncome+
#                     aolIntakeFreq+no.daysModeAct2+education+loweringINS+
#                     masld*prs1,ukb8)
# summary(fit_cox1)
# #p=0.05
# 
# fit_cox1 <- coxph(Surv(followUp_days_T2D,T2D)~sex+age+smokingStatus+averageTotalIncome+
#                     aolIntakeFreq+no.daysModeAct2+education+loweringINS+
#                     masld*prs234,ukb8)
# summary(fit_cox1)
# #p=0.5;###no interaction effect
# 
# fit_cox1 <- coxph(Surv(followUp_days_T2D,T2D)~sex+age+smokingStatus+averageTotalIncome+
#                     aolIntakeFreq+no.daysModeAct2+education+loweringINS+
#                     prs234+
#                     masld,ukb8,subset = prs1G=='G2')
# summary(fit_cox1)
# #4.79 (4.19-5.47);4.71(4.55-4.89);4.75(4.05-5.57) 
# fit_cox1 <- coxph(Surv(followUp_days_T2D,T2D)~sex+age+smokingStatus+averageTotalIncome+
#                     aolIntakeFreq+no.daysModeAct2+education+loweringINS+
#                     masld,ukb8,subset = prs234G=='G3')
# summary(fit_cox1)
# #3.66(3.06-4.38);3.65(3.21-4.14);3.96(3.34-4.69)

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

####T2D subtypes
# t2d_s1 <- fread('GCST90026412_buildGRCh37.tsv.gz')
# t2d_s1 <- t2d_s1 %>% filter(variant_ID %in% nafld$rsid)
# t2ds1_mr <- format_data(as.data.frame(t2d_s1),
#                        type = 'outcome',
#                        snp_col =  'variant_ID',
#                        beta_col = 'beta',
#                        se_col = 'standard_error',
#                        eaf_col = 'effect_allele_frequency',
#                        effect_allele_col = 'effect_allele',
#                        other_allele_col = 'other_allele',
#                        pval_col = 'p_value',
#                        chr_col = 'chromosome',
#                        pos_col = 'base_pair_location')
# h_dat <- harmonise_data(nafld_mr,t2ds1_mr)
# h_dat <- h_dat %>% mutate(mr_keep=T)
# mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
#                                   'mr_egger_regression',
#                                   'mr_weighted_median',
#                                   'mr_weighted_mode'))
# h_test = mr_heterogeneity(h_dat)
# p_test = mr_pleiotropy_test(h_dat)
# df <- data.frame(
#     exposure = 'Overall',
#     outcome  = 'Severe autoimmune T2D',
#     method   = mr_res$method,
#     beta     = mr_res$b,
#     se       = mr_res$se,
#     p        = mr_res$pval,
#     p_Qval   = c(h_test$Q_pval[2],NA,NA,NA),
#     p_pleio  = c(p_test$pval,NA,NA,NA)
#   )
# h_dat <- harmonise_data(prs1,t2ds1_mr)
# h_dat <- h_dat %>% mutate(mr_keep=T)
# mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
#                                   'mr_egger_regression',
#                                   'mr_weighted_median',
#                                   'mr_weighted_mode'))
# h_test = mr_heterogeneity(h_dat)
# p_test = mr_pleiotropy_test(h_dat)
# df <- df %>% bind_rows(
#   data.frame(
#   exposure = 'PRS1',
#   outcome  = 'Severe autoimmune T2D',
#   method   = mr_res$method,
#   beta     = mr_res$b,
#   se       = mr_res$se,
#   p        = mr_res$pval,
#   p_Qval   = c(h_test$Q_pval[2]),
#   p_pleio  = c(p_test$pval)
# ))
# 
# h_dat <- harmonise_data(prs2,t2ds1_mr)
# h_dat <- h_dat %>% mutate(mr_keep=T)
# mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
#                                   'mr_egger_regression',
#                                   'mr_weighted_median',
#                                   'mr_weighted_mode'))
# 
# h_test = mr_heterogeneity(h_dat)
# p_test = mr_pleiotropy_test(h_dat)
# df <- df %>% bind_rows(
#   data.frame(
#     exposure = 'PRS2',
#     outcome  = 'Severe autoimmune T2D',
#     method   = mr_res$method,
#     beta     = mr_res$b,
#     se       = mr_res$se,
#     p        = mr_res$pval,
#     p_Qval   = c(h_test$Q_pval[2]),
#     p_pleio  = c(p_test$pval)
#   ))
# 
# h_dat <- harmonise_data(prs3,t2ds1_mr)
# h_dat <- h_dat %>% mutate(mr_keep=T)
# mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
#                                   'mr_egger_regression',
#                                   'mr_weighted_median',
#                                   'mr_weighted_mode'))
# 
# h_test = mr_heterogeneity(h_dat)
# p_test = mr_pleiotropy_test(h_dat)
# df <- df %>% bind_rows(
#   data.frame(
#     exposure = 'PRS3',
#     outcome  = 'Severe autoimmune T2D',
#     method   = mr_res$method,
#     beta     = mr_res$b,
#     se       = mr_res$se,
#     p        = mr_res$pval,
#     p_Qval   = c(h_test$Q_pval[2]),
#     p_pleio  = c(p_test$pval)
#   ))
# 
# h_dat <- harmonise_data(prs4,t2ds1_mr)
# h_dat <- h_dat %>% mutate(mr_keep=T)
# mr_res = mr_wald_ratio(b_exp = h_dat$beta.exposure,
#                        b_out = h_dat$beta.outcome,
#                        se_exp = h_dat$se.exposure,
#                        se_out = h_dat$se.outcome)
# 
# df <- df %>% bind_rows(
#   data.frame(
#     exposure = 'PRS4',
#     outcome  = 'Severe autoimmune T2D',
#     method   = 'Wald ratio',
#     beta     = mr_res$b,
#     se       = mr_res$se,
#     p        = mr_res$pval,
#     p_Qval   = NA,
#     p_pleio  = NA
#   ))
# 
# t2d_s2 <- fread('GCST90026413_buildGRCh37.tsv.gz')
# t2d_s2 <- t2d_s2 %>% filter(variant_ID %in% nafld$rsid)
# t2ds2_mr <- format_data(as.data.frame(t2d_s2),
#                         type = 'outcome',
#                         snp_col =  'variant_ID',
#                         beta_col = 'beta',
#                         se_col = 'standard_error',
#                         eaf_col = 'effect_allele_frequency',
#                         effect_allele_col = 'effect_allele',
#                         other_allele_col = 'other_allele',
#                         pval_col = 'p_value',
#                         chr_col = 'chromosome',
#                         pos_col = 'base_pair_location')
# h_dat <- harmonise_data(nafld_mr,t2ds2_mr)
# h_dat <- h_dat %>% mutate(mr_keep=T)
# mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
#                                   'mr_egger_regression',
#                                   'mr_weighted_median',
#                                   'mr_weighted_mode'))
# h_test = mr_heterogeneity(h_dat)
# p_test = mr_pleiotropy_test(h_dat)
# df <- df %>% bind_rows(data.frame(
#   exposure = 'Overall',
#   outcome  = 'Severe insulin-deficient T2D',
#   method   = mr_res$method,
#   beta     = mr_res$b,
#   se       = mr_res$se,
#   p        = mr_res$pval,
#   p_Qval   = c(h_test$Q_pval[2],NA,NA,NA),
#   p_pleio  = c(p_test$pval,NA,NA,NA)
# ))
# h_dat <- harmonise_data(prs1,t2ds2_mr)
# h_dat <- h_dat %>% mutate(mr_keep=T)
# mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
#                                   'mr_egger_regression',
#                                   'mr_weighted_median',
#                                   'mr_weighted_mode'))
# h_test = mr_heterogeneity(h_dat)
# p_test = mr_pleiotropy_test(h_dat)
# df <- df %>% bind_rows(
#   data.frame(
#     exposure = 'PRS1',
#     outcome  = 'Severe insulin-deficient T2D',
#     method   = mr_res$method,
#     beta     = mr_res$b,
#     se       = mr_res$se,
#     p        = mr_res$pval,
#     p_Qval   = c(h_test$Q_pval[2]),
#     p_pleio  = c(p_test$pval)
#   ))
# 
# h_dat <- harmonise_data(prs2,t2ds2_mr)
# h_dat <- h_dat %>% mutate(mr_keep=T)
# mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
#                                   'mr_egger_regression',
#                                   'mr_weighted_median',
#                                   'mr_weighted_mode'))
# 
# h_test = mr_heterogeneity(h_dat)
# p_test = mr_pleiotropy_test(h_dat)
# df <- df %>% bind_rows(
#   data.frame(
#     exposure = 'PRS2',
#     outcome  = 'Severe insulin-deficient T2D',
#     method   = mr_res$method,
#     beta     = mr_res$b,
#     se       = mr_res$se,
#     p        = mr_res$pval,
#     p_Qval   = c(h_test$Q_pval[2]),
#     p_pleio  = c(p_test$pval)
#   ))
# 
# h_dat <- harmonise_data(prs3,t2ds2_mr)
# h_dat <- h_dat %>% mutate(mr_keep=T)
# mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
#                                   'mr_egger_regression',
#                                   'mr_weighted_median',
#                                   'mr_weighted_mode'))
# 
# h_test = mr_heterogeneity(h_dat)
# p_test = mr_pleiotropy_test(h_dat)
# df <- df %>% bind_rows(
#   data.frame(
#     exposure = 'PRS3',
#     outcome  = 'Severe insulin-deficient T2D',
#     method   = mr_res$method,
#     beta     = mr_res$b,
#     se       = mr_res$se,
#     p        = mr_res$pval,
#     p_Qval   = c(h_test$Q_pval[2]),
#     p_pleio  = c(p_test$pval)
#   ))
# 
# h_dat <- harmonise_data(prs4,t2ds2_mr)
# h_dat <- h_dat %>% mutate(mr_keep=T)
# mr_res = mr_wald_ratio(b_exp = h_dat$beta.exposure,
#                        b_out = h_dat$beta.outcome,
#                        se_exp = h_dat$se.exposure,
#                        se_out = h_dat$se.outcome)
# 
# df <- df %>% bind_rows(
#   data.frame(
#     exposure = 'PRS4',
#     outcome  = 'Severe insulin-deficient T2D',
#     method   = 'Wald ratio',
#     beta     = mr_res$b,
#     se       = mr_res$se,
#     p        = mr_res$pval,
#     p_Qval   = NA,
#     p_pleio  = NA
#   ))
# 
# 
# t2d_s3 <- fread('GCST90026414_buildGRCh37.tsv.gz')
# t2d_s3 <- t2d_s3 %>% filter(variant_ID %in% nafld$rsid)
# t2ds3_mr <- format_data(as.data.frame(t2d_s3),
#                         type = 'outcome',
#                         snp_col =  'variant_ID',
#                         beta_col = 'beta',
#                         se_col = 'standard_error',
#                         eaf_col = 'effect_allele_frequency',
#                         effect_allele_col = 'effect_allele',
#                         other_allele_col = 'other_allele',
#                         pval_col = 'p_value',
#                         chr_col = 'chromosome',
#                         pos_col = 'base_pair_location')
# h_dat <- harmonise_data(nafld_mr,t2ds3_mr)
# h_dat <- h_dat %>% mutate(mr_keep=T)
# mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
#                                   'mr_egger_regression',
#                                   'mr_weighted_median',
#                                   'mr_weighted_mode'))
# h_test = mr_heterogeneity(h_dat)
# p_test = mr_pleiotropy_test(h_dat)
# df <- df %>% bind_rows(data.frame(
#   exposure = 'Overall',
#   outcome  = 'Severe insulin-resistant T2D',
#   method   = mr_res$method,
#   beta     = mr_res$b,
#   se       = mr_res$se,
#   p        = mr_res$pval,
#   p_Qval   = c(h_test$Q_pval[2],NA,NA,NA),
#   p_pleio  = c(p_test$pval,NA,NA,NA)
# ))
# h_dat <- harmonise_data(prs1,t2ds3_mr)
# h_dat <- h_dat %>% mutate(mr_keep=T)
# mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
#                                   'mr_egger_regression',
#                                   'mr_weighted_median',
#                                   'mr_weighted_mode'))
# h_test = mr_heterogeneity(h_dat)
# p_test = mr_pleiotropy_test(h_dat)
# df <- df %>% bind_rows(
#   data.frame(
#     exposure = 'PRS1',
#     outcome  = 'Severe insulin-resistant T2D',
#     method   = mr_res$method,
#     beta     = mr_res$b,
#     se       = mr_res$se,
#     p        = mr_res$pval,
#     p_Qval   = c(h_test$Q_pval[2]),
#     p_pleio  = c(p_test$pval)
#   ))
# 
# h_dat <- harmonise_data(prs2,t2ds3_mr)
# h_dat <- h_dat %>% mutate(mr_keep=T)
# mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
#                                   'mr_egger_regression',
#                                   'mr_weighted_median',
#                                   'mr_weighted_mode'))
# 
# h_test = mr_heterogeneity(h_dat)
# p_test = mr_pleiotropy_test(h_dat)
# df <- df %>% bind_rows(
#   data.frame(
#     exposure = 'PRS2',
#     outcome  = 'Severe insulin-resistant T2D',
#     method   = mr_res$method,
#     beta     = mr_res$b,
#     se       = mr_res$se,
#     p        = mr_res$pval,
#     p_Qval   = c(h_test$Q_pval[2]),
#     p_pleio  = c(p_test$pval)
#   ))
# 
# h_dat <- harmonise_data(prs3,t2ds3_mr)
# h_dat <- h_dat %>% mutate(mr_keep=T)
# mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
#                                   'mr_egger_regression',
#                                   'mr_weighted_median',
#                                   'mr_weighted_mode'))
# 
# h_test = mr_heterogeneity(h_dat)
# p_test = mr_pleiotropy_test(h_dat)
# df <- df %>% bind_rows(
#   data.frame(
#     exposure = 'PRS3',
#     outcome  = 'Severe insulin-resistant T2D',
#     method   = mr_res$method,
#     beta     = mr_res$b,
#     se       = mr_res$se,
#     p        = mr_res$pval,
#     p_Qval   = c(h_test$Q_pval[2]),
#     p_pleio  = c(p_test$pval)
#   ))
# 
# h_dat <- harmonise_data(prs4,t2ds3_mr)
# h_dat <- h_dat %>% mutate(mr_keep=T)
# mr_res = mr_wald_ratio(b_exp = h_dat$beta.exposure,
#                        b_out = h_dat$beta.outcome,
#                        se_exp = h_dat$se.exposure,
#                        se_out = h_dat$se.outcome)
# 
# df <- df %>% bind_rows(
#   data.frame(
#     exposure = 'PRS4',
#     outcome  = 'Severe insulin-resistant T2D',
#     method   = 'Wald ratio',
#     beta     = mr_res$b,
#     se       = mr_res$se,
#     p        = mr_res$pval,
#     p_Qval   = NA,
#     p_pleio  = NA
#   ))
# 
# t2d_s4 <- fread('GCST90026415_buildGRCh37.tsv.gz')
# t2d_s4 <- t2d_s4 %>% filter(variant_ID %in% nafld$rsid)
# t2ds4_mr <- format_data(as.data.frame(t2d_s4),
#                         type = 'outcome',
#                         snp_col =  'variant_ID',
#                         beta_col = 'beta',
#                         se_col = 'standard_error',
#                         eaf_col = 'effect_allele_frequency',
#                         effect_allele_col = 'effect_allele',
#                         other_allele_col = 'other_allele',
#                         pval_col = 'p_value',
#                         chr_col = 'chromosome',
#                         pos_col = 'base_pair_location')
# h_dat <- harmonise_data(nafld_mr,t2ds4_mr)
# h_dat <- h_dat %>% mutate(mr_keep=T)
# mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
#                                   'mr_egger_regression',
#                                   'mr_weighted_median',
#                                   'mr_weighted_mode'))
# h_test = mr_heterogeneity(h_dat)
# p_test = mr_pleiotropy_test(h_dat)
# df <- df %>% bind_rows(data.frame(
#   exposure = 'Overall',
#   outcome  = 'Mild obesity-related T2D',
#   method   = mr_res$method,
#   beta     = mr_res$b,
#   se       = mr_res$se,
#   p        = mr_res$pval,
#   p_Qval   = c(h_test$Q_pval[2],NA,NA,NA),
#   p_pleio  = c(p_test$pval,NA,NA,NA)
# ))
# h_dat <- harmonise_data(prs1,t2ds4_mr)
# h_dat <- h_dat %>% mutate(mr_keep=T)
# mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
#                                   'mr_egger_regression',
#                                   'mr_weighted_median',
#                                   'mr_weighted_mode'))
# h_test = mr_heterogeneity(h_dat)
# p_test = mr_pleiotropy_test(h_dat)
# df <- df %>% bind_rows(
#   data.frame(
#     exposure = 'PRS1',
#     outcome  = 'Mild obesity-related T2D',
#     method   = mr_res$method,
#     beta     = mr_res$b,
#     se       = mr_res$se,
#     p        = mr_res$pval,
#     p_Qval   = c(h_test$Q_pval[2]),
#     p_pleio  = c(p_test$pval)
#   ))
# 
# h_dat <- harmonise_data(prs2,t2ds4_mr)
# h_dat <- h_dat %>% mutate(mr_keep=T)
# mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
#                                   'mr_egger_regression',
#                                   'mr_weighted_median',
#                                   'mr_weighted_mode'))
# 
# h_test = mr_heterogeneity(h_dat)
# p_test = mr_pleiotropy_test(h_dat)
# df <- df %>% bind_rows(
#   data.frame(
#     exposure = 'PRS2',
#     outcome  = 'Mild obesity-related T2D',
#     method   = mr_res$method,
#     beta     = mr_res$b,
#     se       = mr_res$se,
#     p        = mr_res$pval,
#     p_Qval   = c(h_test$Q_pval[2]),
#     p_pleio  = c(p_test$pval)
#   ))
# 
# h_dat <- harmonise_data(prs3,t2ds4_mr)
# h_dat <- h_dat %>% mutate(mr_keep=T)
# mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
#                                   'mr_egger_regression',
#                                   'mr_weighted_median',
#                                   'mr_weighted_mode'))
# 
# h_test = mr_heterogeneity(h_dat)
# p_test = mr_pleiotropy_test(h_dat)
# df <- df %>% bind_rows(
#   data.frame(
#     exposure = 'PRS3',
#     outcome  = 'Mild obesity-related T2D',
#     method   = mr_res$method,
#     beta     = mr_res$b,
#     se       = mr_res$se,
#     p        = mr_res$pval,
#     p_Qval   = c(h_test$Q_pval[2]),
#     p_pleio  = c(p_test$pval)
#   ))
# 
# h_dat <- harmonise_data(prs4,t2ds4_mr)
# h_dat <- h_dat %>% mutate(mr_keep=T)
# mr_res = mr_wald_ratio(b_exp = h_dat$beta.exposure,
#                        b_out = h_dat$beta.outcome,
#                        se_exp = h_dat$se.exposure,
#                        se_out = h_dat$se.outcome)
# 
# df <- df %>% bind_rows(
#   data.frame(
#     exposure = 'PRS4',
#     outcome  = 'Mild obesity-related T2D',
#     method   = 'Wald ratio',
#     beta     = mr_res$b,
#     se       = mr_res$se,
#     p        = mr_res$pval,
#     p_Qval   = NA,
#     p_pleio  = NA
#   ))
# 
# t2d_s5 <- fread('GCST90026416_buildGRCh37.tsv.gz')
# t2d_s5 <- t2d_s5 %>% filter(variant_ID %in% nafld$rsid)
# t2ds5_mr <- format_data(as.data.frame(t2d_s5),
#                         type = 'outcome',
#                         snp_col =  'variant_ID',
#                         beta_col = 'beta',
#                         se_col = 'standard_error',
#                         eaf_col = 'effect_allele_frequency',
#                         effect_allele_col = 'effect_allele',
#                         other_allele_col = 'other_allele',
#                         pval_col = 'p_value',
#                         chr_col = 'chromosome',
#                         pos_col = 'base_pair_location')
# h_dat <- harmonise_data(nafld_mr,t2ds5_mr)
# h_dat <- h_dat %>% mutate(mr_keep=T)
# mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
#                                   'mr_egger_regression',
#                                   'mr_weighted_median',
#                                   'mr_weighted_mode'))
# h_test = mr_heterogeneity(h_dat)
# p_test = mr_pleiotropy_test(h_dat)
# df <- df %>% bind_rows(data.frame(
#   exposure = 'Overall',
#   outcome  = 'Mild age-related T2D',
#   method   = mr_res$method,
#   beta     = mr_res$b,
#   se       = mr_res$se,
#   p        = mr_res$pval,
#   p_Qval   = c(h_test$Q_pval[2],NA,NA,NA),
#   p_pleio  = c(p_test$pval,NA,NA,NA)
# ))
# h_dat <- harmonise_data(prs1,t2ds5_mr)
# h_dat <- h_dat %>% mutate(mr_keep=T)
# mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
#                                   'mr_egger_regression',
#                                   'mr_weighted_median',
#                                   'mr_weighted_mode'))
# h_test = mr_heterogeneity(h_dat)
# p_test = mr_pleiotropy_test(h_dat)
# df <- df %>% bind_rows(
#   data.frame(
#     exposure = 'PRS1',
#     outcome  = 'Mild age-related T2D',
#     method   = mr_res$method,
#     beta     = mr_res$b,
#     se       = mr_res$se,
#     p        = mr_res$pval,
#     p_Qval   = c(h_test$Q_pval[2]),
#     p_pleio  = c(p_test$pval)
#   ))
# 
# h_dat <- harmonise_data(prs2,t2ds5_mr)
# h_dat <- h_dat %>% mutate(mr_keep=T)
# mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
#                                   'mr_egger_regression',
#                                   'mr_weighted_median',
#                                   'mr_weighted_mode'))
# 
# h_test = mr_heterogeneity(h_dat)
# p_test = mr_pleiotropy_test(h_dat)
# df <- df %>% bind_rows(
#   data.frame(
#     exposure = 'PRS2',
#     outcome  = 'Mild age-related T2D',
#     method   = mr_res$method,
#     beta     = mr_res$b,
#     se       = mr_res$se,
#     p        = mr_res$pval,
#     p_Qval   = c(h_test$Q_pval[2]),
#     p_pleio  = c(p_test$pval)
#   ))
# 
# h_dat <- harmonise_data(prs3,t2ds5_mr)
# h_dat <- h_dat %>% mutate(mr_keep=T)
# mr_res = mr(h_dat,method_list = c('mr_ivw_fe',
#                                   'mr_egger_regression',
#                                   'mr_weighted_median',
#                                   'mr_weighted_mode'))
# 
# h_test = mr_heterogeneity(h_dat)
# p_test = mr_pleiotropy_test(h_dat)
# df <- df %>% bind_rows(
#   data.frame(
#     exposure = 'PRS3',
#     outcome  = 'Mild age-related T2D',
#     method   = mr_res$method,
#     beta     = mr_res$b,
#     se       = mr_res$se,
#     p        = mr_res$pval,
#     p_Qval   = c(h_test$Q_pval[2]),
#     p_pleio  = c(p_test$pval)
#   ))
# 
# h_dat <- harmonise_data(prs4,t2ds5_mr)
# h_dat <- h_dat %>% mutate(mr_keep=T)
# mr_res = mr_wald_ratio(b_exp = h_dat$beta.exposure,
#                        b_out = h_dat$beta.outcome,
#                        se_exp = h_dat$se.exposure,
#                        se_out = h_dat$se.outcome)
# 
# df <- df %>% bind_rows(
#   data.frame(
#     exposure = 'PRS4',
#     outcome  = 'Mild age-related T2D',
#     method   = 'Wald ratio',
#     beta     = mr_res$b,
#     se       = mr_res$se,
#     p        = mr_res$pval,
#     p_Qval   = NA,
#     p_pleio  = NA
#   ))
# write.csv(df,'MR_T2D_subtypes.csv',row.names = F)
# 
# ####UKB-PPP
# rs1260326 <- read.csv('rs1260326-Protein-UKBPPP.csv')
# rs28601761 <- read.csv('rs28601761-Proteins-UKBPPP.csv')
# olink <- fread('F:\\E\\UkbData\\olink_protein_map_3k_v1.tsv')
# 
# rs1260326 <- rs1260326 %>% 
#   left_join(olink %>% select(UniProt,Assay,chr,gene_start,gene_end),
#             by = c('Target.1'='UniProt'))
# 
# rs1260326 <- rs1260326 %>% mutate(
#   anno = ifelse(chr.y==2 & (gene_start>26459497 & gene_start<28556649),
#                 'cis','trans')
# )
# ggplot(rs1260326,aes(beta,-log10(p)))+
#   geom_point()+
#   theme_bw()

####cumulative risk of T2D in PRS subgroups
# fit_cox1 <- coxph(Surv(followUp_days_T2D,T2D)~sex+age+smokingStatus+averageTotalIncome+
#                     aolIntakeFreq+no.daysModeAct2+education+loweringINS+
#                     prs1,ukb6)
# summary(fit_cox1)
# 

###pPRS and protein levels
olink <- fread('F:\\E\\UkbData\\UKB-PPP-analyzingData.csv')
ukb7 <- olink %>% left_join(
  ukb6 %>% select(eid,sex,age,smokingStatus,averageTotalIncome,
                  aolIntakeFreq2,no.daysModeAct2,education,BMI,
                  prs1,prs2,prs3,prs4,rs1260326), 
  by = c('eid.63726'='eid')
)
ukb7 <- ukb7 %>% mutate(
  prs234 = prs2+prs3+prs4
)
proteins <- names(olink)[-1]
df2 <- data.frame()
for(i in 1:length(proteins)){
  pro = get(proteins[i],ukb7)
  fit1 = lm(pro~sex+age+smokingStatus+averageTotalIncome+aolIntakeFreq2+
              no.daysModeAct2+education+BMI+scale(prs1),data = ukb7)
  fit2 = lm(pro~sex+age+smokingStatus+averageTotalIncome+aolIntakeFreq2+
              no.daysModeAct2+education+BMI+scale(prs234),data = ukb7)
  x1   = broom::tidy(fit1)
  x2   = broom::tidy(fit2)
  df = data.frame(
    protein = proteins[i],
    PRS     = c('PRS1','PRS2-4'),
    beta    = c(x1$estimate[12],
                x2$estimate[12]),
    se      = c(x1$std.error[12],
                x2$std.error[12]),
    p       = c(x1$p.value[12],
                x2$p.value[12])
  )
  df2 = rbind(df2,df)
}

olink2 <- fread('F:\\E\\UkbData\\olink_protein_map_3k_v1.tsv')

df3 <- df2 %>% left_join(olink2 %>% select(HGNC.symbol,chr,
                                          gene_start,gene_end,
                                          UniProt) %>%
                           distinct(HGNC.symbol,.keep_all = T),
                         by = c('protein'='HGNC.symbol'))
df3 <- df3 %>% filter(!chr == 'X',!is.na(chr)) %>%
  mutate(chr=as.numeric(chr))

write.csv(df3,'PRS-protein.csv',row.names = F)


proteins <- names(olink)[-1]
df2 <- data.frame()
for(i in 1:length(proteins)){
  pro = get(proteins[i],ukb7)
  fit1 = lm(pro~sex+age+smokingStatus+averageTotalIncome+aolIntakeFreq2+
              no.daysModeAct2+education+BMI+rs1260326,data = ukb7)
  x1   = broom::tidy(fit1)
  df = data.frame(
    protein = proteins[i],
    PRS     = c('rs1260326'),
    beta    = c(x1$estimate[12]),
    se      = c(x1$std.error[12]),
    p       = c(x1$p.value[12])
  )
  df2 = rbind(df2,df)
}
df3 <- df2 %>% left_join(olink2 %>% select(HGNC.symbol,chr,
                                           gene_start,gene_end,
                                           UniProt) %>%
                           distinct(HGNC.symbol,.keep_all = T),
                         by = c('protein'='HGNC.symbol'))
df3 <- df3 %>% filter(!chr == 'X',!is.na(chr)) %>%
  mutate(chr=as.numeric(chr))
write.csv(df3,'rs1260326-protein-UKBB.csv',row.names = F)


df <-read.csv('PRS-protein.csv')
# Convert chromosome to factor and ensure numeric conversion for positioning
df$chr <- as.factor(df$chr)
df$gene_start <- as.numeric(df$gene_start)
df$p <- as.numeric(df$p)

# Define genome-wide significance threshold
data <- df

# Bonferroni校正阈值
bonferroni_threshold <- 0.05 / nrow(data)
#8.896797e-06

# 显著性筛选 (Bonferroni校正)
prs1_significant <- data %>% filter(PRS == "PRS1" & p < bonferroni_threshold)
prs2_4_significant <- data %>% filter(PRS == "PRS2-4" & p < bonferroni_threshold)

# 交集蛋白
common_proteins <- intersect(prs1_significant$protein, prs2_4_significant$protein)

# 打印交集蛋白数量
cat("Number of common proteins significantly associated with both PRS1 and PRS2-4 (Bonferroni corrected):", length(common_proteins), "\n")
#140


# 火山图绘制函数
volcano_plot <- function(data, title) {
  ggplot(data, aes(x = beta, y = -log10(p))) +
    geom_point(aes(color = p < bonferroni_threshold), size = 2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    scale_color_manual(values = c("gray", "red")) +
    theme_minimal() +
    ggtitle(title) +
    xlab("Beta Coefficient") +
    ylab("-log10(p-value)") +
    geom_hline(yintercept = -log10(bonferroni_threshold), linetype = "dashed", color = "blue") +
    geom_text(aes(label = ifelse(p < bonferroni_threshold, as.character(protein), "")), 
              vjust = 1.5, hjust = 0.5, size = 2.5)
}

# PRS1火山图
prs1_volcano <- volcano_plot(prs1_significant, "PRS1-associated Proteins (Bonferroni)")

# PRS2-4火山图
prs2_4_volcano <- volcano_plot(prs2_4_significant, "PRS2-4-associated Proteins (Bonferroni)")

# 显示火山图
print(prs1_volcano)
print(prs2_4_volcano)



###曼哈顿图
genecard <- read.csv('GeneCards-glucoseMetabolism.csv')
df <- df %>% left_join(genecard %>% select(Uniprot.ID,Relevance.score,Description
                                           ),
                       by = c('UniProt'='Uniprot.ID'))

df2 <- df %>% filter(p < 8.89e-6) %>% distinct(protein,.keep_all = T)
write.csv(df2,'Significant-protein-genecard.csv',row.names = F)
# Compute -log10(P-value)
# df$logP <- -log10(df$p)
# 
# # Sort data by chromosome and gene start position
# df <- df %>% arrange(chr, gene_start)
# 
# # Get the **maximum gene start position** per chromosome
# chr_offsets <- df %>%
#   group_by(chr) %>%
#   summarize(chr_max = max(gene_start, na.rm=TRUE)) %>%
#   ungroup()
# 
# # Compute **cumulative genome position**
# chr_offsets <- chr_offsets %>%
#   mutate(chr_start = cumsum(lag(chr_max, default = 0) + 1e6))  # Add 1Mb spacing
# 
# # Merge back to the dataset to adjust positions
# df <- df %>%
#   left_join(chr_offsets, by = "chr") %>%
#   mutate(Pos_Cumulative = gene_start + chr_start)
# 
# # Compute chromosome midpoints for proper x-axis spacing
# chr_ticks <- df %>%
#   group_by(chr) %>%
#   summarize(center = mean(Pos_Cumulative))
# 
# # Assign color logic:
# df$fill_group <- case_when(
#   df$p <= bonferroni_threshold & df$PRS == "PRS1" ~ "PRS1 Significant",    # Blue for PRS1 significant
#   df$p <= bonferroni_threshold & df$PRS == "PRS2-4" ~ "PRS234 Significant",  # Red for PRS234 significant
#   TRUE ~ "Not Significant"  # Gray for non-significant
# )
# # '#ED0000','#00468b','#0099b4',"#42b540"
# # # Define color scheme:
# # fill_palette <- c("Not Significant" = "gray70",
# #                    "PRS1 Significant" = "#42b540",
# #                    "PRS234 Significant" = "#ED0000")
# df$shape_group <- case_when(
#   df$p <= bonferroni_threshold & df$beta > 0 ~ "Up",   # Up-triangle for β > 0
#   df$p <= bonferroni_threshold & df$beta < 0 ~ "Down", # Down-triangle for β < 0
#   TRUE ~ "Circle"  # Keep circles for non-significant points
# )
# # Define shape scheme:
# shape_palette <- c("Circle" = 21, "Up" = 24, "Down" = 25)
# # Create Manhattan Plot with **correct chromosome spacing**
# a <- ggplot(df, aes(x = Pos_Cumulative, y = logP, fill = fill_group, 
#                shape = shape_group,color=fill_group)) +
#   geom_point(alpha = 0.8, size = 1, stroke = 0) +  # Shape 16 = Circle for all points
#   scale_fill_manual(values = c('gray70','#56B4E9','#ec347c')) +  # Assign colors
#   scale_color_manual(values = c('gray70','#56B4E9','#ec347c')) +
#   scale_shape_manual(values = c(21,25,24))+
#   geom_hline(yintercept = -log10(bonferroni_threshold), 
#              linetype = "dashed", color = "black") +  # Significance threshold line
#   scale_x_continuous(breaks = chr_ticks$center, labels = chr_ticks$chr,
#                      expand = expansion(add=1e7)) +  # Properly spaced chromosome labels
#   theme_minimal() +
#   ggrepel::geom_text_repel(aes(label = protein),
#                            data = df%>%filter(Relevance.score>=25,
#                                               p<=bonferroni_threshold),
#                            size=3,show.legend = F)+
#   labs(x = "Chromosome",
#        y = "-log10(P-value)",
#        color = "Association Type") +
#   guides(fill = "none",
#          shape= 'none') +
#   theme(legend.position = "top",
#         axis.text.x = element_text(angle = 90, hjust = 1))

#####环形图
# ==== 加载必要包 ====
library(circlize)
library(scales)

# ==== 读取数据 ====
df <- read_csv("circularPlot.data.csv", show_col_types = FALSE)

# 添加 logP 和 significance
df <- df %>%
  mutate(
    logP = -log10(p),
    significant = p < (0.05 / 5620),
    pos = (gene_start + gene_end) / 2
  )

# 最外层 label 信息筛选
label_df <- df %>%
  filter(significant, Relevance.score > 25) %>%
  group_by(protein) %>%
  slice(1) %>%  # 去重
  mutate(label_group = case_when(
    n_distinct(PRS) > 1 ~ "Shared",
    PRS == "PRS1" ~ "PRS1",
    PRS == "PRS2-4" ~ "PRS2-4"
  ))

# 设置环图的扇区
chromosomes <- unique(df$chr)
par(mar = rep(0.000001, 4),xpd=NA)  # 极限压缩图像边距
circos.clear()
circos.par(
  start.degree = 90,
  track.height = 0.17,               # 稍大点让每层更饱满
  track.margin = c(0.003, 0.003),    # 层间间隔更紧凑
  gap.degree = 3,
  cell.padding = c(0, 0, 0, 0),      # 取消每层内部留白
  circle.margin = 0.00001              # 环图与画布边界压紧
)
circos.initialize(factors = df$chr, x = df$pos)

# ==== 最外层：标签 ====
# === 构造标签数据（位置 + 标签 + 颜色） ===
# 提前提取 Shared 蛋白名列表
shared_proteins <- df %>%
  filter(significant, Relevance.score > 25) %>%
  group_by(protein) %>%
  summarise(n = n_distinct(PRS), .groups = "drop") %>%
  filter(n > 1) %>%
  pull(protein)

# 构造标签数据（每个蛋白只保留一行）
label_df_unique <- df %>%
  filter(significant, Relevance.score > 25) %>%
  group_by(protein) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    label_group = case_when(
      protein %in% shared_proteins ~ "Shared",
      PRS == "PRS1" ~ "PRS1",
      PRS == "PRS2-4" ~ "PRS2-4"
    ),
    label_color = case_when(
      label_group == "Shared" ~ "black",
      label_group == "PRS1" ~ "#56B4E9",     # 蓝色
      label_group == "PRS2-4" ~ "#ec347c"    # 粉红色
    )
  )


# === 绘图代码：最外层 ===
label_bed <- label_df_unique %>%
  select(chr, gene_start, gene_end, protein, label_color) %>%
  mutate(chr = as.character(chr)) %>%    # 如果你用的是字符型扇区名
  drop_na()
circos.genomicLabels(
  label_bed,
  labels.column = 4,
  side = "outside",
  line_col = label_bed$label_color,
  col = label_bed$label_color,
  cex = 0.8,
  connection_height = mm_h(3),
  facing = "clockwise",
  niceFacing = TRUE
)



# 第二层：染色体分区
circos.trackPlotRegion(
  ylim = c(0, 1), 
  #bg.col = rep(c("gray90", NA), length.out = 22),
  bg.border = NA,
  track.height = 0.08,
  panel.fun = function(x, y) {
    chr = CELL_META$sector.index
    xlim = CELL_META$xlim
    circos.text(mean(xlim), 0.5, labels = chr, 
                cex = 0.8, facing = "inside", niceFacing = TRUE)
  }
)

# 第三层：PRS1
circos.trackPlotRegion(
  ylim = c(0, max(df$logP, na.rm = TRUE)), bg.border = NA, 
  track.height = 0.2,
  bg.col = rep(c("gray90", NA), length.out = 22),
  panel.fun = function(region, value, ...) {
    chr_current <- CELL_META$sector.index
    subdf <- df %>% filter(chr == chr_current, PRS == "PRS1")
    
    rel_score <- replace_na(subdf$Relevance.score, 0)
    rel_score <- pmin(rel_score, 50)
    cexs <- 0.4
    col <- ifelse(subdf$significant, "#56B4E9", "gray")
    
    circos.points(subdf$pos, subdf$logP, col = col, cex = cexs, pch = 16)
  }
)

# 第四层：PRS2-4
circos.trackPlotRegion(
  ylim = c(0, max(df$logP, na.rm = TRUE)), bg.border = NA, 
  track.height = 0.2,
  bg.col = rep(c("gray90", NA), length.out = 22),
  panel.fun = function(region, value, ...) {
    chr_current <- CELL_META$sector.index
    subdf <- df %>% filter(chr == chr_current, PRS == "PRS2-4")
    
    rel_score <- replace_na(subdf$Relevance.score, 0)
    rel_score <- pmin(rel_score, 50)
    cexs <- 0.4
    col <- ifelse(subdf$significant, "#ec347c", "gray")
    
    circos.points(subdf$pos, subdf$logP, col = col, cex = cexs, pch = 16)
  }
)
circos.yaxis(side = "left", at = c(0, 50,100),
             sector.index = "1", track.index = 4, labels.cex = 0.4)
circos.yaxis(side = "left", at = c(0, 50,100), sector.index = "1", 
             track.index = 5, labels.cex = 0.4)

# 图例
legend("topleft", legend = c("PRS1", "PRS2-4", "Shared"),
       col = c("#56B4E9", "#ec347c", "black"), pch = 16, cex = 0.8, bty = "n")

###upset
# library(UpSetR)
# df <- read_csv("PRS-protein.csv")
# 
# # 设置显著性阈值
# threshold <- 0.05 / 5620
# 
# # 过滤显著结果
# df_sig <- df %>% filter(p < threshold)
# 
# # 创建 wide 格式的 presence-absence 矩阵
# df_matrix <- df_sig %>%
#   mutate(sig = 1) %>%
#   select(protein, PRS, sig) %>%
#   pivot_wider(names_from = PRS, values_from = sig, values_fill = 0) %>%
#   column_to_rownames("protein")
# 
# # 绘制 upset 图
# upset(df_matrix, sets = c("PRS1", "PRS2-4"),
#       order.by = "freq", keep.order = TRUE,
#       sets.bar.color = "steelblue",
#       mainbar.y.label = "Intersection size",
#       sets.x.label = "Significant proteins per PRS")

# 拆分两个 PRS 的结果
prs1 <- df %>% filter(PRS == "PRS1") %>%
  select(protein, beta_PRS1 = beta, p_PRS1 = p)

prs2 <- df %>% filter(PRS == "PRS2-4") %>%
  select(protein, beta_PRS2_4 = beta, p_PRS2_4 = p)

# 合并
merged <- inner_join(prs1, prs2, by = "protein")

# 标记类别
merged <- merged %>%
  mutate(
    opposite_significant = (p_PRS1 < threshold) & (p_PRS2_4 < threshold) & (beta_PRS1 * beta_PRS2_4 < 0),
    consistent_significant = (p_PRS1 < threshold) & (p_PRS2_4 < threshold) & (beta_PRS1 * beta_PRS2_4 > 0),
    PRS1_only = (p_PRS1 < threshold) & (p_PRS2_4 >= threshold),
    PRS2_4_only = (p_PRS2_4 < threshold) & (p_PRS1 >= threshold)
  )
write.csv(merged,'PRS1_PRS2_4_PROTEINS.csv',row.names = F)
# 统计数量
summary_df <- tibble(
  Category = c("Opposite & both significant", "PRS1 only significant", "PRS2-4 only significant"),
  Count = c(
    sum(merged$opposite_significant),
    sum(merged$PRS1_only),
    sum(merged$PRS2_4_only)
  )
)

# 绘制条形图
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = c("steelblue", "forestgreen", "orange")) +
  theme_minimal() +
  labs(
    title = "Differential PRS-protein associations",
    x = NULL,
    y = "Number of significant proteins"
  ) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))



###enrichment analysis
library(simplifyEnrichment)
data <- read.csv('C:\\Users\\Liu Zhenqiu\\Downloads\\gProfiler_hsapiens_2025-6-24 15-05-19__intersections.csv')
# --- 1. 筛选只保留GO:BP条目 ---
data_gobp <- data %>%
  filter(source == "GO:BP")

# --- 2. 准备GO term列表 ---
go_ids <- data_gobp$term_id

# --- 3. 计算GO terms相似性矩阵 ---
# 直接用simplifyEnrichment封装好的接口
mat <- GO_similarity(go_ids)

# --- 4. 聚类 + 去冗余 ---
# 自动cut分组，保留每组代表
clusters <- simplifyGO(mat)

# 先把聚类信息和富集结果表连接起来
data_with_cluster <- data_gobp %>%
  inner_join(clusters, by = c("term_id" = "id"))

# 在每个cluster内，选p值最小的GO term作为代表
representative_terms <- data_with_cluster %>%
  group_by(cluster) %>%
  slice_min(order_by = adjusted_p_value, n = 1, with_ties = FALSE) %>%
  ungroup()

representative_terms <- representative_terms %>%
  mutate(log10_p = -log10(adjusted_p_value))

ggplot(representative_terms, aes(x = log10_p, y = reorder(term_name, log10_p))) +
  geom_point(aes(size = intersection_size, color = log10_p)) +
  scale_color_gradient(low = "#a8ddb5", high = "#084594") +
  scale_size_continuous(range = c(2, 6)) +
  theme_minimal(base_size = 12) +
  labs(
    #title = "Simplified GO:BP Terms",
    x = "-log10(p-value)",
    y = NULL,
    size = "Gene Count",
    color = "-log10(p-value)"
  ) +
  theme(
    axis.text.y = element_text(size = 11,color='black'),
    legend.position = "right"
  ) +
  guides(color = "none")  # 如果你不想显示色条legend



####mediation analysis: PRS --> protein --> T2D
ukb7 <- olink %>% dplyr::select(eid.63726,opposite_proteins$protein) %>%
  left_join(
    ukb6 %>% dplyr::select(eid,sex,age,smokingStatus,averageTotalIncome,
                           aolIntakeFreq2,no.daysModeAct2,education,BMI,loweringINS,
                           prs1,prs2,prs3,prs4,T2D,followUp_days_T2D), 
    by = c('eid.63726'='eid')
  ) %>% mutate(prs234 = prs2+prs3+prs4)

pros <- opposite_proteins$protein
df2 <- data.frame()
for(i in 1:61){
  pro = get(pros[i],ukb7)
  fit_cox <- coxph(Surv(followUp_days_T2D,T2D)~sex+age+averageTotalIncome+aolIntakeFreq2+
                     BMI+education+smokingStatus+no.daysModeAct2+loweringINS+
                     scale(pro),
                   data = ukb7)
  x = broom::tidy(fit_cox)
  df <- data.frame(protein = pros[i],
                   outcome = 'T2D',
                   beta    = x$estimate[12],
                   se      = x$std.error[12],
                   p       = x$p.value[12])
  df2 <- rbind(df2,df)
}
df2 <- df2 %>% mutate(
  hr = exp(beta),
  lwr= exp(beta-1.96*se),
  upr= exp(beta+1.96*se),
  sig= ifelse(p<0.05/61,'P < 0.05/61','P >= 0.05/61')
)
write.csv(df2,'opposite_protein_T2D.csv',row.names = F)


df <- read.csv('opposite_protein_T2D.csv')
df <- df %>%
  mutate(sig_flag = sig == "P < 0.05/61")

# 按 HR 排序蛋白并设置为 factor
df <- df %>%
  arrange(hr) %>%
  mutate(protein = factor(protein, levels = protein))

# 绘图
ggplot(df, aes(x = hr, y = protein)) +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0, color = "black") +
  geom_point(aes(shape = sig_flag), size = 2, fill = "#e64b35") +
  scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 1)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  labs(x = "Hazard Ratio (HR)", y = NULL, shape = "Significance") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 9),
    axis.text = element_text(size=9,color = 'black'),
    legend.position = "none"
  )

df2 <- df2 %>% rename(
  beta.protein.t2d = beta,
  se.protein.t2d   = se,
  p.protein.t2d    = p
)

df3 <- read.csv('PRS-protein.csv')
oppo_pro_data <- df3 %>% filter(protein %in% df2$protein)
oppo_pro_data <- oppo_pro_data %>% dplyr::select(protein,beta,PRS)%>%
  pivot_wider(names_from = PRS,values_from = c(beta))
oppo_pro_data <- oppo_pro_data %>% rename(
  beta.PRS1.protein = PRS1,
  beta.PRS2_4.protein=`PRS2-4`
)

df2 <- df2 %>% left_join(oppo_pro_data,by='protein')
df2 <- df2 %>% mutate(
  PRS1_med = ifelse(beta.protein.t2d*beta.PRS1.protein<0&sig == 'P < 0.05/61',
                    'yes','no'),
  PRS234_med = ifelse(beta.protein.t2d*beta.PRS2_4.protein>0&sig == 'P < 0.05/61',
                      'yes','no')
)
ukb7 <- ukb7 %>% mutate(
  prs1_scale = as.numeric(scale(prs1))
)
library(mediation)
df4 <- data.frame()
for(i in 1:61){
  pro = get(pros[i],ukb7)
  if(pros[i]%in%df2$protein[df2$PRS1_med=='yes']){
    model_m <- lm(pro~sex+age+averageTotalIncome+aolIntakeFreq2+
                    BMI+education+smokingStatus+no.daysModeAct2+prs1_scale,
                  data = ukb7)
    model_y <- glm(T2D~sex+age+averageTotalIncome+aolIntakeFreq2+
                     BMI+education+smokingStatus+no.daysModeAct2+prs1_scale+pro, 
                   data = ukb7, family = binomial())
    med <- tryCatch({
      mediate(
        model.m = model_m,
        model.y = model_y,
        treat = "prs1_scale",
        mediator = 'pro',
        boot = TRUE,
        sims = 500
      )
    }, error = function(e) {
      message("❌ 中介分析失败：", conditionMessage(e))
      return(NULL)
    }) 
    df <- data.frame(X = 'PRS1',
                     M = pros[i],
                     Y = 'T2D',
                     ACME = med$d0,
                     ADE = med$z0,
                     TE = med$tau.coef,
                     PM = med$n0,
                     p_ACME = med$d0.p,
                     p_ADE = med$z0.p) 
  } else {
    df <- data.frame(X = 'PRS1',
                     M = pros[i],
                     Y = 'T2D',
                     ACME = NA,
                     ADE = NA,
                     TE = NA,
                     PM = NA,
                     p_ACME = NA,
                     p_ADE = NA) 
  }
  df4 <- rbind(df4,df)
  print(paste(i,'protein finished!!!'))
}
write.csv(df4,'PRS1_opposite_protein_mediation.csv',row.names = F)

df4 <- data.frame()
for(i in 1:61){
  pro = get(pros[i],ukb7)
  if(pros[i]%in%df2$protein[df2$PRS234_med=='yes']){
    model_m <- lm(pro~sex+age+averageTotalIncome+aolIntakeFreq2+
                    BMI+education+smokingStatus+no.daysModeAct2+ prs234_scale,
                  data = ukb7)
    model_y <- glm(T2D~sex+age+averageTotalIncome+aolIntakeFreq2+
                     BMI+education+smokingStatus+no.daysModeAct2+ prs234_scale+pro, 
                   data = ukb7, family = binomial())
    med <- tryCatch({
      mediate(
        model.m = model_m,
        model.y = model_y,
        treat = "prs234_scale",
        mediator = 'pro',
        boot = TRUE,
        sims = 500
      )
    }, error = function(e) {
      message("❌ 中介分析失败：", conditionMessage(e))
      return(NULL)
    }) 
    df <- data.frame(X = 'PRS234',
                     M = pros[i],
                     Y = 'T2D',
                     ACME = med$d0,
                     ADE = med$z0,
                     TE = med$tau.coef,
                     PM = med$n0,
                     p_ACME = med$d0.p,
                     p_ADE = med$z0.p) 
  } else {
    df <- data.frame(X = 'PRS234',
                     M = pros[i],
                     Y = 'T2D',
                     ACME = NA,
                     ADE = NA,
                     TE = NA,
                     PM = NA,
                     p_ACME = NA,
                     p_ADE = NA) 
  }
  df4 <- rbind(df4,df)
  print(paste(i,'protein finished!!!'))
}
write.csv(df4,'PRS234_opposite_protein_mediation.csv',row.names = F)

df_hr <- read_csv("opposite_protein_T2D.csv")
df_prs1 <- read_csv("PRS1_opposite_protein_mediation.csv")
df_prs234 <- read_csv("PRS234_opposite_protein_mediation.csv")

# 按照hr排序蛋白，并设为因子水平
ordered_proteins <- df_hr %>%
  arrange(hr) %>%
  pull(protein)

# 增加统一列名方便拼接
df_prs1 <- df_prs1 %>%
  mutate(protein = factor(M, levels = ordered_proteins),
         PRS = "PRS1") %>%
  mutate(PM = ifelse(PM > 0.5 | p_ACME >  0.05/61 | PM < 0, NA, PM))

df_prs234 <- df_prs234 %>%
  mutate(protein = factor(M, levels = ordered_proteins),
         PRS = "PRS234") %>%
  mutate(PM = ifelse(PM > 1 | p_ACME >  0.05/61 | PM < 0, NA, PM))

# 合并数据并转为热图矩阵格式
df_combined <- bind_rows(df_prs1, df_prs234) %>%
  mutate(protein = factor(protein, levels = ordered_proteins),
         PRS = factor(PRS, levels = c("PRS1", "PRS234")))

ggplot(df_combined, aes(x = PRS, y = protein, fill = PM)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "#e64b35", na.value = "grey90", name = "PM") +
  theme_void() +
  scale_x_discrete(expand = expansion(add=0))+
  theme(axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10),
        axis.title = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))




library(ggvenn)
# Venn图展示PRS1和PRS2-4交集蛋白
# 提取蛋白名称
set1 <- unique(prs1_significant$protein)
set2 <- unique(prs2_4_significant$protein)

# 创建数据列表（确保是向量）
venn_data <- list(
  `FBG-lowering PRS` = set1,
  `FBG-increasing PRS` = set2
)

# 绘制改进后的 Venn 图
venn_plot <- ggvenn(
  venn_data,
  fill_color = c('#56B4E9','#ec347c'),
  stroke_size = 0.5,
  text_size = 5,
  set_name_size = 6
) + 
  #ggtitle("Common Proteins between PRS1 and PRS2-4 (Bonferroni)") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# 找到交集蛋白
common_proteins <- intersect(prs1_significant$protein, prs2_4_significant$protein)

# 提取共同蛋白数据
common_data <- data %>% filter(protein %in% common_proteins)

# 转换为宽格式，便于比较
common_data_wide <- common_data %>%
  dplyr::select(protein, PRS, beta) %>%
  spread(key = PRS, value = beta)

# 计算效应值差异，用于颜色区分
common_data_wide <- common_data_wide %>%
  mutate(effect_diff = PRS1 * `PRS2-4`,
         regulation = case_when(
           PRS1 > 0 & `PRS2-4` < 0 ~ "Opposite",
           PRS1 < 0 & `PRS2-4` > 0 ~ "Opposite",
           TRUE ~ "Same"
         ))

# 绘制散点图
library(ggrepel)

# 仅筛选效应差异较大的 Opposite 蛋白进行标注
label_data <- common_data_wide %>%
  filter(regulation == "Opposite" & abs(PRS1 - `PRS2-4`) > 0.03)

# 绘制优化后的散点图
scatterPlot <- ggplot(common_data_wide, aes(x = PRS1, y = `PRS2-4`, color = regulation)) +
  # 显示散点
  geom_point(aes(size = abs(PRS1 - `PRS2-4`))) +
  
  # 仅标注效应差异较大的 Opposite 蛋白
  geom_text_repel(
    data = label_data,
    aes(label = protein),
    size = 3, color = "black", 
    box.padding = 0.5, max.overlaps = 50
  ) +
  
  # 辅助线
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  
  # 颜色设置：Opposite 蓝色，Same 浅灰色
  scale_color_manual(values = c("Opposite" = "#e64b35", "Same" = "#91d1c2")) +
  scale_size_continuous(name='Effect size\ndifference') +
  
  # 主题和标题
  labs(
    x = "Effect size (FBG-lowring PRS)",
    y = "Effect size (FBG-increasing PRS)",
    color = "Effect Direction"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

opposite_proteins <- common_data_wide %>% filter(regulation == 'Opposite')
write.csv(opposite_proteins,'opposite_proteins.csv',row.names = F)

# library(clusterProfiler)
# library(org.Hs.eg.db)  # Human gene annotation
# library(enrichplot)  # For visualization
# 
# entrez_ids <- bitr(opposite_proteins$protein, fromType = "SYMBOL", 
#                    toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# 
# # GO生物过程富集分析
# # GO富集分析（BP、MF、CC）
# go_bp <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
#                   pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
# go_mf <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "MF",
#                   pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
# go_cc <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "CC",
#                   pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
# 
# 
# # 转换为数据框并添加分类标签
# bp_df <- as.data.frame(go_bp) %>% mutate(Ontology = "BP")
# mf_df <- as.data.frame(go_mf) %>% mutate(Ontology = "MF")
# cc_df <- as.data.frame(go_cc) %>% mutate(Ontology = "CC")
# 
# # 合并三个GO分析结果
# go_df <- bind_rows(bp_df, mf_df, cc_df)
# 
# # 按 p 值排序，选择 top 10
# go_df <- go_df %>%
#   group_by(Ontology) %>%
#   top_n(-10, p.adjust) %>%
#   ungroup()
# 
# # 修改条目名称，便于展示
# go_df$Description <- factor(go_df$Description, levels = rev(unique(go_df$Description)))
# 
# 
# # 绘制三分类GO富集条形图
# ggplot(go_df, aes(x = -log10(p.adjust), y = Description, fill = Ontology)) +
#   geom_bar(stat = "identity", position = "identity",
#            alpha=.7) +
#   
#   # 添加标签到条柱内部，并从x=0开始
#   geom_text(aes(x=0,label = Description), 
#             color = "black", 
#             hjust = 0,    # 从x=0开始对齐
#             size = 3, 
#             nudge_x = 0.01) +  # 微调防止贴边
#   theme_minimal() +
#   # 调整颜色
#   scale_fill_manual(values = c("BP" = "#1f77b4", "MF" = "#ff7f0e", "CC" = "#2ca02c")) +
#   
#   # 去除图例和纵轴标签
#   theme(legend.position = "none") +
#   theme(axis.text.y = element_blank()) +  # 去除纵轴标签
#   
#   # 调整X轴起点
#   scale_x_continuous(expand = c(0, 0)) +
#   
#   # 标题和轴标签
#   labs(
#     x = "-log10(Adjusted P-value)",
#     y = ""
#   ) +
#   
#   # 分面展示三类GO
#   facet_wrap(~ Ontology, scales = "free_y", ncol = 1) +
#   
#   # 主题设置
#   
#   theme(
#     plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
#     axis.title.x = element_text(size = 12),
#     axis.title.y = element_text(size = 12),
#     axis.text = element_text(size = 12),
#     strip.text = element_text(size = 13, face = "bold"),
#     panel.grid.major.y = element_blank()  # 去除水平网格线
#   )
# 
# # KEGG通路富集分析
# kegg_enrich <- enrichKEGG(
#   gene = entrez_ids$ENTREZID,
#   organism = "hsa",
#   pAdjustMethod = "BH",
#   pvalueCutoff = 0.05,
#   qvalueCutoff = 0.05
# )
# 
# # 将富集结果转换为数据框
# kegg_df <- as.data.frame(kegg_enrich)
# 
# # 筛选显著性通路
# kegg_df <- kegg_df %>%
#   top_n(-10, p.adjust)
# 
# # 修改条目名称，便于展示
# kegg_df$Description <- factor(kegg_df$Description, levels = rev(kegg_df$Description))
# 
# dotplot(kegg_enrich, showCategory = 20, title = "KEGG Pathway Enrichment Analysis")




####MR analysis of opposite proteins on T2D
t2d <- fread('Mahajan.NG2018b.T2D.European.txt.gz')
pqtl <- read.csv('F:\\E\\gwas of protein(Iceland+fenland)\\pqtl.databases2.csv')
opposite_proteins <- common_data_wide$protein[common_data_wide$regulation=='Opposite']
sum(opposite_proteins %in% pqtl$Protein)

df <- data.frame()
for(i in 1:61){
  dat = pqtl %>% filter(Protein == opposite_proteins[i])
  dat = dat %>% filter(pos.exposure < start_position+1024*1024,
                       pos.exposure > start_position-1024*1024)
  if(nrow(dat)==0){
   print('no eligible pQTL!')
    next
  }
  dat_mr = format_data(as.data.frame(dat),
                       type = 'exposure',
                       snp_col = 'SNP',
                       beta_col = 'beta.exposure',
                       se_col = 'se.exposure',
                       eaf_col = 'eaf.exposure',
                       effect_allele_col = 'effect_allele.exposure',
                       other_allele_col = 'other_allele.exposure',
                       pval_col = 'pval.exposure',
                       chr_col = 'chr.exposure',
                       pos_col = 'pos.exposure')
  t2d2 = t2d %>% filter(name %in% dat$SNP)
  t2d_mr <- format_data(as.data.frame(t2d2),
                        type = 'outcome',
                        snp_col =  'name',
                        beta_col = 'Beta',
                        se_col = 'SE',
                        eaf_col = 'EAF',
                        effect_allele_col = 'EA',
                        other_allele_col = 'NEA',
                        pval_col = 'Pvalue',
                        chr_col = 'Chr',
                        pos_col = 'Pos')
  try1 <- try({hdat <- harmonise_data(dat_mr,t2d_mr)},
              silent = T)
  if(inherits(try1,'try-error')|nrow(hdat)==0){
    print(paste(i,'has been skipped due to no eligible SNP!'))
    next
  } else if(length(hdat$SNP)==1){
    mr_res = mr_wald_ratio(b_exp=hdat$beta.exposure,
                           b_out=hdat$beta.outcome,
                           se_exp=hdat$se.exposure,
                           se_out=hdat$se.outcome)
    df <- df %>% bind_rows(
      data.frame(exposure = opposite_proteins[i],
                 outcome  = 'T2D',
                 method = 'Wald ratio',
                 nsnp = 1,
                 beta = mr_res$b,
                 se = mr_res$se,
                 pval = mr_res$pval,
                 p_Qval   = NA,
                 p_pleio  = NA
      ))
  } else {
    hdat$mr_keep=T
    mr_res <- mr(hdat,method_list = c('mr_ivw_fe',
                                     'mr_egger_regression',
                                     'mr_weighted_median',
                                     'mr_weighted_mode'))
    h_test = mr_heterogeneity(hdat)
    p_test = mr_pleiotropy_test(hdat)
    df <- df %>% bind_rows(data.frame(
      exposure = opposite_proteins[i],
      outcome  = 'T2D',
      method   = mr_res$method,
      nsnp     = mr_res$nsnp,
      beta     = mr_res$b,
      se       = mr_res$se,
      p        = mr_res$pval,
      p_Qval   = c(h_test$Q_pval[2],NA,NA,NA),
      p_pleio  = c(p_test$pval,NA,NA,NA)
    ))
  }
  print(paste(i,'was finished :)'))
}
write.csv(df,'opposite_protein_MR_t2d.csv',row.names = F)



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
loc <- link_LD(loc, token = Sys.getenv("LDLINK_TOKEN"))
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
loc <- link_LD(loc, token = Sys.getenv("LDLINK_TOKEN"))
summary(loc)
locus_plot(loc, labels = c("index"),
           label_x = c(4, -5))


chr2 <- pdff %>% dplyr::filter(chr==8,pos > 126500031-1024*1024,
                               pos<126500031+1024*1024)
chr2 <- chr2 %>% dplyr::select(chr,pos,snp,beta,se,pval) %>% as.data.frame() 
loc <- locus(data = chr2, index_snp =  'rs28601761', flank = 1e5,
             LD = 'r2',
             ens_db = "EnsDb.Hsapiens.v75")
loc <- link_LD(loc, token = Sys.getenv("LDLINK_TOKEN"))
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
loc <- link_LD(loc, token = Sys.getenv("LDLINK_TOKEN"))
summary(loc)
locus_plot(loc, labels = c("index"),
           label_x = c(4, -5))



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

# P50053 KHK
# Q99674 CGREF1
# Q6UW56 ATRAID
# Q9H477 RBKS

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
  # 绘制误差条，调整宽度，确保分组一致
  geom_errorbarh(aes(xmin = beta - se, xmax = beta + se), height = 0.05, 
                 size = 0.8, position = position_dodge(width = 0.3)) +
  
  # 绘制点，使用相同的分组和错位调整
  geom_point(size = 3, position = position_dodge(width = 0.3)) +
  
  # 添加p值注释，稍微偏移，防止重叠
  geom_text(aes(label = paste0("p=", formatC(p, format = "e", digits = 2))), 
            hjust = 0.5, vjust = -1.2, size = 4, 
            position = position_dodge(width = 0.3),show.legend = F) +
  
  # 颜色设置
  scale_color_brewer(palette = "Set1") +
  
  # 标题和轴标签
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
