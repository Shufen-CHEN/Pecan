library(vcfR)
library(tidyverse)
setwd("~/Desktop/Pecan-update/Missing rate")
d_all <- read.csv("DPec23-8107_MADC_rename_rm10plusHaps-markerlevel.csv",check.names = F)
d_sum <- aggregate(d_all[,-c(1:3)],by=d_all[2],sum)
rownames(d_sum)<- d_sum[,1]
d_sum <- d_sum[,-1]
f1_ID <- read.table("sample_ids_F1.txt")
dive_ID <- read.table("sample_ids_diversity_rmF1parents.txt", sep = "\t", header = FALSE)

d_sum_f1 <- d_sum[,which(colnames(d_sum) %in% f1_ID$V1)]
d_sum_dive <- d_sum[,which(colnames(d_sum) %in% dive_ID$V1)]

##sample-level missing rate
missrate <- c()
for (i in 1:ncol(d_sum)) {
  missrate[i] <-length(which(d_sum[,i] < 10))/nrow(d_sum)
}
boxplot(missrate)

mr <- data.frame("ID"=colnames(d_sum), "missrate"= missrate)
mr_f1 <- mr[grep(pattern = "2016|2017|Lakota|87MX3_2.11",mr$ID),]
mr_dive <-  mr[-grep(pattern = "2016|2017|Lakota|87MX3_2.11",mr$ID),]

mr_all <- rbind.data.frame(mr,mr_f1,mr_dive)
mr_all$Population <- rep(c("Entire population", "F1 population","Diverse population"),c(376,190,186))
mr_all$Population <- factor(mr_all$Population,levels = c("Entire population", "F1 population","Diverse population"))


#write.csv(mr_all,"missingrate for all samples_update.csv",row.names = F)

give.n<- function(y) {
  return( 
    data.frame(
      y =1.1*max(y),  #may need to modify this depending on your data
      label = paste('n = ', length(y), '\n',
                    'mean =', round(mean(y), 1), '\n'
      )
    )
  )
}


png("line_missrate_pecan.png",height = 2000,width = 3000,res = 300)
ggplot(mr_all,aes(x=mr_all$Population,y=mr_all$missrate*100,
                                     fill=mr_all$Population))+geom_boxplot(width=0.5)+ 
  theme_bw()+ ylim(c(0,120))+
  #ggtitle("Sample-based missing rate")+
  theme(axis.text.x=element_text(angle=0,hjust=0.5),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(face = "bold"), axis.title.y = element_text(face = "bold"),
        legend.position = "none",text = element_text(size = 20))+
  xlab("Population")+ ylab("Missing rate (%)")+
  #facet_wrap(~factor(na_count_rate2$Taxon,levels = c(mean_rate_order$Taxon)),nrow = 5,scales ="free" )+
  #stat_summary(fun = mean,geom = "point",size=2.5,shape=17,col="black", position = position_dodge(0.750))+
  # stat_summary(fun = mean, geom = "text", col = "black",   
  #             vjust = 1.5,aes(label = paste("", round(..y.., digits = 2))))+
  stat_summary(fun.data = give.n, geom = "text", fun = median)
dev.off()


## marker-level missing rate

missrate_site <- c()
for (i in 1:nrow(d_sum)) {
  missrate_site [i] <-length(which(d_sum[i,] < 10))/ncol(d_sum)
}

mr_site <- data.frame(SNP=rownames(d_sum),missrate=missrate_site)
length(which(mr_site$missrate >=0.95))

boxplot(missrate_site)

d_ID <- fread("Pecan_microhaplo.botloci",header = F)
setdiff(subset(mr_site,mr_site$missrate>=0.95)$SNP,d_ID$V1)

##
counts <- matrix(nrow = nrow(d_sum))
for (i in 1:nrow(d_sum)) {
  counts[i]<- length(which(d_sum[i,] >=10))
}
count_d <- data.frame(counts)
rownames(count_d) <- rownames(d_sum)

length(which(counts < 10))

#f1
missrate_site_f1 <- c()
for (i in 1:nrow(d_sum_f1)) {
  missrate_site_f1[i] <-length(which(d_sum_f1[i,] < 10))/ncol(d_sum_f1)
}

mr_site_f1 <- data.frame(SNP=rownames(d_sum_f1),missrate=missrate_site_f1)
length(which(mr_site_f1$missrate >=0.95))

#dive
missrate_site_dive <- c()
for (i in 1:nrow(d_sum_dive)) {
  missrate_site_dive[i] <-length(which(d_sum_dive[i,] < 10))/ncol(d_sum_dive)
}

mr_site_dive <- data.frame(SNP=rownames(d_sum_dive),missrate=missrate_site_dive)
length(which(mr_site_dive$missrate >=0.95))

mr_site_all <- rbind.data.frame(mr_site,mr_site_f1,mr_site_dive)
mr_site_all$Population <- rep(c("Entire population", "F1 population","Diverse population"),each=3100)

mr_site_all$Population <- factor(mr_site_all$Population,levels = c("Entire population", "F1 population","Diverse population"))


#write.csv(mr_site_all,"missing rate at marker level_after samplemissrate-update.csv",row.names = F)

give.n<- function(y) {
  return( 
    data.frame(
      y =1.1*max(y),  #may need to modify this depending on your data
      label = paste('n = ', length(y), '\n',
                    'mean =', round(mean(y), 1), '\n'
      )
    )
  )
}


png("Marker_missrate_pecan.png",height = 2000,width = 3000,res = 300)
ggplot(mr_site_all,aes(x=mr_site_all$Population,y=mr_site_all$missrate*100,
                  fill=mr_site_all$Population))+geom_boxplot(width=0.5)+ 
  theme_bw()+ ylim(c(0,120))+
  scale_fill_manual(values = c("#fc8d62","#66c2a5","#8da0cb"))+
  #ggtitle("Marker-based missing rate")+
  theme(axis.text.x=element_text(angle=0,hjust=0.5),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(face = "bold"), axis.title.y = element_text(face = "bold"),
        legend.position = "none",text = element_text(size = 20))+
  xlab("Population")+ ylab("Missing rate (%)")+
  #facet_wrap(~factor(na_count_rate2$Taxon,levels = c(mean_rate_order$Taxon)),nrow = 5,scales ="free" )+
  #stat_summary(fun = mean,geom = "point",size=2.5,shape=17,col="black", position = position_dodge(0.750))+
  # stat_summary(fun = mean, geom = "text", col = "black",   
  #             vjust = 1.5,aes(label = paste("", round(..y.., digits = 2))))+
  stat_summary(fun.data = give.n, geom = "text", fun = median)
dev.off()

d <- read.csv("missing rate at marker level2_after samplemissrate.csv")
d$POS <- str_pad(d$POS,width = 9,side = "left",pad = "0")
d$SNP <- paste0(d$Chr,"_",d$POS)
write.csv(d,"missing rate at marker level2_after samplemissrate-update.csv")
