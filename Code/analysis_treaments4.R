#R analysis of the experimental data

# Install We Anderson palettes
#install.packages("wesanderson")
# Load


remove(list = ls())
library(dplyr)
library(xtable)
library(ggplot2)
library(haven)
library(tidyr)
library(sjlabelled)
library(qualtRics)
library(magrittr)
library(wesanderson)

tempdir()
dir.create(tempdir())
setwd("C:/Tornado_warnings/Experiment")

#Signal descriptions:
totgremlins<-c(2,3,3,3,5,5,5)
blsw_gremlins<-c(0,1,0,1,1,0,1)
wsw_gremlins<-c(0,0,1,1,0,1,1)

#Probabilities of black balls
treatm_probs<-c(0.1,0.2,0.3,0.5)


nsignals=length(totgremlins)
nprobs=length(treatm_probs)


ntreatments=nsignals*nprobs
snames=sapply(c(1:ntreatments),as.character)

exp_treatments<-data.frame(snames=as.numeric(snames), belong=rep(1,each=ntreatments), p=rep(treatm_probs,each=nsignals), tot_gr=rep(totgremlins,times=nprobs), bl_gr=rep(blsw_gremlins,times=nprobs), w_gr=rep(wsw_gremlins,times=nprobs))

print(exp_treatments)

#Save to import into the data analysis:
saveRDS(exp_treatments, file = "exp_treatments.Rdata")
write.csv(exp_treatments,'exp_treatments_pilot.csv')

exp_treatments$snames<-rep('',each=ntreatments)
exp_treatments %<>% mutate(tpr=(tot_gr-w_gr)/tot_gr,fpr=bl_gr/tot_gr, ppv=p*(tot_gr-w_gr)/(p*(tot_gr-w_gr)+(1-p)*bl_gr))
#exp_treatments %<>% mutate(ppv=(1-p)*fpr/(p*tpr+(1-p)*fpr), fdr=(1-p)*fpr/(p*tpr+(1-p)*fpr))
print(exp_treatments)

real_signals<-data.frame(snames=c("COVID PCR","Tornado Warning", "Screening Mammography"), belong=rep(0,each=3), p=c(0.05,0.1,0.0086),  tpr=c(0.7, 0.7, 0.92), fpr=c(0.02, 0.25, 0.17))

print(real_signals)

signals<-rbind(exp_treatments[,c("snames","belong","p","tpr","fpr")],real_signals)
signals %<>% mutate(ppv=(1-p)*fpr/(p*tpr+(1-p)*fpr), fdr=(1-p)*fpr/(p*tpr+(1-p)*fpr))

print(signals)

#signals %>% filter(snames%in% c("I","II","III","IV","COVID PCR","Tornado Warning")) %>%
#  mutate(snames=recode(snames,"I"="I/V","II"="II/VI","III"="III/VII","IV"="IV/VIII"))->signals2

#Graph of treatments (ROC curve)
ROC_curve<-ggplot(signals, aes(x=fpr, y=tpr, label=snames, fill=factor(belong))) + geom_point(size=3, shape=21)+
  labs(x="False positive rate",y="True positive rate", title = "Experiment Treatments vs Other Signals: ROC")+
  scale_x_continuous(limits=c(0,0.55),labels = scales::percent)+
  scale_y_continuous(limits=c(0.5,1.05),labels = scales::percent)+
  #scale_color_manual(values=c("#E69F00", "#56B4E9"))+
  scale_fill_manual(values=wes_palette(n=2, name="Royal1"))+
  theme(text = element_text(size=18), legend.position = "none")+
  geom_text(hjust=-0.1,vjust=0.9,size=5)

print(ROC_curve)
ggsave("ROC_curve2.png")


#Graph of treatments (detection error tradeoff)
PR_curve<-ggplot(signals, aes(x=tpr, y=ppv, label=snames, fill=factor(belong))) + geom_point(size=3, shape=21)+
  labs(x="True Positive Rate (Recall)",y="Precision (PPV)", title = "Experiment Treatments vs Other Signals (recall-precision)")+
  scale_x_continuous(limits=c(0.6,1.05),labels = scales::percent)+
  scale_y_continuous(limits=c(0,1.05),labels = scales::percent)+
  scale_fill_manual(values=wes_palette(n=2, name="Royal1"))+
  theme(text = element_text(size=18), legend.position = "none")+
  geom_text(hjust=0.3,vjust=-0.7,size=5)

print(PR_curve)
ggsave("PR_curve2.pdf")






signals$fnr<-1-signals$tpr

signals %>% filter(snames%in% c("I","II","III","IV","VI","VIII","COVID PCR","Tornado Warning")) %>%
  mutate(snames=recode(snames,"I"="I/V","III"="III/VII"))->signals2

#Graph of treatments (false alarm rate vs false negative rate)
FAR_curve<-ggplot(signals, aes(x=fnr, y=fdr, label=snames, fill=factor(belong))) + geom_point(size=3, shape=21)+
  labs(x="False Negative Rate=Miss Rate",y="False Alarm Rate", title = "Treatments vs Other Signals: Performance Diagram")+
  scale_x_continuous(limits=c(0,0.5),labels = scales::percent)+
  scale_y_continuous(limits=c(0,1),labels = scales::percent)+
  scale_fill_manual(values=wes_palette(n=2, name="Royal1"))+
  theme(text = element_text(size=18), legend.position = "none")+
  geom_text(hjust=0.3,vjust=-0.7,size=5)

print(FAR_curve)
ggsave("FAR_curve2.pdf")
