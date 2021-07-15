#Analyze the optimal alerts experiment data (first attempt)
remove(list = ls())
library(huxtable)
#library(foreign)
library(plyr)
library(dplyr)
library(xtable)
#library(ggpubr)
library(ggplot2)
library(haven)
library(tidyr)
library(sjlabelled)
library(httr)
library(survey)

tempdir()
dir.create(tempdir())


unloadNamespace("memisc")

#Put your own data folder here:
setwd("C:/Tornado_warnings/Experiment/Experiment_Analysis")
require("haven")

#httr::GET('cran.r-project.org/faqs.html')


#Importing the results:
fulldat<-
print(head(fulldat))
print(names(fulldat))

nrounds=6


names(fulldat)[names(fulldat) == "Sequence"] <- "sequence"

saveRDS(fulldat,file="firstpilot.Rdata")

#identify correct quiz answers
fulldat$blind_correct<-(fulldat$Q157==2)+(fulldat$Q158==3)
fulldat$informed_correct<-(fulldat$Q120==2)+(fulldat$Q119==2)+(fulldat$Q121==1)+(fulldat$Q118==3)+(fulldat$Q117==4)
table(fulldat$informed_correct)
fulldat$wtp_correct<-(fulldat$Q164==3)+(fulldat$Q166==3)
table(fulldat$wtp_correct)
fulldat$Ncorrect<-fulldat$blind_correct+fulldat$informed_correct+fulldat$wtp_correct
table(fulldat$Ncorrect)

unique.levels <- sort(unique(data$levels))
count <- table(fulldat$Ncorrect)
count.df <- data.frame(unique.levels, count)

setwd("C:/Tornado_warnings/Experiment")
#Understanding of the instructions:
COMPRplot<-fulldat %>% select(Ncorrect) %>% filter(complete.cases(.)) %>%
  ggplot(aes(Ncorrect))+
  geom_bar()+
  scale_fill_manual(values = c("navy"))+
  #coord_flip()+
  labs(x="Correct answers", y="N respondents", title = "N of correct quiz responses")
  #scale_y_continuous(limits=c(0,1),labels = scales::percent)+
  #theme(text = element_text(size=15), plot.title = element_text(size=15), legend.position = "none")+
  #scale_x_discrete(limits = rev(levels(fulldat$Q117)))
print(COMPRplot)
ggsave("Graphs/COMPRplot.pdf")

#Self-reported understanding:
und_vars<-c("Q107_1","Q107_2","Q107_3","Q107_4")
task_names<-c("Blind_Protection","Informed_Protection","Belief_Elicitation","Value_Elicitation")

fulldat %>% select(und_vars) %>%filter(complete.cases(.)) %>%
  sapply(function(x) ifelse(x==9,1,0)) %>% as.matrix() %>% colSums %>% as.matrix ->nagree
resp<-data.frame(task_names,nagree)
resp$pagree<-resp$nagree/19

uplot1<-resp %>%
  ggplot(aes(x=task_names,y=pagree, fill="blue"))+
  geom_bar(stat='identity')+
  scale_fill_manual(values = c("navy"))+
  scale_y_continuous(limits=c(0,1),labels = scales::percent)+
  labs(x="Task", y="% agree", title = "I understood the instructions:")+
  theme(text = element_text(size=15), plot.title = element_text(size=15), legend.position = "none")+
  coord_flip()
print(uplot1)
ggsave("Graphs/Uplot1.pdf")


und_vars_payoff<-c("Q112_1","Q112_2","Q112_3","Q112_4")
fulldat %>% select(und_vars_payoff) %>%filter(complete.cases(.)) %>%
  sapply(function(x) ifelse(x==9,1,0)) %>% as.matrix() %>% colSums %>% as.matrix ->nagree
resp2<-data.frame(task_names,nagree)
resp2$pagree<-resp2$nagree/19


uplot2<-resp2 %>%
  ggplot(aes(x=task_names,y=pagree, fill="blue"))+
  geom_bar(stat='identity')+
  scale_fill_manual(values = c("navy"))+
  scale_y_continuous(limits=c(0,1),labels = scales::percent)+
  labs(x="Task", y="% agree", title = "I understood how my choices affected the payoff in:")+
  theme(text = element_text(size=15), plot.title = element_text(size=15), legend.position = "none")+
  coord_flip()
print(uplot2)
ggsave("Graphs/Uplot2.pdf")


#Filtering only the responses with good understanding of instructions:
fulldat$goodquiz<-(fulldat$Ncorrect>6)
fulldat<-subset(fulldat,goodquiz==TRUE)

#merge treatment sequences:
setwd("C:/Tornado_warnings/Experiment/")
exp_treatments<-readRDS(file="exp_treatments.Rdata") #reading the 

names(exp_treatments)<-c("treatN","belong","p","tot_gr","bl_gr","w_gr")
exp_treatments<-subset(exp_treatments,select= -c(belong))
write.csv(exp_treatments,'exp_treatments_pilot.csv')



treat_sequences<-read.csv(file='treatment_sequences.csv')
names(treat_sequences)<-c("sequence","r1","r2","r3","r4","r5","r6")

#treatments dataframe has a panel structure: sequence and round in two first columns
#then a prior prob, nhonest, nwhite, nblack gremlins, then posterior white, posterior black, WTP

treat_sequences %>% gather(round, treatN, r1:r6, factor_key=TRUE) %>% 
  mutate(round=as.numeric(substr(round,2,2))) %>% arrange(sequence, round) ->treat_sequences

treatments<-merge(treat_sequences,exp_treatments,by=c("treatN"),all.x = TRUE)

treatments %>% mutate(phintWB=(w_gr/tot_gr), phintWW=(tot_gr-bl_gr)/tot_gr) %>%
  mutate(phintBB=1-phintWB, phintBW=1-phintWW) %>%
  mutate(post_probW=p*phintWB/(p*phintWB+(1-p)*phintWW), post_probB=p*phintBB/(p*phintBB+(1-p)*phintBW)) -> treatments
print(treatments)


#Saving the text information
fileConn<-file("Feedback2.txt","w")
writeLines(c("Which task did you find the most confusing?"), fileConn)
writeLines(paste(fulldat$Q104, collapse = "\n"), fileConn)
writeLines(c("**************************************************\n \n"), fileConn)
writeLines(c("Please explain the strategy you used for Task 2 (Informed Protection)? This is the task in which you see a hint and then decide to protect or not."), fileConn)
writeLines(paste(fulldat$Q105, collapse = "\n"), fileConn)
writeLines(c("**************************************************\n \n"), fileConn)
writeLines(c("Please explain the strategy you used in Task 3 (Measuring Chances)? This is the task in which you use sliders."), fileConn)
writeLines(paste(fulldat$Q106, collapse = "\n"), fileConn)
writeLines(c("**************************************************\n \n"), fileConn)
writeLines(c("Do you see any ways to improve the experiment?"), fileConn)
writeLines(paste(fulldat$Q114, collapse = "\n"), fileConn)
writeLines(c("**************************************************\n \n"), fileConn)
writeLines(c("Anything else you would like to say about the experiment?"), fileConn)
writeLines(paste(fulldat$Q115, collapse = "\n"), fileConn)
close(fileConn)


#Prepare the data on blind protection choices:
blprotN<-c("1_B1","2_B1", "3_B1", "4_B1", "5_B1")
#create better names for informed protection variables
IP_W<-names(fulldat)[substr(names(fulldat),2,6)=="_IP_W"]
IP_B<-names(fulldat)[substr(names(fulldat),2,6)=="_IP_B"]

print(fulldat[,blprotN])


#Prepare the data on informed protection choices
IPdata<-fulldat[,c("Participant_ID","sequence",IP_W,IP_B)]


h=1
for (s in c("W","B")) {
  for (r in 1:nrounds) {
    old_choicename=sprintf("%d_IP_%s",r,s)
    print(old_choicename)
    new_choicename=sprintf("IPchoice_%d_%d", h,r)
    names(IPdata)[names(IPdata) == old_choicename] <-sprintf("IP_%d_%d", h,r)
  
  }
  h=h+1
  
}
print(names(IPdata))
print(lapply(IPdata,class))



IPdata<-IPdata[complete.cases(IPdata),]
IPdata %>% 
  pivot_longer(
    cols=IP_1_1:IP_2_6,
    names_to = c(".value","Signal","round"),
    names_pattern="(IP)_(.)_(.)",
    values_drop_na = TRUE)->IPtransf

IPtransf<-merge(IPtransf,treatments,by=c("sequence","round"),all.x=TRUE)
IPtransf%>%mutate(post_prob=post_probW*(2-as.numeric(Signal))+post_probB*(as.numeric(Signal)-1)) -> IPtransf




#Preparing reported beliefs data:
belief_repn<-names(fulldat)[substr(names(fulldat),2,6)=="_BE1_"] #select variable names

belief_dat<-fulldat[,c("Participant_ID","sequence",belief_repn)] #select the data

#reformatting variable names to BE_(hint color)_(round)
h=1
for (s in c("W","B")) {
  for (r in 1:nrounds) {
    new_choicename=sprintf("IPchoice_%d_%d", h,r)
    names(belief_dat)[names(belief_dat) == sprintf("%d_BE1_%s_1",r,s)] <-sprintf("BE_%d_%d", h,r)
    names(belief_dat)[names(belief_dat) == sprintf("%d_BE1_%s_15",r,s)] <-sprintf("BE_%d_%d", h,r)
  }
  h=h+1
}
belief_dat<-belief_dat[complete.cases(belief_dat),]

#reshaping to a long format (panel with sequence, round in columns)
belief_dat %>% 
  pivot_longer(
    cols=BE_1_1:BE_2_6,
    names_to = c(".value","Signal","round"),
    names_pattern="(BE)_(.)_(.)",
    values_drop_na = TRUE)->belief_dat

belief_dat<-merge(belief_dat,treatments,by=c("sequence","round"),all.x=TRUE) #merging round characteristics
belief_dat%>%mutate(post_prob=post_probW*(2-as.numeric(Signal))+post_probB*(as.numeric(Signal)-1), belief=0.01*BE) -> belief_dat #calculate posterior prob
print(belief_dat)


#Preparing WTP data
WTPmatrix<-matrix(0, nrow=dim(fulldat)[1], ncol = nrounds)

for (r in 1:6){
  
  WTPnames<-names(fulldat)[substr(names(fulldat),1,6)==sprintf("%d_WTP1",r)] 
  print(WTPnames)
  fulldat %>% select(WTPnames) %>% sapply(function(x) ifelse(x==-99,0,x)) %>% as.matrix %>% rowSums->WTP
  WTPmatrix[,r]<-0.5*(WTP-1)
  
}
print(WTPmatrix)
WTPdat<-data.frame(fulldat$Participant_ID, fulldat$sequence,WTPmatrix)
names(WTPdat)<-c("Participant_ID","sequence","r_1","r_2","r_3","r_4","r_5","r_6")

WTPdat %>% filter(complete.cases(.)) %>% pivot_longer(
  cols=r_1:r_6,
  names_to = c(".value","round"),
  names_pattern="(r)_(.)",
  values_drop_na = TRUE) ->WTPdat
names(WTPdat)<-c("Participant_ID","sequence","round","WTP")

WTPdat<-merge(WTPdat,treatments,by=c("sequence","round"),all.x=TRUE)

#Calculating the signal's value for a risk-neutral agent
Loss=20
cost=5

WTPdat %>% mutate(Value=-(p*phintWB*Loss-(p*phintBB+(1-p)*phintBW)*cost))->WTPdat
print(table(WTPdat$Value))



#Protection choices based on probability (boxplot):
BlChoiceNames<-names(fulldat)[endsWith(names(fulldat), '_B1')]

fulldat[,BlChoiceNames] %>%  pivot_longer(
  cols='1_B1':'5_B1',
  names_to = c("round",".value"),
  names_sep = "_",
  #names_pattern="(IPchoice|IPhint)_(.)_(.)",
  values_drop_na = TRUE)->Blind_Dat

Blind_Dat$p<-0.1*as.numeric(Blind_Dat$round)

#aggregate to calculate trials and successes
Blind_Dat %>% group_by(p) %>% dplyr::summarize(average = mean(B1), totprotect=sum(B1), n = n()) ->props

#calculate confidence intervals
props_CI<-mapply(function(x,y) prop.test(x,y, alternative = "two.sided", conf.level = 0.95)$conf.int, props$totprotect,props$n)
vnames<-names(props)


props<-bind_cols(props,as.data.frame(t(props_CI)))#add confidence bounds back
names(props)<-c(vnames,"L_bound","U_bound")

props$Prob_black<-factor(props$p,levels = c(0.1,0.2,0.3,0.4,0.5),labels = c("0.1","0.2","0.3","0.4","0.5"),ordered=TRUE)
#making the graph

BL_plot<-props %>%
  ggplot(aes(x=Prob_black,y=average, fill="blue"))+
  geom_crossbar(aes(ymin=L_bound,ymax=U_bound),width=0.5)+
  scale_fill_manual(values = c("grey"))+
  labs(x="Probability of Black Ball", y="% protecting", title = "Blind Protection vs Probability")+
  theme(legend.position = "none")+
  theme(text = element_text(size=18))+
  scale_y_continuous(limits=c(0,1),labels = scales::percent)
print(BL_plot) 
ggsave("Graphs/BLProt_plot2.pdf")




#Understanding of the gremlins:
fulldat$Q117<-factor(fulldat$Q117,levels = c(1,2,3,4),labels = c("The Ball is white","The Ball is black", "We don't know", "This is impossible"),ordered=TRUE)
fulldat$Q118<-factor(fulldat$Q118,levels = c(1,2,3,4),labels = c("The Ball is white","The Ball is black", "We don't know", "This is impossible"),ordered=TRUE)

GRdat<-data.frame(fulldat[,c("Q117","Q118")])
GRdat<-GRdat[complete.cases(GRdat),]
GRplot1<-GRdat %>%
  ggplot()+
  geom_bar(mapping = aes(x = Q117, y=..prop..,group=1,fill="blue"), stat="count", na.rm=TRUE)+
  scale_fill_manual(values = c("navy"))+
  coord_flip()+
  labs(x="", y="% mentioned", title = "If a black-swamp gremlin says\n that the ball is white__")+
  scale_y_continuous(limits=c(0,1),labels = scales::percent)+
  theme(text = element_text(size=15), plot.title = element_text(size=15), legend.position = "none")+
  scale_x_discrete(limits = rev(levels(fulldat$Q117)))
print(GRplot1)

GRplot2<-GRdat %>%
  ggplot()+
  geom_bar(mapping = aes(x = Q118, y=..prop..,group=1,fill="blue"), stat="count", na.rm=TRUE)+
  scale_fill_manual(values = c("navy"))+
  coord_flip()+
  labs(x="", y="% mentioned", title = "If a white-swamp gremlin says\n that the ball is white__")+
  scale_y_continuous(limits=c(0,1),labels = scales::percent)+
  theme(text = element_text(size=15), plot.title = element_text(size=15), legend.position = "none")+
  scale_x_discrete(limits = rev(levels(fulldat$Q117)))
print(GRplot2)
library(ggpubr)
fig_response_sd <- ggarrange(GRplot1,GRplot2,widths=c(2,2), common.legend = TRUE,legend="none")+
  theme(legend.title=element_blank()) 
print(fig_response_sd)
ggsave("Graphs/GR_responses2.pdf",width = 10, height = 6)
#Q117, Q118

#Filtering data based on responses
fulldat$correct<-(fulldat$Q117=="This is impossible")&(fulldat$Q118=="We don't know")


#Protection choices based on the signal received #calculate average protection rate based on the posterior
IPtransf %>% group_by(post_prob) %>% summarize(average = mean(IP), totprotect=sum(IP), n = n()) ->IPprops

IP_plot<-IPprops %>% ggplot(aes(post_prob,average),color = "#1F3552")+
  geom_line(color = "#1F3552",size=2)+
  geom_point(color = "#1F3552",size=4)+
  labs(x="Posterior probability", y="% protecting", title = "Protection response to signals/hints")+
  geom_vline(xintercept=0.25, linetype="dashed", color = "darkgray", size=1)+
  theme(text = element_text(size=18))+
  annotate("text", x = 0.32, y = 0.05, label = "Cost/loss ratio")
print(IP_plot)
ggsave("Graphs/IP_plot.pdf")






#Scatter plot reported beliefs vs posteriors
UPD_curve<-ggplot(belief_dat, aes(x=post_prob, y=belief)) + geom_point(colour = 'dodgerblue4',fill="dodgerblue4", size=3, shape=21)+
  labs(x="Posterior Probability",y="Reported Belief", title = "Bayesian Updating")+
  scale_x_continuous(limits=c(0,1),labels = scales::percent)+
  scale_y_continuous(limits=c(0,1),labels = scales::percent)+
  theme(text = element_text(size=18))+
  geom_jitter(width = 0.1, height = 0.1, colour='dodgerblue4')
  #geom_text(hjust=0.3,vjust=-0.7,size=5)

print(UPD_curve)
ggsave("Graphs/UPD_curve2.pdf")



#Scatter plot: WTP vs signal value
WTP_curve2<-ggplot(WTPdat, aes(x=Value, y=WTP)) + geom_jitter(width=.05, colour = 'dodgerblue4',fill="dodgerblue4", size=3, shape=21)+
  labs(x="Value",y="Willingness-to-pay (USD)", title = "WTP for a signal vs predicted value")+
  scale_x_continuous(limits=c(-0.5,5))+
  scale_y_continuous(limits=c(-0.5,5))+
  theme(text = element_text(size=18))
#geom_text(hjust=0.3,vjust=-0.7,size=5)

print(WTP_curve2)
ggsave("Graphs/WTP_curve2.pdf")

