sacc_accuracy <- read.csv('sacc_accuracy.csv')
sacc_latency <- read.csv('sacc_latency.csv')
sacc_pvel <- read.csv('sacc_pvel.csv')

sacc <-melt(sacc_pvel, id = "pid", measured = c("pvel_c1", "pvel_c2", "pvel_c3"))
names(sacc)<-c("Participant", "Condition", "Pvel")
sacc$Condition<-factor(sacc$Condition, labels = c("C1", "C2", "C3"))
sacc<-sacc[order(sacc$Participant),]

by(sacc$Pvel, sacc$Condition, stat.desc, basic = FALSE)

#Set Contrasts
C1vsC2C3<-c(1, -0.5, -0.5)
C2vsC3<-c(0, -1, 1)
contrasts(sacc$Condition)<-cbind(C1vsC2C3, C2vsC3)

baseline<-lme(Pvel ~ 1, random = ~1|Condition/Pvel, data = sacc, method = "ML")
saccModel<-lme(Pvel ~ Condition, random = ~1|Participant/Condition, data = sacc, method = "ML")
summary(baseline)
summary(saccModel)
print(anova(baseline, saccModel))

lowerCI <- c(1.8242752-0.6048310, 5.4095274-1.9110249, 8.2451789-1.9369098)
upperCI <- c(1.8242752+0.6048310, 5.4095274+1.9110249, 8.2451789+1.9369098)
ci3 <- 1.9369098 

saccBar <- ggplot(sacc, aes(Condition, TargetOff))
saccBar + stat_summary(fun.y=mean, geom ='bar', fill=c("#E69F00", "#56B4E9", "#009E73"), colour=c("#000000", "#000000", "#000000")) + 
stat_summary(fun.data="mean_c1_boot", colour="red")



