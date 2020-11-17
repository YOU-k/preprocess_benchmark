library(lubridate)
setwd("/Users/you.y/Documents/preprocess/computation")

read.csv("cpu_drop.csv",header = TRUE) ->cpu_drop
as.character(cpu_drop$method) ->cpu_drop$method
cpu_drop$method[cpu_drop$method=="kallisto"]<-"kallisto_bustools"
cpu_all <- data.frame()
library(reshape)
for (i in c("a","b","c")){
  cpu_drop %>% filter(rep==i) -> new
  new$rep <- NULL
  methodnames <- new$method
  new <- new[,-1]
  matrix(period_to_seconds(hms(t(new))),nrow=ncol(new))/60 ->min
  colnames(min) <- methodnames
  min <- as.data.frame(min)
  min$reads <- c(8,30,70,100,130)*1000000
  melt(min,id="reads")->min
  min$rep=i
  rbind(cpu_all,min) -> cpu_all
}


library(reshape)

library(ggplot2)

#ggplot(cpu_drop,aes(x=reads,y=value,col=variable))+geom_point()+geom_line()


read.csv("walltime_drop.csv",header = TRUE) ->walltime_drop
as.character(walltime_drop$X) ->walltime_drop$X
walltime_drop$X[walltime_drop$X=="kallisto"]<-"kallisto_bustools"
walltime_all <- data.frame()
library(reshape)
for (i in c("a","b","c")){
  walltime_drop %>% filter(rep==i) -> new
  new$rep <- NULL
  methodnames <- new$X
  new <- new[,-1]
  matrix(period_to_seconds(hms(t(new))),nrow=ncol(new))/60 ->min
  colnames(min) <- methodnames
  min <- as.data.frame(min)
  min$reads <- c(8,30,70,100,130)*1000000
  melt(min,id="reads")->min
  min$rep=i
  rbind(walltime_all,min) -> walltime_all
}

library(plyr)
ddply(walltime_all, c("reads",  "variable"), summarise,
      mean = mean(value), sd = sd(value)) ->walltime_draw

library(ggplot2)
names(methods_colors)<-colornames
out.path <-"/Volumes/MattLab/Yue/preprocess/results/computation"

pdf("walltime_drop.pdf")
ggplot(walltime_draw,aes(x=reads,y=mean,col=variable))+geom_point()+geom_line()+theme_bw()+
  #theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_colour_manual(values = methods_colors) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05))+
  labs(y="run time(min)",x="number of reads",col="preprocess_pipelines")
dev.off()       



cpu_all$value <- cpu_all$value/walltime_all$value
ddply(cpu_all, c("reads",  "variable"), summarise,
      mean = mean(value), sd = sd(value)) ->cpu_draw

pdf("cpu_drop.pdf")
ggplot(cpu_draw,aes(x=reads,y=mean,col=variable))+geom_point()+geom_line()+theme_bw()+
  #theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_colour_manual(values = methods_colors) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05))+
  labs(y="CPU utilisation",x="number of reads",col="preprocess_pipelines")
dev.off() 

read.csv("mem_drop.csv",header = TRUE) ->mem_drop
as.character(mem_drop$X) ->mem_drop$X
mem_drop$X[mem_drop$X=="kallisto"]<-"kallisto_bustools"
mem_all <- data.frame()
library(reshape)
for (i in c("a","b","c")){
  mem_drop %>% filter(rep==i) -> new
  new$rep <- NULL
  methodnames <- new$X
  new <- new[,-1]
  t(new)->min
  colnames(min) <- methodnames
  min <- as.data.frame(min)
  min$reads <- c(8,30,70,100,130)*1000000
  melt(min,id="reads")->min
  min$rep=i
  rbind(mem_all,min) -> mem_all
}

library(plyr)
ddply(mem_all, c("reads",  "variable"), summarise,
      mean = mean(value), sd = sd(value)) ->mem_draw
pdf("mem_drop.pdf")
ggplot(mem_draw,aes(x=reads,y=mean,col=variable))+geom_point()+geom_line()+theme_bw()+
  #theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_colour_manual(values = methods_colors) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05))+
  labs(y="Memory(kb)",x="number of reads",col="preprocess_pipelines")
dev.off()


library(ggplot2)
library(cowplot)
as.numeric(as.character(mem_drop$value)) ->mem_drop$value
pdf(file.path(out.path,"drop.pdf"))

p1<-ggplot(mem_drop,aes(x=reads,y=value,col=variable))+geom_point()+geom_line()+theme_bw()+
  #theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_colour_manual(values = methods_colors) +
  labs(y="Memory(kb)",x="number of reads",col="preprocess_pipelines")

p2<-ggplot(walltime_drop,aes(x=reads,y=value,col=variable))+geom_point()+geom_line()+theme_bw()+
  #theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_colour_manual(values = methods_colors) +
  labs(y="run time(min)",x="number of reads",col="preprocess_pipelines")

p3<-ggplot(walltime_drop,aes(x=reads,y=cpu,col=variable))+geom_point()+geom_line()+theme_bw()+
  #theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_colour_manual(values = methods_colors) +
  labs(y="CPU utilisation",x="number of reads",col="preprocess_pipelines")
plot_grid(p1,p2,p3,nrow=3)
dev.off()
