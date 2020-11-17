methods_colors <-c("#E69F00" ,"#56B4E9"  , "#009E73"  ,"#F0E442"  ,  "#0072B2","#D55E00"  , "#CC79A7"  , "yellowgreen"  , "#91D1C2B2" )
names(methods_colors) <- c("scpipe","celseq2","scruff","optimus","dropseqpipe","zumis","cellranger","alevin","kallisto_bustools")
library(lubridate)
setwd("Documents/preprocess/computation/")
library(tidyverse)
read.csv("cpu_plate.csv",header = TRUE) ->cpu_plate
as.character(cpu_plate$method) ->cpu_plate$method
cpu_plate$method[cpu_plate$method=="kallisto"]<-"kallisto_bustools"

cpu_all <- data.frame()
library(reshape)
for (i in c("a","b","c")){
  cpu_plate %>% filter(rep==i) -> new
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

library(plyr)

library(ggplot2)

#ggplot(cpu_plate,aes(x=reads,y=value,col=variable))+geom_point()+geom_line()


read.csv("walltime_plate.csv",header = TRUE) ->walltime_plate
as.character(walltime_plate$X) ->walltime_plate$X
walltime_plate$X[walltime_plate$X=="kallisto"]<-"kallisto_bustools"

walltime_all <- data.frame()
library(reshape)
for (i in c("a","b","c")){
  walltime_plate %>% filter(rep==i) -> new
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


library(reshape)

library(ggplot2)

out.path <-"/Volumes/MattLab/Yue/preprocess/results/computation"
#setwd(out.path)
pdf("walltime_plate.pdf",width = 8,height = 6)
library(scales)
ggplot(walltime_draw,aes(x=reads,y=mean,col=variable))+geom_point()+geom_line()+theme_bw()+
  #theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_colour_manual(values = methods_colors) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),width = 0)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides="b")  + theme(text=element_text(size=20)) +
  labs(y="run time(min)",x="number of reads",col="preprocess_workflows")
dev.off()       
       


cpu_all$value <- cpu_all$value/walltime_all$value
ddply(cpu_all, c("reads",  "variable"), summarise,
      mean = mean(value), sd = sd(value)) ->cpu_draw

pdf("cpu_plate.pdf",width = 8,height = 6)
ggplot(cpu_draw,aes(x=reads,y=mean,col=variable))+geom_point()+geom_line()+theme_bw()+
  #theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_colour_manual(values = methods_colors) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),width = 0)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides="b")  + theme(text=element_text(size=20))+
  labs(y="CPU utilisation",x="number of reads",col="preprocess_workflows")
dev.off() 

read.csv("mem_plate.csv",header = TRUE) ->mem_plate
as.character(mem_plate$X) ->mem_plate$X
mem_plate$X[mem_plate$X=="kallisto"]<-"kallisto_bustools"

mem_all <- data.frame()
library(reshape)
for (i in c("a","b","c")){
  mem_plate %>% filter(rep==i) -> new
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
pdf("mem_plate.pdf",width = 8,height = 6)
ggplot(mem_draw,aes(x=reads,y=mean,col=variable))+geom_point()+geom_line()+theme_bw()+
  #theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_colour_manual(values = methods_colors) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),width = 0)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides="lb") + theme(text=element_text(size=20)) +
  labs(y="Memory(kb)",x="number of reads",col="preprocess_workflows")
dev.off()








library(ggplot2)
library(cowplot)
as.numeric(as.character(mem_plate$value)) ->mem_plate$value
pdf(file.path(out.path,"plate.pdf"))

p1<-ggplot(mem_draw,aes(x=reads,y=value,col=variable))+geom_point()+geom_line()+theme_bw()+
  #theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_colour_manual(values = methods_colors) +
  labs(y="Memory(kb)",x="number of reads",col="preprocess_pipelines")

p2<-ggplot(walltime_draw,aes(x=reads,y=value,col=variable))+geom_point()+geom_line()+theme_bw()+
  #theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_colour_manual(values = methods_colors) +
  labs(y="run time(min)",x="number of reads",col="preprocess_pipelines")

p3<-ggplot(cpu_draw,aes(x=reads,y=cpu,col=variable))+geom_point()+geom_line()+theme_bw()+
  #theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_colour_manual(values = methods_colors) +
  labs(y="CPU utilisation",x="number of reads",col="preprocess_pipelines")
plot_grid(p1,p2,p3,nrow=3)
dev.off()
