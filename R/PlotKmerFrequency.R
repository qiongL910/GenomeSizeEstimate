PlotKmerFrequency<-function(file,kmer_len,start_point,peak,end_point){
  plot_title=substitute(paste("Genome size estimate based on ",italic("k"),"-mer anlysis (",kmer_len,"mer)",sep=""))
  singleC<-sum(as.numeric(file[start_point:end_point,1]*file[start_point:end_point,2]))/peak
  # Generate a data frame to draw segment line in plot
  df<-data.frame(x1=peak,x2=peak,y1=0,y2=file[peak,2])
  # this value is used for generating y axis label
  c<-as.integer(file[peak,2]/1000000)
  #print(c)
  
  ggplot(data=file[(start_point-8):(end_point+20),],mapping=aes(x=frequency,y=counts))+geom_line(size=1)+
    labs(title=plot_title)+
    geom_area(mapping=aes(x=ifelse(frequency>=start_point & frequency<=end_point,frequency,0)),fill="lightskyblue1",alpha=0.4)+
    geom_area(mapping=aes(x=ifelse(frequency>=(start_point-8) & frequency<=start_point,frequency,0)),fill="coral",alpha=0.4)+
    theme_classic(base_size=15)+
    labs(x=substitute(paste(italic("k"),"-mer frequency (k=",kmer_len,")",sep="")),
         y=substitute(paste(italic("k"),"-mer count",sep="")))+
    geom_line(y=dpois((start_point-8):(end_point+20),peak)*singleC,lty=2,color="orangered",lwd=1)+
    scale_y_continuous(limit=c(0,c*3*10^6),
                       breaks=seq(0,c*3*10^6,c*10^6),
                       labels=c("0M",paste(c,"M",sep=""),paste(2*c,"M",sep=""),paste(3*c,"M",sep="")))+
    scale_x_continuous(breaks=c(0,start_point,peak,end_point,end_point+20))+
    geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2),data=df,color="blue",lty=2,lwd=1)+
    theme(plot.title = element_text(hjust = 0.5),
          plot.margin=unit(c(1,1,1,1),"cm"),
          axis.text.x=element_text(size=14,color="black"),
          axis.text.y=element_text(size=14,color="black"),
          axis.title=element_text(size=16))
}