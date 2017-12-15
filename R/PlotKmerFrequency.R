PlotKmerFrequency<-function(file,kmer_len,start_point,peak,end_point){
  # create a variable as plot title
  plot_title=substitute(paste("Genome size estimate based on ",italic("k"),"-mer anlysis (",kmer_len,"mer)",sep=""))
  # A variable to hold the value of single copy region
  singleC<-sum(as.numeric(file[start_point:end_point,1]*file[start_point:end_point,2]))/peak
  # Generate a data frame to draw segment line in plot
  x1<-c(peak)
  x2<-c(peak)
  y1<-c(0)
  y2<-c(file[peak,2])
  df<-data.frame(x1,x2,y1,y2)
  # this variable will be used for generating y axis label
  c<-as.integer(file[peak,2]/1000000)
  #print(c)
  
  # Using ggplot function to plot the distribution of kmer frequency
  ggplot(data=file[(start_point-8):(end_point+20),],mapping=aes_string(x='frequency',y='counts'))+
    # set the size of line
    geom_line(size=1)+
    # set the plot title
    labs(title=plot_title)+
    # color the trusted kmer area as blue
    geom_area(mapping=aes(x=ifelse(frequency>=start_point & frequency<=end_point,frequency,0)),fill="lightskyblue1",alpha=0.4)+
    # color the untrusted kmer area as red
    geom_area(mapping=aes(x=ifelse(frequency>=(start_point-8) & frequency<=start_point,frequency,0)),fill="coral",alpha=0.4)+
    # set the plot theme
    theme_classic(base_size=15)+
    # set x and y axis label
    labs(x=substitute(paste(italic("k"),"-mer frequency (k=",kmer_len,")",sep="")),
         y=substitute(paste(italic("k"),"-mer count",sep="")))+
    # add a poisson distribution in the plot indicating the theoretical distribution of kmer frequency
    geom_line(y=dpois((start_point-8):(end_point+20),peak)*singleC,lty=2,color="orangered",lwd=1)+
    # set y axis scale and tick label
    scale_y_continuous(limits=c(0,c*3*10^6),
                       breaks=seq(0,c*3*10^6,c*10^6),
                       labels=c("0M",paste(c,"M",sep=""),paste(2*c,"M",sep=""),paste(3*c,"M",sep="")))+
    # set x axis scale and tick label
    scale_x_continuous(breaks=c(0,start_point,peak,end_point,end_point+20))+
    # add a segment line indicating the mean coverage of kmer
    geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2),data=df,color="blue",lty=2,lwd=1)+
    # some parameters for plot style
    theme(plot.title = element_text(hjust = 0.5),
          plot.margin=unit(c(1,1,1,1),"cm"),
          axis.text.x=element_text(size=14,color="black"),
          axis.text.y=element_text(size=14,color="black"),
          axis.title=element_text(size=16))
}