GenomeEstimate<-function(file,kmer_len){
  # Get the start point of trusted kmer
  for (q in 1:length(file$counts)){
    if (file[q,2]<file[q+1,2]){
      # Creat a variale hold the value of starting point of trusted k-mer frequency
      start_point<-q
      # Stop at the first position that meets the argument
      break
    }
  }
  # Look for mean coverage (peak) starting from trusted kmer
  for (h in start_point:length(file$counts)){
    if(file[h,2]>file[h+1,2]){
      # Identify the position at the peak of bell shape and store it
      mean_coverage<-h
      # print out the mean coverage of kmer in trusted kmers
      cat(paste("The mean coverage of ",kmer_len,"-mer is ",mean_coverage,"\n",sep=""))
      break
    }
  }
  
  # Estimate genome size
  # Calculate the total number of trusted kmer
  trusted_kmer<-sum(as.numeric(file$frequency[start_point:length(file$frequency)])*as.numeric(file$counts[start_point:length(file$counts)]))
  # Estimated genome size from total number of trusted divided by mean kmer coverage 
  genome_size<-round(trusted_kmer/(mean_coverage*10^9),3)
  cat(paste("Estimated genome size based on ",kmer_len,"-mer analysis: ",genome_size," GB\n",sep=""))
  
  # Estimate single copy region
  # start at somewhere after peak 
  i<-mean_coverage+start_point
  subset_file<-file$counts[i:length(file$counts)]
  # Generate a decrease vector hold the counts different between each value
  decrease<-numeric(length(subset_file))
  for (h in 1:length(subset_file)){
    decrease[h]<-(subset_file[h]-subset_file[h+1])/subset_file[h]
  }
  # Evaluate the value in decrease vector based on 0.03 threshold
  for (k in 1:length(decrease)){
    if (decrease[k]<0.03){
      stop_point<-k
      break
    }
  }
  # create a variable to hold the stop point of single copy region and print out
  single_copy_stop<-i+stop_point
  cat(paste("Single copy region end point: ",single_copy_stop,"\n",sep=""))
  # Calculate the single copy region and its percentage in genome
  single_copy_region<-round(sum(as.numeric(file$frequency[start_point:single_copy_stop])*
                                  as.numeric(file$counts[start_point:single_copy_stop]))/(mean_coverage*1000000000),3)
  single_copy_percentage<-round(single_copy_region*100/genome_size,2)
  # Output result
  cat("Estimated single copy region:",single_copy_region," Gb\n",sep="")
  cat("Percentage of single copy region in genome: ",single_copy_percentage,"%\n",sep="")
}