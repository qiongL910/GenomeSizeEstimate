FindTrustKmer<-function(file,kmer_len){
  # Identify the starting point of trusted k-mer
  for (q in 1:length(file$counts)){
    if (file[q,2]<file[q+1,2]){
      # Creat a variale hold the value of starting point of trusted k-mer frequency
      start_point<-q
      # Print out the result of the value of starting point of trusted k-mer frequency
      cat(paste("The trusted ",kmer_len,"-mers start at frequency of ",q,"\n\n",sep=""))
      # Stop at the first position that meets the argument
      break
    }
  }
  
  # Calculate the total amount of kmer in the dataset
  Total_kmer<-sum(as.numeric(file$frequency)*as.numeric(file$counts))
  # Calcualte the total amount of trusted and untrusted kmer in the dataset, as well as their percentage 
  untrusted_kmer<-sum(as.numeric(file$frequency[1:start_point-1])*as.numeric(file$counts[1:start_point-1]))
  portion_untrusted<-round(untrusted_kmer*100/Total_kmer,3)
  trusted_kmer<-sum(as.numeric(file$frequency[start_point:length(file$frequency)])*as.numeric(file$counts[start_point:length(file$counts)]))
  portion_trusted<-round(trusted_kmer*100/Total_kmer,3)
  # Create an empty vector to hold value
  kmer_report<-numeric(5)
  kmer_report[1:5]<-c(Total_kmer,trusted_kmer,portion_trusted,untrusted_kmer,portion_untrusted)
  # Add names for each value in the vector and return it
  names(kmer_report)<-c(paste("Total ",kmer_len,"-mer",sep=""),
                        paste("Trusted ",kmer_len,"-mer",sep=""),
                        paste("Trusted ",kmer_len,"-mer percentage",sep=""),
                        paste("Untrusted ",kmer_len,"-mer",sep=""),
                        paste("Untrusted ",kmer_len,"-mer percentage",sep=""))
  
  options("scipen"=100, "digits"=3)
  return(kmer_report)
  
}