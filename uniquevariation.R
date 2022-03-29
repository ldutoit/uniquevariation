
library("vcfR")


###Parameters
reference_list <- "background_list.txt"# one individual sample name per line
output_folder<- "output_unique" #need to exists, if not existing, create it.
vcf <- "populations.snps.vcf" # each individual in the vcf will be compared against a background list

#read vcf and background list
all_vcf<-read.vcfR(vcf)
background_inds<-readLines(reference_list)


background_vcf <-all_vcf@gt[,which(colnames(all_vcf@gt)%in%background_inds)]# a background of reference individuals
length(background_inds)==dim(background_vcf)[2] ## HAS TO BE TRUE

#Number of uniques individuals  
number_uniques<-cbind(colnames(all_vcf@gt),rep(0,dim(all_vcf@gt)[2]))

for (i in 2:dim(all_vcf@gt)[2]){
  ind_table=matrix(rep(NA,dim(background_vcf)[1]*3),ncol=3)
  colnames(ind_table)<-c("rare_allele_shared_with","common_allele_shared_with","unique")
  sub_ind =(all_vcf@gt)[,i]
  #sub_ind_gt<-c()
  print(colnames(all_vcf@gt)[i])
    for (j in 1:length(sub_ind)){
      if(j%%1000==0){print(c("SNP ",j))}
      allele1<-unlist(strsplit(unlist(strsplit(sub_ind[j],":"))[1],"/"))[1]
      allele2<-unlist(strsplit(unlist(strsplit(sub_ind[j],":"))[1],"/"))[2]
      if (allele1=="."){ 
        ind_table[j,1:3]<-c("NA","NA","NA")
        
      }else{
      if  (!(colnames(all_vcf@gt)[i]%in%background_inds)){ # if not in reference list, add your own ind alleles
        all_alleles_to_compare_to<-c(allele1,allele2)
      }
      if (colnames(all_vcf@gt)[i]%in%background_inds){  # ind in reference list, do not ind alleles
      all_alleles_to_compare_to<-c()
      }
      for (gen in background_vcf[j,] ){ # store all alles in backgroud
        all_alleles_to_compare_to[length(all_alleles_to_compare_to)+1]=unlist(strsplit(unlist(strsplit(gen,":"))[1],"/"))[1]
        all_alleles_to_compare_to[length(all_alleles_to_compare_to)+1]=unlist(strsplit(unlist(strsplit(gen,":"))[1],"/"))[2]
      }
      number_of_shared_allele_1<-length(which(allele1==all_alleles_to_compare_to))-1 # the -1 is for the self allele. if the ind is homozygote we only remove 1 still
      number_of_shared_allele_2<-length(which(allele2==all_alleles_to_compare_to))-1 # the -1 is for allele
      #unique?
      unique=0
      if(length(levels(as.factor((c(number_of_shared_allele_1,number_of_shared_allele_2)))))==1){
        if(levels(as.factor((c(number_of_shared_allele_1,number_of_shared_allele_2))))=="1"){
        unique=1
       }}
      if("0"%in%c(number_of_shared_allele_1,number_of_shared_allele_2)){unique=1}
      if(unique==1){
      print(c("unique, alleles shared with:",number_of_shared_allele_1,number_of_shared_allele_2))}
    ###populate table  
  ind_table[j,1:3]<-c(min(c(number_of_shared_allele_1,number_of_shared_allele_2)),max((c(number_of_shared_allele_1,number_of_shared_allele_2))),unique)
}}

write.table(ind_table,paste(output_folder,"/",colnames(all_vcf@gt)[i],"ind_unique.txt",sep=""),row.names=F,sep="\t",quote=F)
numeric_table<-apply(ind_table,c(1,2),as.numeric)
number_uniques[i,2]<-apply(numeric_table,2,sum,na.rm=T)[3]
}
#write the count of uniques to a table
colnames(number_uniques)<-c("ind","number_unique_pos")
number_uniques<-number_uniques[-1,]
write.table(number_uniques,paste(output_folder,"/number_uniques.txt",sep=""),row.names=F,sep="\t",quote=F)
