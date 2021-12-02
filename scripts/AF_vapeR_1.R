library(vcfR)
library(ggplot2)
library(gggenes)
devtools::load_all("~/Exeter/afvaper/")

# Test inputs
vcf_test <- "~/Exeter/VCFs/FIBR_STAR_final_noGL_maf01.vcf.gz"
#vcf_input <- read.vcfR(vcf_test)

# Popmap
inds = system(paste0("bcftools query -l ",vcf_test),intern = T)
test_popmap <- data.frame(inds=inds)
test_popmap$pop <- gsub('[0-9]+', '', test_popmap$inds)

# Set dummies
test_vectors=list(c("GH","UL"),c("GH","LL"),c("GH","C"),c("GH","T"))
names(test_vectors) <- c("UL","LL","C","T")

# Global test variables
window_snps <- 200
core_N <- 6
permutations <- 10000
cutoff <- 0.99

# Get per chrom perms for a total of N permutations
chr_sizes <- read.table("~/Exeter/Genomes/STAR.chromosomes.release.fasta.fai")
chr_sizes <- chr_sizes[chr_sizes$V2 > 1000000,]
chr_props <- chr_sizes$V2/sum(chr_sizes$V2)
chr_perms <- ceiling(chr_props * permutations)

# Run over the following chr
chrs <- chr_sizes$V1
names(chr_perms) <- chrs 

# # Dummies
# vcf = vcf_input
# window_size = window_snps
# popmap = test_popmap
# vectors = test_vectors
# n_cores = core_N

#chr_res <- lapply(chrs,function(chr){
for(chr in chrs){
  
  message(paste0("STARTING CHR ",which(chrs == chr)," of ",length(chrs)))
  
  # Get VCF
  system(paste0("bcftools view -r ",chr," ",vcf_test," > outputs/tmp.vcf"))
  vcf_in <- read.vcfR("outputs/tmp.vcf")
  system("rm -f outputs/tmp.vcf")
  
  
  # Calculate allele frequency matrices
  AF_input <- afvaper::calc_AF_vectors(vcf = vcf_in,
                                       window_size = window_snps,
                                       popmap = test_popmap,
                                       vectors = test_vectors,
                                       n_cores = core_N)
  
  # Calculate null frequencies
  null_input <- afvaper::calc_AF_vectors(vcf = vcf_in,
                                         window_size = window_snps,
                                         popmap = test_popmap,
                                         vectors = test_vectors,
                                         n_cores = core_N,
                                         null_perms = chr_perms[chr])
  
  # Calculate Eigen statistics
  eigen_res <- lapply(AF_input,afvaper::eigen_analyse_vectors)
  
  # Get null matrix
  null_cutoffs <- afvaper::find_null_cutoff(null_res = null_input,
                                            cutoffs = c(0.95,0.99,0.999))
  
  # Test plot
  #eigenval_plot(eigen_res)[[1]]
  
  # Save both
  saveRDS(list(AF_input,null_input,eigen_res,null_cutoffs),
          paste0("outputs/fibr_",chr,"_AF_eigen_res_windsize_",window_snps,".rds"))
}

# Read in all the results...
all_chr_res <- lapply(chrs,function(x) readRDS(paste0("outputs/fibr_",x,"_AF_eigen_res_windsize_",window_snps,".rds")))

# Merge them all together
all_AF_vectors <- unlist(lapply(all_chr_res,'[[',1),recursive = F)
all_null_vectors <- unlist(lapply(all_chr_res,'[[',2),recursive = F)
all_eigen_res <- unlist(lapply(all_chr_res,'[[',3),recursive = F)

# Fetch cutoffs
all_eigen_cutoff <-  afvaper::find_null_cutoff(null_res = all_null_vectors,cutoffs = c(0.99,0.999,1))

# Test
genome_plots <- afvaper::eigenval_plot(eigen_res=all_eigen_res,null_vectors = all_null_vectors,plot.pvalues = T)
genome_plots[[2]]

# chr15
chr15_plots <- afvaper::eigenval_plot(eigen_res=all_eigen_res[grep("chr15",names(all_eigen_res),value=T)],null_vectors = all_null_vectors,plot.pvalues = T)
chr15_plots[[1]]
chr15_plots[[1]]+xlim(4000000,6000000)

# Calc and save pvals
fibr_pvals <- afvaper::eigen_pvals(all_eigen_res,all_null_vectors)
saveRDS(fibr_pvals,"~/Exeter/tmp/fibr_eigen_pvals.rds")


# Make Supp Tables of loading and signif windows --------------------------
# Get the signif windows
outliers <- afvaper::signif_eigen_windows(all_eigen_res,all_eigen_cutoff[,2])

# Eig1 summary
eig1_summary <- data.frame(rbindlist(lapply(outliers[[1]],afvaper::summarise_window_parallelism,eigen_res=all_eigen_res)))

# And fetch all the loadings
loading_pops <- c("LL","UL","C","T")
eig1_loadings <- matrix(nrow=length(outliers[[1]]),ncol=4)
colnames(eig1_loadings) <- loading_pops
for(i in 1:length(outliers[[1]])){
  for(pop in loading_pops){
    eig1_loadings[i,pop] <- all_eigen_res[[outliers[[1]][i]]]$eigenvecs[pop,1]
  }
}

# Eig2 summary
eig2_summary <- data.frame(rbindlist(lapply(outliers[[2]],afvaper::summarise_window_parallelism,eigen_res=all_eigen_res,eigenvector=2)))

# And fetch all the loadings
eig2_loadings <- matrix(nrow=length(outliers[[2]])*2,ncol=4)
colnames(eig2_loadings) <- loading_pops
outliers2 <- rep(outliers[[2]],each=2)
for(i in seq(1,length(outliers2),2)){
  for(pop in loading_pops){
    eig2_loadings[c(i,i+1),pop] <- all_eigen_res[[outliers2[i]]]$eigenvecs[pop,c(1,2)]
  }
}

# Write to outputs
eig1_out <- cbind(eig1_summary,eig1_loadings)
eig1_out <- eig1_out[order(-eig1_out$eigenvalue),]
write.table(eig1_out,
            "~/Exeter/fibr_simulations/tables/TableS16_eig1_windows.txt",
            quote = F,row.names = F,sep = "\t")

eig2_out <- cbind(eig2_summary,eig2_loadings)
eig2_out <- eig2_out[order(-eig2_out$eigenvalue_sum),]
write.table(eig2_out,
            "~/Exeter/fibr_simulations/tables/TableS17_eig2_windows.txt",
            quote = F,row.names = F,sep = "\t")


##### Chr15 region A matrices #####
# Get the signif windows
outliers <- afvaper:::signif_eigen_windows(all_eigen_res,all_eigen_cutoff[,3])
outlier_sum <- data.frame(rbindlist(lapply(outliers[[1]],summarise_window_parallelism,
                                           eigen_res = all_eigen_res,
                                           loading_cutoff = 0.3,
                                           eigenvector = 1)))

# Get focal regions
chr15_focal <- outlier_sum$window_id[39:55]
chr15_focal

# For each of these, get the A matrices
chr15_Amats <- data.frame(rbindlist(lapply(chr15_focal,function(wind){
  mat_tmp <- all_eigen_res[[wind]]$A_matrix
  out <- data.frame(Eig1=mat_tmp[,1],
                    pos=as.integer(sapply(strsplit(rownames(mat_tmp),"_"),'[[',2)))
})))

# And plot
Amat_fig <- ggplot(chr15_Amats,aes(pos,abs(Eig1)))+
  geom_point()+
  geom_hline(yintercept = c(0.3),linetype="dashed",colour="red2")+
  xlim(4950000,5500000)+
  theme_minimal()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=14))+
  labs(x="Chromosome 15 BP",y="Eigenvector 1 Score (Abs)")


# Fetch the genes
chr15_genes <- read.csv("data/fibr_chr15_genes.csv")
chr15_genes <- unique(chr15_genes[,1:4])
chr15_genes$gene_F <- factor(chr15_genes$gene,levels=unique(chr15_genes$gene))

# Remove vax1 and atp6v1ba which get dropped from figure anyway
chr15_genes <- chr15_genes[!(chr15_genes$gene %in% c("vax1","atp6v1ba")),]
chr15_gene_fig <- ggplot(chr15_genes, aes(xmin = start, xmax = end, y = gene_F)) +
  geom_gene_arrow(arrowhead_height = unit(0, "mm"), arrowhead_width = unit(0, "mm"),fill="orange") +
  theme_genes()+
  xlim(4950000,5500000)+
  theme(axis.title = element_blank(),
        axis.text.x = element_blank())

# Combine
cowplot::plot_grid(chr15_gene_fig,Amat_fig,
                   align="v",axis = "tblr",rel_heights = c(1,1.5),
                   ncol=1,nrow=2)

##### Chr8 region A matrices #####
chr8_plots <- afvaper::eigenval_plot(eigen_res=all_eigen_res[grep("chr8",names(all_eigen_res),value=T)],null_vectors = all_null_vectors,plot.pvalues = T)
chr8_plots[[1]]
chr8_plots[[1]]+xlim(24000000,26000000)

# Get the signif windows
outliers <- afvaper:::signif_eigen_windows(all_eigen_res,all_eigen_cutoff[,3])
outlier_sum <- data.frame(rbindlist(lapply(outliers[[1]],summarise_window_parallelism,
                                           eigen_res = all_eigen_res,
                                           loading_cutoff = 0.3,
                                           eigenvector = 1)))

# Get focal regions
chr8_focal <- grep("chr8:",outlier_sum$window_id,value = T)
chr8_focal

# For each of these, get the A matrices
chr8_Amats <- data.frame(rbindlist(lapply(chr8_focal,function(wind){
  mat_tmp <- all_eigen_res[[wind]]$A_matrix
  out <- data.frame(Eig1=mat_tmp[,1],
                    pos=as.integer(sapply(strsplit(rownames(mat_tmp),"_"),'[[',2)))
})))

# And plot
Amat_fig <- ggplot(chr8_Amats,aes(pos,(Eig1)))+
  geom_point()+
  geom_hline(yintercept = c(0.3),linetype="dashed",colour="red2")+
  theme_minimal()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=14))+
  labs(x="Chromosome 8 BP",y="Eigenvector 1 Score (Abs)")

# Combine
cowplot::plot_grid(chr8_gene_fig,Amat_fig,
                   align="v",axis = "tblr",rel_heights = c(1,1.5),
                   ncol=1,nrow=2)