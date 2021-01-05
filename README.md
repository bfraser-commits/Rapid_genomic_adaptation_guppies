# FIBR
To accompany FIBR paper

This repository contains various scripts that were used for population genomics analyses in  
**Title of paper**

Raw reads were processed with a set of in-house scripts published in Fraser et al. (2020) - Improved reference genome uncovers novel sex-linked regions in the guppy (Poecilia reticulata)

VCFs were phased using the following script:  
```phase_vcf.sh```  
For this you will need the final filtered bam files.  

Download the scripts into a folder (FIBR_ms) with the following sub-folders: 
 1. scripts (this is where the scripts go)  

 2. data  
  - popgenome  
  - heterozygosity  
  - allele_freq  
  - ROH  
  - PCA  
  - selscan  
    - xpehh  
    - ihh12  

 3. output  
  - popgenome  
  - heterozygosity  
  - allele_freq  
  - ROH  
  - PCA  
  - selscan  
    - xpehh  
    - ihh12 
    
 4. figures  
  - popgenome  
  - heterozygosity  
  - allele_freq  
  - ROH  
  - PCA  
  - selscan  
    - xpehh  
    - ihh12 

### Processing per analysis  
#### Popgenome  
**To obtain raw data**  
Run these scripts on the HPC to get Fst, nuc.div. and Tajima's D:  
FST - ```popgenome.R```  
Tajima'S - ```popgenome.R```  
Nucleotide diversity (π) - ```popgenome.R```  
These scripts can be run by this shell:  ```popgenome.sh```  

**To process the raw data**  
On the HPC, in the folder with the outputs, run the following to concatenate the files per statistic:  
```
for stat in fst pi td; 
   do     
   for file in *${stat}*out;
      do
      awk 'FNR==1 && NR!=1 { while (/^chrom+/) getline; } 1 {print}' *${stat}*out > FIBR_${stat}.txt;
      done;
   done
 ```  
 Now move the output files to the main folder FIBR_ms/data/popgenome  
 
 To get FST means and medians run: ```fst_global.R```  
 To process Tajima's D run: ```TD_windows.R```  
 To process π run: ```Pi_windows.R```  

#### VCFtools  
**To obtain raw data**  
Run these scripts on the HPC:  
Heterozygosity - ```heterozygosity_vcftools.sh```  
Allele frequency - ```allele_freq_vcftools.sh```  

**To process the data**  
*Heterozygosity*  
On the HPC, in the folder with the raw output, run the following:  
```
for POP in GHP GLP IC IT ILL IUL
   do
   cat ${POP}_het.hwe | tr ‘/’ ‘\t’ > ${POP}_het.txt
   done
```  
Now move the output files to the main folder FIBR_ms/data/heterozygosity 
and then run: ```Het_windows.R```  

*Allele frequency*  
On the HPC, in the folder with the raw output, run the following:  
```
for POP in GHP GLP IC IT ILL IUL
   do
   cat ${POP}_AF.frq | tr ':' '\t' | sed '1s/^.*$/CHROM\tPOS\tN_ALLELES\tN_CHR\tREF\tREF_FRQ\tALT\tALT_FRQ/' > ${POP}_AF.txt
   done
```  
Now move the output files to the main folder FIBR_ms/data/allele_freq  
and then run: ```AF_windows.R```  to calculate (absolute) delta allele frequencies.  

To calculate fixed/unchanged AFs between GHP and the LP populations run: ```fixed_AFs.R```  

*Plot genome-wide stats*  
To create the plot seen in fig 1D, run  
``` plot_genomewide_stats.R ```  


#### PCA
**To obtain the raw data**  
Run this scripts on the HPC:  
PCA - ```plink_PCA.sh```  

**To plot the PCA**  
Move the eigenvector and eigenvalue files from the HPC to the local FIBR_ms/data/PCA
Then run ```PCA_plot.R``` to create the plot in figure 1B


#### Runs of homozygosity
This analysis requires phased VCFs  

**To obtain the raw data**  
Run these scripts on the HPC to get ROH per individual per population:  
Runs of homozygosity - ```plink_ROH.sh```  

**To process the raw data**  
Move the output files from the HPC output folder to the local FIBR_ms/data/ROH
To process the files and create the plot seen in figure 1C run: ```ROH_plots.R``` 


#### Selscan genome scans
This analysis requires phased VCFs  

**To obtain the raw data**  
Run this scripts on the HPC to get XP-EHH and iHH12:  
Selscan - ```selscan.sh```  

