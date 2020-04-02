## This script performs GWA using mlmm model (Segura et al., Nat Genet) on a panel of 
## exome-enrichment homozygous SNP, INDEL, and presence-absence VRN1 intron locus markers.
## It extracts the BIC outliers from gwas and all the tightly linked markers (haplotypes). 
## The output of the script is the publication-ready table of thoroughtly annotated GWA signals. 

setwd("~/mounts/cluster/project/projects/PlantGen/artem/gwas/")

#install.packages("pacman")
library(pacman)
p_load(mlmm, corpcor, stringr, grid, UpSetR, dplyr, tidyr, data.table, purrr)

##create gwas panel file from plink files

constel <- data.frame(a=sort(rep(c("wn", "wdn", "all"), 7)),
                      b = rep(c("_2013_chamber_ld_mod", 
                                "nonvernalized_ipk_mod",  
                                "x_2009_shortday", 
                                "_2013_chamber_sd_mod", 
                                "vernalized_ipk_mod", 
                                "_2017_field", 
                                "x_2009_longday"), 3)
)

constel_sub <- constel
constel_sub$b <- gsub("_2013_chamber_ld_mod", "CH_LD", constel$b)  %>%
{ gsub("nonvernalized_ipk_mod", "IPK_NVRN", .) } %>%
{ gsub("x_2009_shortday", "GL_SD", .) } %>%
{ gsub("_2013_chamber_sd_mod", "CH_SD", .) } %>%
{ gsub("vernalized_ipk_mod", "IPK_VRN", .) } %>%
{ gsub("_2017_field", "MPIPZ_VRN", .) } %>%
{ gsub("x_2009_longday", "GL_LD", .) } 
constel_sub$c <- paste0(constel_sub$a,"_",constel_sub$b)

## extract genotypes for each gwas (different sets of phenotypes; not comparable)

lapply(1:21,function(x){
  set <- constel[x,1]
  env <- constel[x,2]
  genot <- read.table(paste0("version3_mlmm_R/mni_v3_maf_trans_",set,".tab"), sep = "\t", header = T)
  phenot <- read.table(paste0("version3_mlmm_R/",env, ".tab"), sep = "\t", header = T)
  colnames(phenot) <- c("V1","V2")
  phenot_sub <- subset(phenot, V1 %in% genot$genotype)
  phenot_no_NA <- subset(phenot_sub, V2 != "NA")
  write.table(phenot_no_NA$V1, 
              file = paste0("~/mounts/cluster/project/projects/PlantGen/artem/gwas/phenotyped_genotypes_",set,env,".list"),
              sep = "\t", 
              col.names = T, 
              row.names = F, 
              quote = F)
  
  print(x)
})

## BASH commands (implement bash execution from R script): creates ped files for each gwas panel; lists of markers to include (filered by maf > 0.03, hets removed, missing < 0.05, r2 < 1, removed 11-12 collapsed markers) 

# for file in phenotyped_genotypes_*list; do plink1.9 -file wild_dom_409.maf.003.miss0.1.nofilt_ld_1.2_miss0.1.no_hets1 --recode --geno 0.05 --keep <(paste $file $file) --not-chr 9,10 --indep-pairwise 25000000 5000000 0.999 --memory 1000000 --threads 2 --out $file.no910.5miss       ; done && rm *.nosex *.log

## genotype IDs for each set of accessions

all_genot <- read.table("all_ids.1col", header = F)
wdn_genot <- read.table("wdn_ids.1col", header = F)
wn_genot <- read.table("wn_ids.1col", header = F)

## import genotypes (prefiltered by maf (for the complete set) and missing data, removed heterozyg)

genot <- read.table("wild_dom_409.maf.003.miss0.1.nofilt_ld_1.2_miss0.1.no_hets1.gwas.out", sep = "\t", header = F)
snps <- read.table("wild_dom_409.maf.003.miss0.1.nofilt_ld_1.2_miss0.1.no_hets1.map")[,2]
colnames(genot) <- c("genotype",as.character(snps))

## converting 0,1,2 to NA minor allele -0 and major allele -2 (removing heterozygous genotypes)

genot[genot == 0] <- NA

genot_mod <- apply(genot[,-1], 2, function(x){
  a <- sum(x == 1, na.rm = TRUE) / sum(!is.na(x)) 
  if(a >= 0.5){
    gsub("1","0",x) %>% as.numeric()
  } else {
    gsub("2","0",x) %>%
      gsub("1","2",.) %>% as.numeric()
  }
})

colnames(genot_mod) <- as.character(snps)
rownames(genot_mod) <- genot[,1]

## run mlmm for all environments gwa; maf 0.04, missing data .01

mygwas_list16 <- lapply(1:21,function(x){
set <- constel[x,1]
env <- constel[x,2]
set_gen <- paste0(set,"_genot")

prune_in <- read.table(paste0("phenotyped_genotypes_",set,env,".list.no910.5miss.prune.in")) ## for gwas 16
prune_in <- as.character(prune_in$V1)

genot_gwa <- genot_mod[as.character(get(set_gen)$V1),prune_in]
phenot <- read.table(paste0("version3_mlmm_R/",env, ".tab"), sep = "\t", header = T)
colnames(phenot) <- c("id","val")
phenot_sub <- subset(phenot, id %in% rownames(genot_gwa))
phenot_no_NA <- subset(phenot_sub, val != "NA")
genot_no_NA <- genot_gwa[as.character(phenot_no_NA$id),]

# loading PCA results; to correct for pop structure as co-factors

pcacof <- read.table(paste0("version3_mlmm_R/",set,".genotypes.v3.pca"))
pcacof$V1 <- gsub("-","_",pcacof$V1)
pcacof_no_NA <- subset(pcacof, V1 %in% phenot_no_NA$id)
pcacof_no_NA_mat <- as.matrix(pcacof_no_NA[,2:ncol(pcacof_no_NA)])

## filtering by maf

b <- apply(genot_no_NA, 2, function(y){
  sum(str_count(y, "2"), na.rm = T) / sum(!is.na(y))
})
c <- b[ b < 0.03] # maf filtering on

#c <- b[b < 0 | b > 1] # maf filtering off
if(length(c) > 0){
  genot_no_NA_maf <- genot_no_NA[,-which(colnames(genot_no_NA) %in% names(c))]
} else
{
  genot_no_NA_maf <- genot_no_NA
}

## imputating missing genotypes

genot_no_NA_mat <- as.matrix(genot_no_NA_maf)

average <- colMeans(genot_no_NA_mat, na.rm = T)

invisible(lapply(1:ncol(genot_no_NA_mat),function(i){
  genot_no_NA_mat[is.na(genot_no_NA_mat[,i]), i] <<- average[i]
}))

##### calculating kinship matrix #####

average <- colMeans(genot_no_NA_mat, na.rm = T)
stdev <- apply(genot_no_NA_mat, 2, sd)
genot_stand <- sweep(sweep(genot_no_NA_mat, 2, average, "-"), 2, stdev, "/")
K_mat <- (genot_stand %*% t(genot_stand)) / ncol(genot_stand)
K_mat2 <- make.positive.definite(K_mat)

##### MLMM analysis #####

mygwas <- mlmm_cof(Y = phenot_no_NA$val, X = genot_no_NA_mat, cofs = pcacof_no_NA_mat, K = K_mat2, nbchunks = 2, maxsteps = 500)
mygwas$maf <- b

return(mygwas)

})

## add q values
mygwas_list16_pvals <- lapply(mygwas_list16, function(x){
  fdr <- p.adjust(x$opt_extBIC$out$pval, method = "fdr")
  cbind(x$opt_extBIC$out, fdr)
})

## filling in with maf values

mygwas_list16_pvals <- lapply(1:21, function(x){
  a <- as.data.frame(mygwas_list16[[x]]$maf)
  a$SNP <- rownames(a)
  names(a)[1] <- "maf"
  inner_join(mygwas_list16_pvals[[x]], a)
})

## naming list elements as subset_environment

names(mygwas_list16_pvals) <- constel_sub$c

## adding envinronment values

mygwas_list16_pvals <- lapply(1:21, function(x){
  mygwas_list16_pvals[[x]]$env <- names(mygwas_list16_pvals)[x]
  return(mygwas_list16_pvals[[x]])
})

## outliers by BIC

outliers_mlmm_BIC_list16 <- lapply(1:21, function(x){
  subset(mygwas_list16_pvals[[x]], SNP %in% mygwas_list16[[x]]$opt_extBIC$cof)
})

## outliers by fdr

outliers_mlmm_fdr_list16 <- lapply(1:21, function(x){
  subset(mygwas_list16_pvals[[x]], fdr < 0.05)
})

names(outliers_mlmm_fdr_list16) <- names(mygwas_list16_pvals)
names(outliers_mlmm_BIC_list16) <- names(mygwas_list16_pvals)

## exrtacting overlaps between FDR and BIC outliers

lapply(1:21, function(x){
  a <- length(intersect(outliers_mlmm_fdr_list16[[x]]$SNP, outliers_mlmm_BIC_list16[[x]]$SNP))
  b <- length(outliers_mlmm_fdr_list16[[x]]$SNP)
  c <- length(outliers_mlmm_BIC_list16[[x]]$SNP)
  c(a,b,c)
})

## separating BIC outlier marker ID into chr and position

outliers_mlmm_BIC_list16 <- lapply(1:21,function(x){
  separate(outliers_mlmm_BIC_list16[[x]],"SNP", c("chr","pos"),sep = ":",remove = F)
})

## separating ID into chr and position for all markers

allmarkers_mlmm_BIC_list16 <- lapply(1:21,function(x){
  separate(mygwas_list16[[x]]$opt_extBIC$out,"SNP", c("chr","pos"),sep = ":",remove = F)
})

names(allmarkers_mlmm_BIC_list16) <- constel_sub$c

## adding estimate + std. error 

outliers_mlmm_BIC_list16 <- lapply(1:21, function(x){
  a <- as.data.frame(mygwas_list16[[x]]$opt_extBIC$coef)
  colnames(a)[4] <- "pval"
  left_join(outliers_mlmm_BIC_list16[[x]],a)[,-10]
})

names(outliers_mlmm_BIC_list16) <- constel_sub$c

# removing environments w/o outliers

outliers_mlmm_BIC_list16 <- outliers_mlmm_BIC_list16 %>%
  keep(function(x) nrow(x) > 0)

## exporting .annot files for clumping mlmm data into haplotypes of linked markers

lapply(1:length(allmarkers_mlmm_BIC_list16), function(x){
  c <- data.frame(CHR=allmarkers_mlmm_BIC_list16[[x]]$chr, SNP=allmarkers_mlmm_BIC_list16[[x]]$SNP, BP=allmarkers_mlmm_BIC_list16[[x]]$pos, P=allmarkers_mlmm_BIC_list16[[x]]$pval) 
  write.table(c, 
              file = paste0("~/mounts/cluster/project/projects/PlantGen/artem/gwas/list16.mlmm_",names(allmarkers_mlmm_BIC_list16)[x],".assoc"),
              sep = "\t", 
              col.names = T, 
              row.names = F, 
              quote = F)
})

## extracting p-value thresholds for clumping - different in each gwa

lapply(1:length(outliers_mlmm_BIC_list16), function(x){
  d <- max(outliers_mlmm_BIC_list16[[x]]$pval) + 0.0001* max(outliers_mlmm_BIC_list16[[x]]$pval)
  write.table(d, file = paste0("~/mounts/cluster/project/projects/PlantGen/artem/gwas/list16.mlmm_",names(outliers_mlmm_BIC_list16)[x],".pval"),
              sep = "\t", 
              col.names = T, 
              row.names = F, 
              quote = F)
})

## BASH commands: run gwas_clumping_by_ld gwas_annot.sh

## BASH commands: run these two commands for clumping and parsing clumped files

# for file in list16.mlmm_*.assoc; do pval=`cat ${file%%assoc}pval | sed -n 2p`; con=`echo ${file%%.assoc} | sed 's/list16.mlmm_//'`; con1=`grep $con constel_sub | awk '{print $1}'`; plink1.9 --file phenotyped_genotypes_$con1.list.no910.5miss --clump $file --clump-p1 $pval --clump-p2 1 --clump-r2 0.75 --clump-kb 100000000 --out $file && rm $file.nosex; done
# rm list16.interchr_ld1 && for file in list16.mlmm_*.assoc.clumped; do con=`echo ${file%%.assoc.clumped} | sed 's/list16.mlmm_//'`; con1=`grep $con constel_sub | awk '{print $1}'`; cat $file | awk '{print $3}' | grep -v SNP | awk /./ > temp1 && plink1.9 --file $con1.list.no910.5miss  --r2 inter-chr with-freqs yes-really --ld-snp-list temp1 --ld-window-r2 0.75  --out temp.$file && cat temp.$file.ld | awk '$1 != $5' | grep -v CHR_A | awk '{print "'$con'\t'$con1'\t"$0}' >> list16.interchr_ld1; done ## extracts intrachromosomal LD for the associations
# ls list16.mlmm*.clumped | grep -v p1 | while read file; do cat $file |  awk /./ |  grep -v CHR |  while read line; do orig=`echo $line | grep -v CHR | awk '{print $3}'`;  env=`echo $file`;  echo $line |  awk '{print $3"\n"$NF}' |  sed 's/([1-9])//g' |  tr "," "\n"  |  awk '{print "'${env}'\t'${orig}'\t"$0}' |  sed 's/.assoc.clumped//' |  awk '{if($2 == $3)print $0"\toriginal";else print $0"\tclumped"}'; done; done |  grep -v "NONE" | awk '{print $0"\t"$3}' |  awk 'BEGIN{FS=OFS="\t"}{gsub(":","\t",$NF);print}' | sed 's/list16.mlmm_//' | cat - <(cat list16.interchr_ld1 | awk '{print $1"\t"$5"\t"$9"\tinterchr\t"$7"\t"$8}') | sort -k1,1 -k2,2 > list16.clumped_inter_all
# for file in  temp.list16.mlmm_*.assoc.clumped.ld; do env=`echo $file | sed 's/temp.list16.mlmm_//' | sed 's/.assoc.clumped.ld//'`; cat $file | grep -v "CHR_A" | awk 'BEGIN{FS=OFS="\t"}{print "'$env'\t"$0}'; done | tr " " "\t"  > list16.outliers.ld.all ## extracts LD between outliers + clumped and original

## loading clumped outlier haplotypes

outliers_clumped16 <- fread("~/mounts/cluster/project/projects/PlantGen/artem/gwas/list16.clumped_inter_all")

outliers_clumped_list16 <- split(outliers_clumped16, f = outliers_clumped16$V1)

# sort outliers_clumped by names in annot_outliers to 

outliers_clumped_list16 <- outliers_clumped_list16[names(outliers_mlmm_BIC_list16)]

## ANNOTATING gwa outliers

# loading files for annotation

# allele frequency

majmin <- fread("~/mounts/cluster/project/projects/PlantGen/artem/gwas/minor_major_alleles", header = F) ## output of gwas_annot.sh
names(majmin) <- c("SNP","minor_allele","major_allele","type")

# r2 values between the linked makers within haplotypes

outliers_ld <- read.table("~/mounts/cluster/project/projects/PlantGen/artem/gwas/list16.outliers.ld.all", header = F)  

# homology with arabidopsis

horvu_to_at <- fread("~/mounts/cluster/project/projects/PlantGen/reference/blastx_HighLow_vs_TAIR10.corrected.out", header = F)
horvu_to_at <- horvu_to_at[,c(1,3,16)]
names(horvu_to_at) <- c("HORVU_locus","TAIR_locus","homolog/ortholog")

# arabidopsis annotations

tair_annot <- fread("~/mounts/cluster/project/projects/PlantGen/artem/barley_genome_genes/TAIR_gene_aliases_20130831.txt", header = F)
names(tair_annot)[1] <- "TAIR_locus"

# snpeff annotations; functional characteristics

snpeffind <- fread("~/mounts/cluster/project/projects/PlantGen/artem/important_files_misc/snpeff.indels.table.main_modified", fill = T, header = F)
snpeffsnp <- fread("~/mounts/cluster/project/projects/PlantGen/artem/important_files_misc/snpeff.snps.table.main_modified", fill = T, header = F)
snpeffft1 <- fread("~/mounts/cluster/project/projects/PlantGen/artem/important_files_misc/snpeff.ft1.table.main_modified", fill = F, sep = "\t", header = F)

snpeff <- rbind(snpeffsnp, snpeffind, snpeffft1)
snpeff$SNP <- unite(snpeff[,c(1,2)], "snp" ,sep = ":")
snpeff <- snpeff[,3:5]
names(snpeff)[1:2] <- c("Effect","HORVU_transcript")
snpeff <- subset(snpeff, !Effect %in% "INTERGENIC")

# Pankin et al. New Phyt annotations

horvuseq <- fread("~/mounts/cluster/project/projects/PlantGen/artem/important_files_misc/horvu_to_seq.new", header = F)
names(horvuseq) <- c("HORVU_locus","NewPhyt_seq")

seq_to_gene <- fread("~/mounts/cluster/project/projects/PlantGen/artem/important_files_misc/selection_genes", header = F)
seq_to_gene <- seq_to_gene[,c(1,3,10)]
names(seq_to_gene) <- c("NewPhyt_seq", "selection", "barley_gene")

# tailor-made barley annotations: literature

barley_annot <- fread("~/mounts/cluster/project/projects/PlantGen/reference/160517_Hv_IBSC_PGSB_r1_proteins_HighLowConf_REPR_annotation.annotonly", header = F)
names(barley_annot)[1] <- "HORVU_locus"

## function to merge all possible annotations into data.frame: publication-ready

outliers_annotated <- lapply(1:length(outliers_clumped_list16),function(x){ ###!!!
  
  names(outliers_clumped_list16[[x]])[3] <- "SNP"
  
  a <- outliers_clumped_list16[[x]] %>%
    left_join(outliers_mlmm_BIC_list16[[x]], by = "SNP")
  b <- a[,c(1:6,9:10,13:ncol(a))]
  
  ## extracting distance of clumped to original
  
  c <- split(b, f = b$V2)
  
  q <- sapply(c, function(x){
    y <- subset(x, ! V4 %in% "interchr")
    if (nrow(y) == 1) {
      z <- 0
    } else {
      u <- subset(y, V4 %in% "original")[1,6]
      z <- as.vector(
        apply(x,1,function(b){
          if(!b[4] %in% "interchr"){
            o <- as.numeric(b[6]) - u}
          else { o <- NA}
          return(abs(o))
        }))
    }
    return(z)
  })
  
  ## END extracting distance of clumped to original
  
  b$distance <- as.vector(simplify(q))
  
  ## adding ld between original and clumped markers
  
  ld <- as.vector(apply(b,1, function(x){
    subset(outliers_ld, (V1 %in% x[1] & V4 %in% x[2] & V8 %in% x[3]) | (V1 %in% x[1] & V8 %in% x[2] & V4 %in% x[3]))$V10[1]
  }))
  
  b$ld <-  ld
  
  ## adding major minor alleles 
  
  b <- b %>%
    left_join(majmin, by = "SNP")
  
  tryCatch(b[b$V5 == 11,]$minor_allele <- "1", error = function(e) NULL)
  tryCatch(b[b$V5 == 11,]$major_allele <- "0", error = function(e) NULL)
  tryCatch(b[b$V5 == 11,]$type <- "coverage", error = function(e) NULL)
  
  ## adding MAF
  
  maf <- outliers_ld[,c(1,8:9)] %>%
    setNames(names(outliers_ld[,c(1,4:5)])) %>%
    rbind(outliers_ld[,c(1,4:5)], .) %>%
    unique()
  
  names(maf) <- c("env","SNP","MAF") 
  
  names(b)[1] <- "env"
  
  b1 <- b %>% 
    left_join(., maf, by = c("env", "SNP"))
  
  ## fill in with snpeff annotations
  
  b1 <- b1 %>%
    left_join(snpeff, by = "SNP")
  
  tryCatch(b1[b1$V5 == 11,]$HORVU_transcript <- "HORVU5Hr1G095630.3", error = function(x)NULL)
  tryCatch(b1[b1$V5 == 12,]$HORVU_transcript <- "HORVU7Hr1G024610.1", error = function(x)NULL)
  
  # filling in locus name
  
  b1$HORVU_locus <- gsub("\\..*","",b1$HORVU_transcript)
  
  # filling in seq
  
  b1 <- b1 %>%
    left_join(horvuseq, by = "HORVU_locus")
  
  ## filling in seq to barley gene names
  
  b1 <- b1 %>%
    left_join(seq_to_gene)
  
  ## filling in horvu to tair ids 
  
  b1 <- b1 %>%
    left_join(horvu_to_at, by = "HORVU_locus")
  
  b1$TAIR_locus <- gsub("\\..*","",b1$TAIR_locus)
  
  # filling in TAIR annotations
  
  b1 <- b1 %>% 
    left_join(tair_annot, by = "TAIR_locus")
  
  # filling in barley prot annotations
  
  b1[b1==""] <- "-"
  b1[is.na(b1)] <- "-"
  
  c1 <- b1 %>%
    left_join(barley_annot, by = "HORVU_locus")
  
  c1[c1==""] <- "-"
  c1[is.na(c1)] <- "-"
  c1[c1=="none"] <- "-"
  return(c1)
  
})

outliers_annotated_df <- do.call(rbind.data.frame, outliers_annotated)

outliers_annotated_df <- separate(outliers_annotated_df, env,  c("set","environment"),extra = "merge")

write.table(outliers_annotated_df, file = "~/mounts/cluster/project/projects/PlantGen/artem/gwas/list16_outliers.annotated.txt", quote = F, col.names = F, row.names = F, sep = "\t")

### END