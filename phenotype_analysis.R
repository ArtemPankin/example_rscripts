#install.packages("pacman")
library(pacman)
p_load(data.table, plyr, dplyr, ggplot2, gplots, corrplot, ggridges, ggplotify, gridGraphics, cowplot, patchwork)

# load input data

data.flowering.2017 <- fread(file="~/mounts/cluster/project/projects/PlantGen/artem/gwas/phenotypes/final.2017_field_mod.tab", header = F)
data.flowering.2013LD <- fread(file="~/mounts/cluster/project/projects/PlantGen/artem/gwas/phenotypes/final.2013_chamber_LD_mod.tab", header = F)
data.flowering.2013SD <- fread(file="~/mounts/cluster/project/projects/PlantGen/artem/gwas/phenotypes/final.2013_chamber_SD_mod.tab", header = F)
data.flowering.nonver <- fread(file="~/mounts/cluster/project/projects/PlantGen/artem/gwas/phenotypes/final.nonvernalized_ipk_mod.tab", header = F)
data.flowering.ver <- fread(file="~/mounts/cluster/project/projects/PlantGen/artem/gwas/phenotypes/final.vernalized_ipk_mod.tab", header = F)
data.flowering.2009SD <- fread(file="~/mounts/cluster/project/projects/PlantGen/artem/gwas/phenotypes/final.2009_shortday_mod.tab", header = F)
data.flowering.2009LD <- fread(file="~/mounts/cluster/project/projects/PlantGen/artem/gwas/phenotypes/final.2009_longday_mod.tab", header = F)
data.flowering.2009LD$V2 <- as.factor(data.flowering.2009LD$V2)
data.flowering.2009SD$V2 <- as.factor(data.flowering.2009SD$V2)

data.flowering.list1 <- list(data.flowering.2017, 
                             data.flowering.2013LD, 
                             data.flowering.2013SD, 
                             data.flowering.nonver, 
                             data.flowering.ver, 
                             data.flowering.2009LD,
                             data.flowering.2009SD)

# adjusted no_flowering values = max + 20

data.flowering.list <- lapply(data.flowering.list1,function(x){
  a <- max(as.numeric(as.vector(unlist(x[-1]))), na.rm = T) + 20
  x <- plyr::colwise(function(y)
    revalue(y,c("no_flowering"=a)))(x)
})

a <- c("MPIPZ_VRN", "CH_LD", "CH_SD", "IPK_NVRN", "IPK_VRN", "GL_LD","GL_SD")

data.flowering.list <- lapply(seq(1:length(data.flowering.list)), function(x){
  cbind(exp= a[x],data.flowering.list[[x]])
})

# converting data to date of sowing 

data.flowering.list[[4]] <- data.flowering.list[[4]] %>%
  mutate_at(colnames(data.flowering.list[[4]])[3:6], function(x)
    as.numeric(as.character(x))-96)

data.flowering.list[[5]] <- data.flowering.list[[5]] %>%
  mutate_at(colnames(data.flowering.list[[5]])[3:6], function(x)
    as.numeric(as.character(x))-75)

data.flowering.list[[6]] <- data.flowering.list[[6]] %>%
  mutate_at(colnames(data.flowering.list[[6]])[3], function(x)
    as.numeric(as.character(x)) + 28)

data.flowering.list[[7]] <- data.flowering.list[[7]] %>%
  mutate_at(colnames(data.flowering.list[[7]])[3], function(x)
    as.numeric(as.character(x)) + 28)

# converting data to numeric values

data.flowering.list <- lapply(data.flowering.list, function(x){
  x %>%
    mutate_at(colnames(x)[3:ncol(x)],function(x)as.numeric(as.character(x)))
})

## converting list to data.frame

data.flowering.all <- do.call(bind_rows,data.flowering.list)

# calculaiing mean flowering time across the replicates and correlations; subsets: wild_domesticated; wild; domesticated

data.flowering.wdn <- data.flowering.all %>% 
  subset(V1 %in% wild_nonadm$V1 | V1 %in% dom_nonadm$V1) %>%
  data.table::melt(id = c("exp","V1")) %>%
  dplyr::mutate(value = as.numeric(value))  %>%
  dplyr::group_by(exp,V1) %>%
  dplyr::summarize(val = round(mean(value, na.rm = T), 0)) %>%
  data.table::dcast(V1 ~ exp)

flowering.cor.wdn <- data.flowering.wdn[,-1] %>%
  apply(.,2, function(y){
    apply(., 2, function(x){
      cor.test(y, x, alternative = "two.sided", method = c("pearson"))$estimate
    })
  })

data.flowering.wn <- data.flowering.all %>% 
  subset(V1 %in% wild_nonadm$V1) %>%
  data.table::melt(id = c("exp","V1")) %>%
  dplyr::mutate(value = as.numeric(value))  %>%
  dplyr::group_by(exp,V1) %>%
  dplyr::summarize(val = round(mean(value, na.rm = T), 0)) %>%
  data.table::dcast(V1 ~ exp)

flowering.cor.wn <- data.flowering.wn[,-1] %>%
  apply(.,2, function(y){
    apply(., 2, function(x){
      cor.test(y, x, alternative = "two.sided", method = c("pearson"))$estimate
    })
  })

data.flowering.dn <- data.flowering.all %>% 
  subset(V1 %in% dom_nonadm$V1) %>%
  data.table::melt(id = c("exp","V1")) %>%
  dplyr::mutate(value = as.numeric(value))  %>%
  dplyr::group_by(exp,V1) %>%
  dplyr::summarize(val = round(mean(value, na.rm = T), 0)) %>%
  data.table::dcast(V1 ~ exp)

flowering.cor.dn <- data.flowering.dn[,-1] %>%
  apply(.,2, function(y){
    apply(., 2, function(x){
      cor.test(y, x, alternative = "two.sided", method = c("pearson"))$estimate
    })
  })

## plotting correlation heatmaps

dev.off()

corrplot(flowering.cor.wn, order = "hclust", type = c("lower"), method = c("color"), 
         number.cex = 1.5, tl.cex = 1.2, addCoef.col = "black", cl.cex = 1.3, mar = c(0, 1.5, 0, 0), diag = F, tl.srt = 20, tl.col = "black", cl.lim = c(0,1), rect.col = "white",
         col= colorRampPalette(c("red","white","#E64B35FF"))(20))
grid.echo()
b <- as.ggplot(grid.grab())

corrplot(flowering.cor.dn, order = "hclust", type = c("lower"), method = c("color"), 
         number.cex = 1.5, tl.cex = 1.2, addCoef.col = "black", cl.cex = 1.3, mar = c(0, 1.5, 0, 0), diag = F, tl.srt = 20, tl.col = "black",cl.lim = c(0,1), rect.col = "white",
         col= colorRampPalette(c("red","white","#4DBBD5FF"))(20))    
grid.echo()
c <- as.ggplot(grid.grab())

## END plotting correlation heatmaps

## plotting phenotype distributions; all environments

data.flowering.wdn_toplot <- data.flowering.wdn %>%
  melt() %>%
  left_join(wild_nonadm) %>%
  mutate(V2 = replace_na(V2,"cultivated")) %>%
  mutate(V4 = paste(variable,V2)) %>% 
  subset(!V4 %in% "MPIPZ_VRN cultivated") %>%
  subset(!grepl("IPK_VRN", V4) | !  value == 134) %>% 
  subset(!grepl("IPK_NVRN", V4) | ! value == 113) %>% 
  subset(!grepl("GL_SD", V4) | ! value == 179) %>%
  subset(!grepl("CH_SD", V4) | ! value == 267) %>%
  subset(!grepl("GL_LD", V4) | ! value == 125)


## count and plot number of non-flowering plants

nonflower <- data.flowering.wdn %>%
  melt() %>% 
  left_join(wild_nonadm)  %>%
  mutate(V2 = replace_na(V2,"cultivated")) %>%
  mutate(V4 = paste(variable,V2)) %>% 
  subset(!V4 %in% "MPIPZ_VRN cultivated") %>%
  dplyr::group_by(V4) %>%
  drop_na(value) %>%
  dplyr::count(name = "total") %>%
  separate(V4, c("variable","V2"), sep = " ") %>%
  join(nonflower) %>%
  dplyr::mutate(val = round(n / total,2)) # %>%

nonflower <- nonflower[3:12,] %>%
  rbind(.,data.frame(variable = c("CH_LD","CH_LD","MPIPZVRN"), V2 = c(rep(c("cultivated","wild"),1), "wild"), total = rep(0,3), n = rep(0,3), num = c(1,2,13), val = rep(0,3)))

## extract latest flowering data point for each environment

last_flowering <- data.flowering.wdn_toplot[,c(2,3)] %>%
  dplyr::group_by(variable) %>%
  distinct() %>%
  dplyr::top_n(1,value)

last_flowering$yval <- seq(1,14,2)

## create y-axis breaks

env_breaks <- c("MPIPZ_VRN wild",unique(data.flowering.wdn_toplot[grepl("cultivated",data.flowering.wdn_toplot$V4),5]))

## test significance difference between the wild and cult barley in different environments; wilcoxon

wilc_pval <- lapply(unique(constel_sub$b)[-6], function(x){
  wilc <- wilcox.test(subset(data.flowering.wdn_toplot, variable %in% x & V2 %in% "wild")$value, 
                      subset(data.flowering.wdn_toplot, variable %in% x & V2 %in% "cultivated")$value,
                      alternative = "t")
  data.frame(x, wilc$p.value)
}) %>%
  do.call(rbind.data.frame,.)
wilc_pval$y <- c(1,9,7,3,11,5)
wilc_pval <- subset(wilc_pval, wilc.p.value < 0.01)

## plot flowering time ridges with ggridge package

a <- ggplot(data.flowering.wdn_toplot) +
  geom_density_ridges(aes(x=value, y= V4, fill= V2), alpha = .9, jittered_points = T, quantile_lines = T, quantiles = 2,
                      position = position_points_jitter(width = 0, height = 0),
                      point_shape = '|', point_size = 3, point_alpha = .8, colour = "white", point_colour = "black") +
  scale_fill_cyclical(values = c("#4DBBD5FF","#E64B35FF")) +
  geom_point(data = last_flowering, aes(x=value, y = yval), colour = "black", shape = 18, size = 4) +
  geom_text(data = nonflower, aes(x=280, y = num+0.3, label = val), size = 7) +
  geom_text(data = wilc_pval , aes(x=8, y = y+0.7, label = "***"), size = 6) +
  coord_cartesian(clip = "off") +
  scale_y_discrete(breaks = env_breaks, labels = sub(" cultivated| wild","",env_breaks), expand = c(0, 0)) +
  xlab("Days") +
  scale_x_continuous(breaks = seq(0,250,25), expand = c(0, 10)) +
  theme_ridges(grid = TRUE) +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(hjust = .5, vjust = -2), 
        axis.text.y = element_text(vjust = -2, size = 14), 
        axis.text.x = element_text(vjust = -2, size = 14), 
        #panel.grid.major.x =  element_line(colour="white", size=.7),
        panel.grid.major.y =  element_line(colour="black", size=.7))

## END plotting phenotype distributions; all environments

## saving publication ready figure

ggsave("~/mounts/cluster/project/projects/PlantGen/artem/gwas/figures/main_phenotypes.pdf",
       (a | (b / c)),
       dpi = "retina",
       width = 41,
       height = 26, 
       unit = "cm")
## evaluate plasticity and gxe interactions


# estimating plasticity as CV between traits of each genotype by subspecies

g <- data.flowering.wdn %>%
  melt() %>%
  left_join(wild_nonadm) %>%
  mutate(V2 = replace_na(V2,"cultivated")) %>%
  mutate(V4 = paste(variable,V2)) %>% 
  subset(!V4 %in% "MPIPZ_VRN cultivated") %>%
  dplyr::group_by(V2, V1) %>%
  dplyr::summarize(cv = sd(value, na.rm = T) / mean(value, na.rm = T)*100) 

e <- ggplot(g,aes(x=V2,y=cv, fill = V2)) +
  geom_flat_violin(position = position_nudge(x = .12, y = 0), alpha = 0.7) +
  geom_jitter(width = .02, alpha = .5) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  scale_fill_manual(values = c("#4DBBD5FF","#E64B35FF")) +
  labs(y = "CV, %") +
  theme_bw() +
  raincloud_theme +
  geom_signif(comparisons = list(c("cultivated","wild")), 
              map_signif_level=TRUE, 
              test = "wilcox.test", margin_top = 0.07, textsize = 5)

ggsave("~/mounts/cluster/project/projects/PlantGen/artem/gwas/figures/suppl_cv_plastisity.pdf",
       e,
       dpi = "retina",
       width = 17,
       height = 13, 
       unit = "cm")


## test attempt to plot gxe interactions as radar plot - not very readable

h <- data.flowering.dn %>%
  as_tibble() %>%
  mutate_at(vars(-V1), rescale) %>%
  #head(10) %>%
  #select(c(1,4,6,7,2)) %>%
  drop_na() %>%
  ggradar(legend.position = "", group.point.size = 1, group.line.width = .2, axis.label.size = 5, axis.line.colour = "white", values.radar = NA,
          label.gridline.min = F, label.gridline.mid = F,
          label.gridline.max = F,
          background.circle.transparency = 0,
          gridline.mid.colour = NA)



r <- data.flowering.wn %>%
  as_tibble() %>%
  mutate_at(vars(-V1), rescale) %>%
  #head(10) %>%
  #select(c(1,4,6,7,2)) %>%
  drop_na() %>%
  ggradar(legend.position = "", group.point.size = 1, group.line.width = .2, axis.label.size = 5, axis.line.colour = "white", values.radar = NA,
          label.gridline.min = F, label.gridline.mid = F,
          label.gridline.max = F,
          background.circle.transparency = 0,
          gridline.mid.colour = NA)

plot_grid(h,r, nrow = 1, labels=letters[1:2], label_size = 18) 


## table - descriptive statistics of flowering distributions

data.flowering.wdn %>%
  melt() %>%
  left_join(wild_nonadm) %>%
  mutate(V2 = replace_na(V2,"cultivated")) %>%
  mutate(V4 = paste(variable,V2)) %>% 
  subset(!V4 %in% "MPIPZ_VRN cultivated") %>%
  dplyr::group_by(V2, variable) %>%
  drop_na(value) %>%
  dplyr::summarize(number = n(), 
                   mean = round(mean(value, na.rm = T),2), 
                   median = round(median(value, na.rm = T),2), 
                   sd = round(sd(value, na.rm = T),2), 
                   cv = round(sd(value, na.rm = T) / mean(value, na.rm = T)*100,2)) %>%
  drop_na(cv) %>%
  write.table("~/Downloads/temp.xls")

## anova to quantify GxE

data.flowering.wdn.anova <- data.flowering.all %>% 
  subset(V1 %in% dom_nonadm$V1 | V1 %in% wild_nonadm$V1) %>%
  data.table::melt(id = c("exp","V1")) %>%
  dplyr::mutate(value = as.numeric(value)) %>%
  left_join(wild_nonadm) %>%
  mutate(V2 = replace_na(V2,"cultivated")) %>%
  setNames(c("envir","genot","repl","value","type"))


combenv <-combn(colnames(data.flowering.dn)[-1], 2)
crossover_dom <- lapply(1:ncol(combenv), function(z){
  j= 0
  a <- combenv[,z]
  b <- data.flowering.dn[,c("V1",a)] %>%
    drop_na()
  combgenot <- combn(b$V1, 2)
  
  lapply(1:ncol(combgenot), function(p){
    
    d <- b %>%
      subset(V1 %in% combgenot[,p])
    x <- d[1,2] - d[2,2]
    y <- d[1,3] - d[2,3]
    
    if(sign(x) != sign(y) & abs(x) > 4 & abs(y) > 4){
      j <<- j+1
      
    }
})
  c(a, j, j / (nrow(b)^2 - nrow(b))/2)
})

crossover_dom.dt <- do.call(rbind.data.frame, crossover_dom) %>%
  setNames(c("env1","env2","num","value"))

crossover_wild.dt <- do.call(rbind.data.frame, crossover_wild) %>%
  setNames(c("env1","env2","num","value"))

crossover_dt <- inner_join(crossover_dom.dt, crossover_wild.dt, by = c("env1","env2")) %>%
  unite("env_pair",c("env1","env2")) 

crossover_dt_melt <- crossover_dt[,c(1,3,5)] %>%
  setNames(c("env_pair","dom","wild")) %>%
  data.table::melt(id.vars = "env_pair")

ggplot(crossover_dt_melt) +
  geom_point(aes(x=variable, y = as.numeric(value), colour = env_pair)) +
  geom_line(aes(x=variable, y = as.numeric(value), group = env_pair))

crossover_dt$diff <- as.numeric(as.character(crossover_dt$value.x))-as.numeric(as.character(crossover_dt$value.y))

## gxe plot and significance - pairwise environments

combenv1 <- combenv[,-which(apply(combenv,2,function(y) any(grepl("GL_",y, fixed=TRUE))))]

gxe_plots <- vector(mode = "list", length = 6)

lapply(1:6, function(x){
  
  cenv <- combenv1[,x]
  
  fl_subs <- data.flowering.wdn_toplot %>%
    subset(variable %in% cenv)
  
  fl_incl <- fl_subs %>% 
    dplyr::group_by(V1) %>%
    drop_na(value) %>%
    dplyr::summarise(n()) %>%
    .[.[,2] > 1, ] %>%
    .[,1] %>%
    ungroup()
  
  fl_subs <- subset(fl_subs, V1 %in% fl_incl$V1)
  
  fl_subs_w <- subset(fl_subs, V2 %in% "wild")
  fl_subs_c <- subset(fl_subs, V2 %in% "cultivated")
  
  fl_subs_w$V1_ord <- factor(fl_subs_w$V1, 
                             levels = subset(fl_subs_w, variable %in% cenv[1])$V1[rev(order(fl_subs_w$V2, fl_subs_w$value))])
  fl_subs_c$V1_ord <- factor(fl_subs_c$V1, 
                             levels = subset(fl_subs_c, variable %in% cenv[1])$V1[rev(order(fl_subs_c$V2, fl_subs_c$value))])
  
  ## calculate gxe significance by group for pairwise environments
  
  fl_subs_w_gxe <- data.flowering.wdn.anova %>%
    subset(genot %in% fl_subs_w$V1 & envir %in% cenv) %>%
    aov(value ~ genot*envir, data = .) %>%
    drop1(.,~.,test="F")
  
  fl_subs_c_gxe <- data.flowering.wdn.anova %>%
    subset(genot %in% fl_subs_c$V1 & envir %in% cenv) %>%
    aov(value ~ genot*envir, data = .) %>%
    drop1(.,~.,test="F")
  
  pvals_w <- fl_subs_w_gxe[["Pr(>F)"]]
  sigSymbols_w <- symnum(pvals_w, na = FALSE, 
                         cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                         symbols = c("***", "**", "*", ".", " "))[4]
  
  pvals_c <- fl_subs_w_gxe[["Pr(>F)"]]
  sigSymbols_c <- symnum(pvals_c, na = FALSE, 
                         cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                         symbols = c("***", "**", "*", ".", " "))[4]
  
  # PLOT GxE pairwise environments
  
gxe_plots[[x]] <<- ggplot() +
    geom_segment(data = subset(fl_subs_w, variable %in% cenv[1]), aes(x=0, xend = -value, y = as.numeric(V1_ord), yend = as.numeric(V1_ord)), colour = "#4DBBD5FF") +
    geom_segment(data = subset(fl_subs_w, variable %in% cenv[2]), aes(x=0, xend = value, y = as.numeric(V1_ord), yend = as.numeric(V1_ord)), colour = "#4DBBD5FF") +
    geom_segment(data = subset(fl_subs_c, variable %in% cenv[1]), aes(x=0, xend = -value, y = -as.numeric(V1_ord), yend = -as.numeric(V1_ord)), colour = "#E64B35FF") +
    geom_segment(data = subset(fl_subs_c, variable %in% cenv[2]), aes(x=0, xend = value, y = -as.numeric(V1_ord), yend = -as.numeric(V1_ord)), colour = "#E64B35FF") +
    geom_vline(aes(xintercept=0), alpha = .5) +
    geom_label(aes(x= -5, y = -nrow(fl_subs_c) / 2 - 35, label = cenv[1]), size = 2.5, hjust = "right") +
    geom_label(aes(x= 5, y = -nrow(fl_subs_c) / 2 - 35, label = cenv[2]), size = 2.5, hjust = "left") +
    geom_text(aes(x= -100, y = nrow(fl_subs_w) / 2 + 15, label = sigSymbols_w)) +
    geom_text(aes(x= -100, y = -nrow(fl_subs_c) / 2 - 15, label = sigSymbols_c)) +
    scale_x_continuous(labels=abs) +
    theme_bw() %+replace% theme(axis.title = element_blank(), 
                                axis.title.x = element_blank(), 
                                axis.text.x = element_text(), 
                                axis.line.y = element_blank(), 
                                axis.text.y = element_blank(), 
                                axis.ticks.y = element_blank(), 
                                legend.position = "none", 
                                panel.border = element_blank(), 
                                axis.line = element_line(), 
                                panel.grid = element_blank())
  
})

ggsave("~/mounts/cluster/project/projects/PlantGen/artem/gwas/figures/gxe_pairwise_envir.pdf",
       plot_grid(plotlist = gxe_plots, ncol = 3),
       dpi = "retina",
       width = 13,
       height = 18, 
       unit = "cm")

## END gxe plot - pairwise

## quantifying GxE and H2 heritability

flow.for.lmer_wn <-  
  data.flowering.wdn.anova %>%
  subset(!grepl("GL_",envir)) %>%
  subset(type %in% "wild") %>%
  drop_na(value) %>%
  subset(!repl %in% "V5")

flow.for.lmer_dn <-  
  data.flowering.wdn.anova %>%
  subset(!grepl("GL_",envir)) %>%
  subset(type %in% "cultivated") %>%
  drop_na(value) %>%
  subset(!repl %in% "V5")

gxe_herit_pairwise <- lapply(1:6, function(x){  
  
  cenv <- combenv1[,x]  
  
  wdata <- flow.for.lmer_wn %>%
    subset(envir %in% cenv)
  ddata <- flow.for.lmer_dn %>%
    subset(envir %in% cenv)
  
  my_model_wn <- lmerTest::lmer(value ~ (1|genot) + (1|envir) + (1|genot:envir) + (1|repl), wdata)
  my_model_dn <- lmerTest::lmer(value ~ (1|genot) + (1|envir) + (1|genot:envir) + (1|repl) , ddata)
  
  h <- print(VarCorr(my_model_dn),comp="Variance")
  ged <- as.data.frame(h)$vcov[1] / sum(as.data.frame(h)$vcov[-4])
  h2d <- as.data.frame(h)$vcov[2] / sum(as.data.frame(h)$vcov[-4])
  
  q <- print(VarCorr(my_model_wn),comp="Variance")
  gew <- as.data.frame(q)$vcov[1] / sum(as.data.frame(q)$vcov[-4])
  h2w <- as.data.frame(q)$vcov[2] / sum(as.data.frame(q)$vcov[-4])
  
  data.frame(env = paste0(cenv[1],"_",cenv[2]), round(gew,2),round(ged,2),round(h2w,2),round(h2d,2))
  
}) %>%
  do.call(rbind.data.frame, .) %>%
  write.table("~/Downloads/temp.txt")

## END quantifying GxE and H2 heritability