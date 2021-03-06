# install.packages('pheatmap')

library(pheatmap)
library(plyr)
library(AUCRF)
library(ggplot2)
library(dplyr)
library(tidyr)

# Load integrals, peak ranges, and sample data# Load  
integrals = read.table('integrals_20170323.txt', header=FALSE, sep = ",")
integral_ranges = read.table('coculture_peak_ranges.txt',header=TRUE, sep = ",")
master <- read.table(file="merged_metadata.txt",header=TRUE,sep='\t')

# set column names of integrals using integral_ranges met names
colnames(integrals) <- integral_ranges$met_name
mets = colnames(integrals)
# Set rownames using experiment_id. This should be ordered correctly.
integrals$experiment_id <- master$experiment_id
# merge so that classes and integrals are in one dataframe
all_data <- join(master,integrals,by="experiment_id")

# remove non-growing samples from the data
# 500 and 492 did not grow in experiment 3
# 502 did not grow in experiment 2 or 3
all_data = all_data[!(all_data$Run == 2 & grepl('502',all_data$species)),]
all_data = all_data[!(all_data$Run == 3 & grepl('502',all_data$species)),]
all_data = all_data[!(all_data$Run == 3 & grepl('500',all_data$species)),]
all_data = all_data[!(all_data$Run == 3 & grepl('492',all_data$species)),]
# remove experiment 4 (3-species subset including only 356,361,519)
all_data = all_data[all_data$Run != "4",]

# replace negative values in integrals with 0.
all_data[,colnames(all_data) %in% mets] = 
  (apply(all_data[,colnames(all_data) %in% mets],2,function(x){x[x<0] = 0;x}))
# center the data by subtracting the mean of blanks
blanks = subset(all_data,species == '0')
blank_means = colMeans(blanks[,names(blanks) %in% mets])
centered_integrals <- sweep(all_data[,mets], 2, blank_means, `-`)
all_data[,mets] = centered_integrals

known_mets = mets[-grep("unknown", mets)]
# known_mets = mets
monoculture_only = all_data[all_data$species %in% c('356','360','361','492','500','519'),]
# collapse columns into means
mono356 = colMeans(monoculture_only[monoculture_only$species == '356',known_mets])
mono360 = colMeans(monoculture_only[monoculture_only$species == '360',known_mets])
monoculture_for_plotting = rbind(mono356,mono360)
mono361 = colMeans(monoculture_only[monoculture_only$species == '361',known_mets])
monoculture_for_plotting = rbind(monoculture_for_plotting,mono361)
mono492 = colMeans(monoculture_only[monoculture_only$species == '492',known_mets])
monoculture_for_plotting = rbind(monoculture_for_plotting,mono492)
mono500 = colMeans(monoculture_only[monoculture_only$species == '500',known_mets])
monoculture_for_plotting = rbind(monoculture_for_plotting,mono500)
# Exclude 502 for the paper
#mono502 = colMeans(monoculture_only[monoculture_only$species == '502',known_mets])
#monoculture_for_plotting = rbind(monoculture_for_plotting,mono502)
mono519 = colMeans(monoculture_only[monoculture_only$species == '519',known_mets])
monoculture_for_plotting = rbind(monoculture_for_plotting,mono519)
row.names(monoculture_for_plotting) = c('ASF356','ASF360','ASF361','ASF492','ASF500','ASF519')
# try min/max scaling within each column
monoculture_for_plotting_maxnorm = sweep(monoculture_for_plotting,MARGIN=2,apply(monoculture_for_plotting,2,function(x) max(abs(x), na.rm = TRUE)), FUN="/")

# generate a heatmap of monoculture metabolomes using the centered, max/min normalized integrals

pheatmap::pheatmap(monoculture_for_plotting_maxnorm,filename='monoculture_metabolome.jpg',width=9,height=3.2,cluster_rows=FALSE)

#### Random Forest ####
# 360, 361 and 502 are non shifters
# 356, 492, 500, 519 are shifters

## Known Mets Only ##
monoculture_data = all_data[all_data$species %in% c('356','360','361','492','500','502','519'),]
known_mets = mets[-grep("unknown", mets)]
monoculture_data_slim = monoculture_data[,known_mets]
colnames(monoculture_data_slim)[grep("2-Oxoisocaproate", colnames(monoculture_data_slim))] <- "Oxoisocaproate"
monoculture_data_slim$species = monoculture_data$species

monoculture_data_slim <- mutate(monoculture_data_slim, shifter = 
                                  ifelse(species == "356" | species == "492" | 
                                         species == "500" | species == "519", 1, 0))
monoculture_data_slim$shifter <- as.factor(monoculture_data_slim$shifter)
monoculture_data_slim$species <-  NULL

# monoculture_data_slim_norm = sapply(monoculture_data_slim[,1:51], function(x) max(abs(x), na.rm = TRUE))

monoculture_data_slim_norm = sweep(monoculture_data_slim[1:51], MARGIN = 2,
                          STATS = sapply(monoculture_data_slim[,1:51], function(x) max(abs(x), na.rm = TRUE)), FUN="/")
monoculture_data_slim_norm$shifter = monoculture_data_slim$shifter

fit <- AUCRF(shifter~., data = monoculture_data_slim_norm, ranking=c("MDG","MDA"), ntrees = 2000)
fitcv <- AUCRFcv(fit)
OptimalSet(fitcv)

tiff("aucrf_known_mets_only.tiff", width = 6, height = 5, units = 'in', res = 600)
plot(fitcv)
dev.off()

# Box Plot #
data_meta <- monoculture_data_slim_norm

n = fitcv[["Kopt"]]
# n = 10

product_names <- fitcv[["ranking"]][1:n] %>% names() %>% as.data.frame()
colnames(product_names) <- c("ids")

info_feats <- data_meta[,c(as.character(product_names$ids),'shifter')]
info_feats_long <- gather(info_feats, product, likelihood, 1:n, factor_key = TRUE)

# info_feats_long$likelihood <- info_feats_long$likelihood - 1

tiff("known_mets_boxplot.tiff", width = 7, height = 5, units = 'in', res = 600)
theme_set(theme_bw())
ggplot(data = info_feats_long, aes(x = product, y = likelihood, color = shifter)) + 
  geom_boxplot(aes(x = product, y = likelihood, color = shifter), outlier.size = .25, position=position_dodge(width = 1)) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
  scale_color_manual(labels = c("No", "Yes"), values = c("#00BFC4", "#F8766D")) +
  # coord_cartesian(ylim = c(0,1)) +
  labs(x = "Products", y = "Normalized Metabolite Concentration Relative to Fresh Media", title = "Shifter Products", color = "Causes Shift:")
dev.off()


#### Follow-up analysis ####
info_feats_long <- gather(data_meta, product, likelihood, 1:51, factor_key = TRUE)

tiff("all_known_mets_boxplot.tiff", width = 20, height = 5, units = 'in', res = 600)
theme_set(theme_bw())
ggplot(data = info_feats_long, aes(x = product, y = likelihood, color = shifter)) + 
  geom_boxplot(aes(x = product, y = likelihood, color = shifter), outlier.size = .25, position=position_dodge(width = 1)) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
  scale_color_manual(labels = c("No", "Yes"), values = c("#00BFC4", "#F8766D")) +
  # coord_cartesian(ylim = c(0,1)) +
  labs(x = "Products", y = "Normalized Metabolite Concentration Relative to Fresh Media", title = "Shifter Products", color = "Causes Shift:")
dev.off()

#### AUCRF with testable mets ####

# Process data to aquire only testable mets #
monoculture_data_slim_1 <- filter(monoculture_data_slim, shifter == 1)
monoculture_data_slim_1$shifter <- NULL
monoculture_data_slim_0 <- filter(monoculture_data_slim, shifter == 0)
monoculture_data_slim_0$shifter <- NULL
medians_1 <- apply(monoculture_data_slim_1, 2, function(x) median(x, na.rm=TRUE))
medians_0 <- apply(monoculture_data_slim_0, 2, function(x) median(x, na.rm=TRUE))
met_keepers <- names(medians_1)[c(medians_1>medians_0 & medians_1>0)]
# met_keepers

monoculture_data_slim_norm = sweep(monoculture_data_slim[1:51], MARGIN = 2,
                                   STATS = sapply(monoculture_data_slim[,1:51], function(x) max(abs(x), na.rm = TRUE)), FUN="/")
monoculture_data_slim_norm = monoculture_data_slim_norm[,met_keepers]
monoculture_data_slim_norm$shifter = monoculture_data_slim$shifter

# Boxplot of testable mets #
info_feats <- data_meta[,c(met_keepers,'shifter')]
info_feats_long <- gather(info_feats, product, likelihood, 1:length(met_keepers), factor_key = TRUE)
pub_order <- c("Alanine","Lysine","Valine","Hypoxanthine","Urocanate", # Amino Acids and derivatives
               "Acetate","Butyrate","Formate","Isovalerate","Propionate","Succinate") #SCFAs

info_feats_long$product <- factor(info_feats_long$product, levels = pub_order)

# tiff("testable_known_mets_boxplot.tiff", width = 10, height = 5, units = 'in', res = 600)
theme_set(theme_bw())
p = ggplot(data = info_feats_long, aes(x = product, y = likelihood, color = shifter)) + 
    geom_boxplot(aes(x = product, y = likelihood, color = shifter), outlier.size = .25, position=position_dodge(width = 1)) +
    theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
    scale_color_manual(labels = c("No", "Yes"), values = c("#00BFC4", "#F8766D")) +
    # coord_cartesian(ylim = c(0,1)) +
    labs(x = "", y = paste0("Normalized Metabolite Concentration\nRelative to Fresh Media"), title = "", color = "Causes Shift:")
# dev.off()
ggsave("testable_known_mets_boxplot.svg", p, width = 8, height = 4)

# AUCRF with only testable mets #
set.seed(42)
fit <- AUCRF(shifter~., data = monoculture_data_slim_norm, ranking=c("MDG","MDA"), ntrees = 2000)
fitcv <- AUCRFcv(fit)
OptimalSet(fitcv)

# tiff("up_aucrf_testable_known_mets_only.tiff", width = 6, height = 5, units = 'in', res = 600)
# plot(fitcv, ylim=c(0.75, 1.05), color=c("black"))
# dev.off()
df = data.frame(x = c(fitcv[["AUCcurve"]][["k"]]), y = c(fitcv[["AUCcurve"]][["AUC"]]))
theme_set(theme_bw())
p = ggplot(df, aes(x = x, y = y)) + geom_line(color = "black", size = 1.1) + 
  geom_point(size = 2) + geom_point(x=3, y=0.997, color = "red",size=5,alpha=0.2) +
  coord_cartesian(ylim = c(0.8,1.05),xlim=c(1,11)) +
  scale_x_continuous(breaks=seq(1,11,1)) + 
  scale_y_continuous(breaks=seq(0.8,1.05,.05)) +
  annotate("text",x=4.5,y=0.95,label=paste0("cvAUC = 0.988\nOOB-AUCopt = 0.997")) +
  labs(x = "Number of Features in Model", y = "OOB-AUC", title = "")
p
ggsave("up_aucrf_testable_known_mets_only.svg", p, width = 4, height = 3)


# Box Plot #
data_meta <- monoculture_data_slim_norm

n = fitcv[["Kopt"]]
# n = 10

product_names <- fitcv[["ranking"]][1:n] %>% names() %>% as.data.frame()
colnames(product_names) <- c("ids")

info_feats <- data_meta[,c(as.character(product_names$ids),'shifter')]
info_feats_long <- gather(info_feats, product, likelihood, 1:n, factor_key = TRUE)

# tiff("up_testable_known_mets_boxplot.tiff", width = 7, height = 5, units = 'in', res = 600)
theme_set(theme_bw())
p = ggplot(data = info_feats_long, aes(x = product, y = likelihood, color = shifter)) + 
    geom_boxplot(aes(x = product, y = likelihood, color = shifter), outlier.size = .25, position=position_dodge(width = 1)) +
    theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
    scale_color_manual(labels = c("No", "Yes"), values = c("#00BFC4", "#F8766D")) +
    # coord_cartesian(ylim = c(0,1)) +
    labs(x = "", y = paste0("Normalized Metabolite Concentration\nRelative to Fresh Media"), title = "", color = "Causes Shift:")
# dev.off()
ggsave("up_testable_known_mets_boxplot.svg", p, width = 5, height = 4)

#### Auxiliary Stats ####

Samp   <- c("Ace", "Ace","Ace", "Ace", "Ace", "Ace","Iso","Iso","Iso", "Iso", "Iso", "Iso","Pro","Pro","Pro","Pro", "Pro", "Pro", "Pro", "Pro")
Concen <- c("low", "low","low","high","high","high","low","low","low","high","high","high","low","low","low","low","high","high","high","high")
PD     <- c(0.711,0.0578,0.056,  3.93,  3.72,  3.61, 6.17, 4.99, 5.09,  8.72,  8.11,  7.93, 3.36, 1.99, 3.15,  6.2, 13.11, 11.94, 11.54, 10.54)
low_avg<- c(0.275, 0.275,0.275, 0.275, 0.275, 0.275, 5.42, 5.42, 5.42,  5.42,  5.42,  5.42, 3.68, 3.68, 3.68, 3.68,  3.68,  3.68,   3.68, 3.68)

shifter_df <- data.frame(Samp, Concen, PD, low_avg)
shifter_df <- transform(shifter_df, norm = PD-low_avg)

p_val <- wilcox.test(norm ~ Concen, data=shifter_df) %>% .["p.value"] %>% as.numeric()

p = ggplot(shifter_df, aes(x=Concen, y=norm, color=Samp)) +
  geom_point(aes(color = Samp), 
             position = position_jitterdodge(jitter.width = NULL, jitter.height = 0,
                                             dodge.width = 0.7), 
             size = 2, alpha = 0.8) + 
  labs(x = "Concentration Class", y = "Normalized Phase Shift")
  
p
ggsave("concentration_dependent_stat_visual_mets.tiff", p, width = 3, height = 5)
write.csv(shifter_df,"concentration_dependent_stat_visual_mets.csv")


Samp   <- c("ASF492", "ASF492","ASF492","ASF492","ASF492","ASF492","ASF500","ASF500","ASF500","ASF500","ASF500", "ASF500","ASF356","ASF356","ASF356","ASF356","ASF356","ASF356","ASF519","ASF519","ASF519","ASF519","ASF519","ASF519","ASF519")
Concen <- c(   "low",    "low",   "low",  "high",  "high",  "high",   "low",   "low",   "low",  "high",  "high",   "high",   "low",   "low",   "low",  "high",  "high",  "high",   "low",   "low",   "low",  "high",  "high",  "high",  "high")
PD     <- c(   0.963,    0.922,   0.759,    3.16,    4.33,    3.37,    2.66,    1.71,    2.20,    7.50,    7.67,     8.09,    6.40,    7.09,    9.47,   11.66,  10.886,    8.38,    6.20,    8.62,    8.69,  11.082,   15.08,   17.42,   13.58)
low_avg<- c(   0.881,    0.881,   0.881,   0.881,   0.881,   0.881,    2.19,    2.19,    2.19,    2.19,    2.19,     2.19,    7.65,    7.65,    7.65,    7.65,    7.65,    7.65,    7.84,    7.84,    7.84,    7.84,    7.84,    7.84,    7.84)

shifter_df <- data.frame(Samp, Concen, PD, low_avg)
shifter_df <- transform(shifter_df, norm = PD-low_avg)

p_val <- wilcox.test(norm ~ Concen, data=shifter_df) %>% .["p.value"] %>% as.numeric()

p = ggplot(shifter_df, aes(x=Concen, y=norm, color=Samp)) +
  geom_point(aes(color = Samp), 
             position = position_jitterdodge(jitter.width = NULL, jitter.height = 0,
                                             dodge.width = 0.7), 
             size = 2, alpha = 0.8) + 
  labs(x = "Concentration Class", y = "Normalized Phase Shift")

p
ggsave("concentration_dependent_stat_visual_asf.tiff", p, width = 3, height = 5)
write.csv(shifter_df,"concentration_dependent_stat_visual_asf.csv")

