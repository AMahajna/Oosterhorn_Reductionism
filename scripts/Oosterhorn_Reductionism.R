################################################################################
##Create folder for project organization

if(!dir.exists("input_data")){dir.create("input_data")}
if(!dir.exists("output_data")){dir.create("output_data")}
if(!dir.exists("figures")){dir.create("figures")}
if(!dir.exists("scripts")){dir.create("scripts")}

################################################################################
##load packages 

source(file = "scripts/install_load_packages.r")

################################################################################
##Load data

#Microbiome Data
table_raw <- read_csv(file = "input_data/4717_groupby_filtered_all_summary_pathways_taxonomy_abs.csv", show_col_types = FALSE)

#Read Removal Efficiency file 
RE = read_excel("input_data/Removal_Efficiency.xlsx")

################################################################################
##Tidying and subsetting data

#Merge uniprot name and code to have unique names 
table_raw$uniprot_species_name = paste(table_raw$species_strain,table_raw$uniprot_name,table_raw$uniprot_acc, sep="_")

#subset assay data is the experimental result 
assay_data <-
  table_raw %>%  
  dplyr::select(uniprot_species_name, starts_with("47")) %>% 
  rename_with(~str_sub(., end = 8))

#Cleaning "<DL" values= replace less than detection level with zero
assay_data =  data.frame(lapply(assay_data, function(x) {
  gsub("<DL", 0, x)
}))

#Checking the class of the values in the assay data
#sapply(assay_data,class)

#Changing the class of the assay data to numeric
assay_data[,2:33] <- sapply(assay_data[,2:33],as.numeric)

#df_check = assay_data
#Checking class and dimensionality
#sapply(assay_data,class)
#dim(assay_data)

#convert assay data to datafrme format 
assay_data <- 
  assay_data %>% 
  as_tibble() %>% 
  rename_with(~gsub("X","",.x)) %>% 
  as.data.frame()

#subset rownames is for row data that contains gene data 
rownames(assay_data) = assay_data$uniprot_

#Here we take only numeric data and convert to matrix for assay data 
assay_data = assay_data[,-1] %>% as.matrix()

#subset rowData which contain taxonomy  
row_data_tax <-
  table_raw %>% 
  dplyr::select(uniprot_species_name,taxid:domain)

#extract column data which is corresponding to the NGS sample information
samdat = read_excel("input_data/WWTP_overview_samples_20180125.xlsx") %>%
  filter(NGS_performed == "YES") %>% 
  dplyr::select(sampleid = MonsterCodeBC, everything()) 

#additional data 
sample_data <-
  table_raw %>%  
  dplyr::select(starts_with("47")) %>%
  colnames() %>% 
  data.frame(check.names = T)
colnames(sample_data) = "sampleid"
samdf = data.frame(sample_data)

#reorganizing for next step
table_raw <- table_raw %>%
  dplyr::select(uniprot_species_name, everything())

#splitting additional data 
col_data <- 
  table_raw %>% 
  dplyr::select(uniprot_species_name:info) %>% 
  dplyr::select(uniprot_species_name,everything())

#Functional data of enzymes  
row_data_genes <-
  table_raw %>% 
  dplyr::select(uniprot_species_name,kegg_reaction_number:go_number) %>% 
  mutate_all(funs(str_replace(., "_","NA"))) %>% 
  mutate_all(funs(str_replace(., "-$","NA")))  

# reorder taxonomy 
#Error: Taxonomic ranks are not in order. Please reorder columns, which #correspond to taxonomic ranks like this:
#'domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'.
#rename sub_domain to kingdom

#rename columns of taxonomy 
names(row_data_tax)[names(row_data_tax) == 'sub_domain'] <- "kingdom"
names(row_data_tax)[names(row_data_tax) == 'species'] <- "species_missing"
names(row_data_tax)[names(row_data_tax) == 'species_strain'] <- "species"
#colnames(row_data_tax)

#Reorder taxonomy columns 
row_data_tax_reorder = select(row_data_tax,uniprot_species_name, domain, kingdom, phylum, class, order, family, genus,species, everything())

#Capitalize taxonomy column names
colnames(row_data_tax_reorder) <- c('uniprot_species_name', 'Domain', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus','Species',"taxid","species_missing")

################################################################################
##Merge Data

assays = SimpleList(counts = assay_data)
colData = data.frame(samdat)
colData$Year_Sample<-as.character(colData$Year_Sample)
rowData = data.frame(row_data_tax_reorder)

################################################################################
##Constructing TreeSummarizrExperiment 
tse<- TreeSummarizedExperiment(assays = assays,
                               colData = colData,
                               rowData = rowData
)


#Transform assay from counts to relative abundance 
tse <- transformAssay(tse, method = "relabundance")

################################################################################
##Calculate alpha diversity measures 
##Pre-process data by correlation analysis for each alpha diversity measure

################################################################################
#Richness  
tse <- mia::estimateRichness(tse, 
                             assay.type = "counts", 
                             index =   c("ace", "chao1", "hill", "observed"), 
                             name=  c("ace", "chao1", "hill", "observed"))

#colSums(df_check[c("X4717_003", "X4717_008", "X4717_010")] > 0)

Richness = data.frame(matrix( 
  vector(), 32, 5, dimnames=list(c(), c("sampleid", "ace", "chao1", "hill", "observed"))), 
  stringsAsFactors=F)
Richness$sampleid = colData$sampleid
Richness$ace = tse$ace
Richness$chao1 = tse$chao1
Richness$hill = tse$hill
Richness$observed = tse$observed 

Richness_corr <- round(cor(Richness[,2:5]), 1)
testRes_rich = cor.mtest(Richness[,2:5],conf.level = 0.95)

tiff("figures/corr_richness.tiff", units="in", width=5, height=5, res=1000)
ggcorrplot(Richness_corr, hc.order = TRUE, type = "upper", lab = TRUE, 
           outline.col = "white",
           ggtheme = ggplot2::theme_gray,
           colors = c("#6D9EC1", "white", "#E46726"), p.mat = testRes_rich$p)
dev.off()
#Take: observed  
alpha = Richness[,c("sampleid",'observed')]
#colnames(colData(tse))
#head(tse$observed)

################################################################################
#Diversity
tse <- mia::estimateDiversity(tse, 
                              assay.type = "counts",
                              index =  c("coverage", "fisher", "gini_simpson", "inverse_simpson",
                                         "log_modulo_skewness", "shannon"), 
                              name =  c("coverage", "fisher", "gini_simpson", "inverse_simpson",
                                        "log_modulo_skewness", "shannon"))
#head(tse$log_modulo_skewness)
#colnames(colData(tse))

#Exclude log_modulo_skewness because it's a measure of rarity
Diversity = data.frame(matrix(
  vector(), 32, 6, dimnames=list(c(), c('sampleid',"coverage", "fisher", "gini_simpson", "inverse_simpson",
                                        "shannon"))), stringsAsFactors=F)
Diversity$sampleid = colData$sampleid
Diversity$coverage = tse$coverage
Diversity$fisher = tse$fisher
Diversity$gini_simpson = tse$gini_simpson
Diversity$inverse_simpson = tse$inverse_simpson
Diversity$shannon = tse$shannon

Diversity_corr <- round(cor(Diversity[,2:6]), 2)
testRes_diversity = cor.mtest(Diversity[,2:6],conf.level = 0.95)

tiff("figures/corr_diversity.tiff", units="in", width=5, height=5, res=1000)
ggcorrplot(Diversity_corr, hc.order = TRUE, type= "upper", lab = TRUE, 
           outline.col = "white",
           ggtheme = ggplot2::theme_gray,
           colors = c("#6D9EC1", "white", "#E46726"), p.mat= testRes_diversity$p)
dev.off()
#Take: shannon
alpha = cbind(alpha,Diversity[6])

################################################################################
#Evenness 
tse <- estimateEvenness(tse, 
                        assay.type = "counts", 
                        index=c("camargo", "pielou", "simpson_evenness", "evar", "bulla"),
                        name = c("camargo", "pielou", "simpson_evenness", "evar", "bulla"))
#head(tse$camargo)
#colnames(colData(tse))

Evenness = data.frame(matrix(
  vector(), 32, 6, dimnames=list(c(), c('sampleid',"camargo", "pielou",
                                        "simpson_evenness", "evar", "bulla"))),
  stringsAsFactors=F)
Evenness$sampleid = colData$sampleid
Evenness$camargo = tse$camargo
Evenness$pielou = tse$pielou
Evenness$simpson_evenness = tse$simpson_evenness
Evenness$evar = tse$evar
Evenness$bulla = tse$bulla

Evenness_corr <- round(cor(Evenness[,2:6]), 1)
testRes_even = cor.mtest(Evenness[,2:6],conf.level = 0.95)

tiff("figures/corr_evenness.tiff", units="in", width=5, height=5, res=1000)
ggcorrplot(Evenness_corr, hc.order = TRUE, type= "upper", lab = TRUE, 
           outline.col = "white",
           ggtheme = ggplot2::theme_gray,
           colors = c("#6D9EC1", "white", "#E46726"), p.mat = testRes_even$p)
dev.off()
#Take 3: 
alpha = cbind(alpha,Evenness[4])

################################################################################
#Dominance
tse <- estimateDominance(tse, 
                         assay.type = "counts", 
                         index=c("absolute", "dbp", "core_abundance", "gini", "dmn", "relative",
                                 "simpson_lambda"), 
                         name = c("absolute", "dbp", "core_abundance", "gini", "dmn", "relative",
                                  "simpson_lambda"))
#head(tse$relative)
#colnames(colData(tse))

Dominance = data.frame(matrix(
  vector(), 32, 8, dimnames=list(c(), c('sampleid',"absolute","dbp", "core_abundance",
                                        "gini", "dmn", "relative",
                                        "simpson_lambda"))),
  stringsAsFactors=F)
Dominance$sampleid = colData$sampleid
Dominance$absolute = tse$absolute
Dominance$dbp = tse$dbp
Dominance$core_abundance = tse$core_abundance
Dominance$gini = tse$gini
Dominance$dmn = tse$dmn
Dominance$relative = tse$relative
Dominance$simpson_lambda = tse$simpson_lambda

Dominance_corr <- round(cor(Dominance[,2:8]), 1)
testRes_dominance = cor.mtest(Dominance[,2:8],conf.level = 0.95)

tiff("figures/corr_dominance.tiff", units="in", width=5, height=5, res=1000)
ggcorrplot(Dominance_corr, hc.order = TRUE, type= "upper", lab = TRUE, 
           outline.col = "white",
           ggtheme = ggplot2::theme_gray,
           colors = c("#6D9EC1", "white", "#E46726"), p.mat=testRes_dominance$p)
dev.off()
#take:relative 
alpha = cbind(alpha,Dominance[7])

################################################################################
##Rarity 
tse <- mia::estimateDiversity(tse, 
                              assay.type = "counts",
                              index = "log_modulo_skewness")

#Another calculation for rarity 
df = data.frame(assay(tse,2))
#check
#dim(df)
#sapply(df, sum)
rarity_0.01 = colSums(df * (df <0.01))
#Same as above based on counts/reads and not relative abundance 
#rarity_0.01_n = colSums(df != 0)

#Get rare species
#x = df * (df <0.01)
#y = cbind(rowData, x)
#You can loop to get all rare phylum 
#y %>%
#  filter( X4717_003 > 0) %>%
#  select(Domain) %>%
#  unique()

Rarity = data.frame(matrix(
  vector(), 32, 3, dimnames=list(c(), c('sampleid',"log_modulo_skewness", 
                                        "rarity_0.01"))), stringsAsFactors=F)
Rarity$sampleid = colData$sampleid
Rarity$log_modulo_skewness = tse$log_modulo_skewness
Rarity$rarity_0.01 = rarity_0.01

Rarity_corr <- round(cor(Rarity[,2:3]), 2)
testRes_rarity = cor.mtest(Rarity[,2:3],conf.level = 0.95)

tiff("figures/corr_rarity.tiff", units="in", width=5, height=5, res=1000)
ggcorrplot(Rarity_corr, hc.order = TRUE, type= "upper", lab = TRUE, 
           outline.col = "white",
           ggtheme = ggplot2::theme_gray,
           colors = c("#6D9EC1", "white", "#E46726"), p.mat = testRes_rarity$p)
dev.off()
alpha = cbind(alpha,Rarity[ , c("log_modulo_skewness","rarity_0.01")])

################################################################################
##Divergence 

tse <- estimateDivergence(tse)
#colnames(colData(tse))

Divergence  = data.frame(matrix(
  vector(), 32, 2, dimnames=list(c(), c('sampleid',"divergence"))), stringsAsFactors=F)
Divergence$sampleid = colData$sampleid
Divergence$divergence = tse$divergence 
alpha = cbind(alpha,Divergence[2])

################################################################################
#Checking correlation in the resulting alpha diversity measure to check that 
# other types of alpha diversity don't correlate to each other 

alpha_corr <- round(cor(alpha[,2:length(colnames(alpha))]), 1)
testRes_alpha = cor.mtest(alpha[,2:length(colnames(alpha))],conf.level = 0.95)

tiff("figures/corr_alpha_all.tiff", units="in", width=5, height=5, res=1000)
ggcorrplot(alpha_corr,hc.order = TRUE,lab = TRUE,lab_size = 2.5,
           outline.col = "white",
           ggtheme = ggplot2::theme_gray,
           colors = c("#6D9EC1", "white", "#E46726"),
           tl.cex = 10, p.mat = testRes_alpha$p)
dev.off()
alpha_clean = subset(alpha, select = -c(rarity_0.01, relative) )

################################################################################
#cross correlation analysis 

for (i in 2:length(alpha_clean)){ 
  #Add name of the alpha diversity used in CCF
  n=colnames(alpha_clean)
  text = paste("Cross Correlation Function for", n[i], "\n" , "and Removal Efficiency")
  #CCF
  plot <- ggCcf(alpha_clean[, i],
                RE$RE,
                type = "correlation",
                na.action = na.contiguous) +
    theme_minimal() + 
    scale_x_continuous(limits = c(-7, 0), breaks = seq(-7,0,1)) +
    scale_y_continuous(limits = c(-1, 1), breaks = seq(-1,1,0.2)) +
    labs(title=text,
         x ="Lag", y = "Correlation")
  #ggplotly(plot)
  #Save plot 
  cff_name = paste('figures/','ccf_', n[i],".tiff")
  tiff(cff_name, units="in", width=5, height=5, res=1000)
  print(plot)
  dev.off()
}





