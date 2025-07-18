#- decotam nochmal nur mit "empty probe"
#  - nochmal mit empty --> very low concenctraion artificially



library(metadeconfoundR)
library(ggplot2)
library(cowplot)
library(vegan)
library(here)
library(rtk)
library(stringr)
library(dplyr)
library(reshape2)
library(ggrepel)
library(latex2exp)
library(ggpubr)
library(decontam)
library(phyloseq)
library(tidyr)
library(metacal) # devtools::install_github("mikemc/metacal")
library(lme4)
library(lmtest)
source(here::here("../snippets/ggsaveAll.R")) # on Till's local Dropbpx

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

source("~/Dropbox/Package/metadeconfoundR/R/CliffsDelta.R")

#####                                                                #####
##########                                                      ##########
###############             read in features and metadata  ###############
##########                                                      ##########
#####                                                                #####

# genus level raw results (containing contaminant reads)

gen_raw <-
  read.table(
    here::here("input/lotus2_SLV_raw/higherLvl/Genus.txt"),
    sep = "\t",
    header = T,
    row.names = 1
  )

# genus level decontaminated results

gen_decont <-
  read.table(
    here::here("input/lotus2_SLV_cont_rm_raw/higherLvl/Genus.txt"),
    sep = "\t",
    header = T,
    row.names = 1
  )

gen_decontLotus <-
  read.table(
    here::here("input/lotus2_SLV_raw_rmhumanlotus/higherLvl/Genus.txt"),
    sep = "\t",
    header = T,
    row.names = 1
  )

# read in theoretic bacterial compositin and match names to 16S- assigned genera
theoComp <-
  read.table(
    here::here("input/microbial_composition.txt"),
    sep = "\t",
    row.names = 1
  )
theoComp$rel <- theoComp$V2/ sum(theoComp$V2)
theoComp$names <- str_extract(string = rownames(theoComp), pattern = "^[[:alpha:]]+ ")
theoComp$names <- sub(x = theoComp$names, pattern = " ", replacement = "")

theoComp$names[theoComp$names == "Lactobacillus"] <- "Limosilactobacillus"
theoComp$names[theoComp$names == "Escherichia"] <- "Escherichia-Shigella"
rownames(theoComp) <- theoComp$names
theoComp <- as.data.frame(t(theoComp))
theoComp <- theoComp[2, ]
theoComp$Saccharomyces <- NULL
theoComp$Cryptococcus <- NULL
colnames(gen_decont) <- str_match(pattern = "P01.+$", string =  colnames(gen_decont))
colnames(gen_decont) <- sub("_hgrm", "", colnames(gen_decont))
colnames(gen_decont) <- stringr::str_replace_all(pattern = "\\.", 
                                                 replacement = "-",
                                                 string = colnames(gen_decont))

colnames(gen_decont) <- stringr::str_replace_all(pattern = "P01.E06.82", 
                                                 replacement = "P01.E06.82.empty",
                                                 string = colnames(gen_decont))

colnames(gen_decont) <- stringr::str_replace_all(pattern = "-", 
                                                 replacement = ".",
                                                 string = colnames(gen_decont))
colnames(gen_raw) <- colnames(gen_decont)
colnames(gen_decontLotus) <- colnames(gen_decont)

# relAbun <- cbind((rowSums(gen_raw)/sum(rowSums(gen_raw)))*100,
# (rowSums(gen_decont)/sum(rowSums(gen_decont)))*100)


metadata <-
  read.table(
    here::here("input/METADATA.csv"),
    sep = ",",
    header = T#,
    #row.names = 1
  )
metadata$Sample_name <- NULL
rownames(metadata) <- metadata$Sample_identifier
metadata$Sample_identifier <- NULL

rownames(metadata) <- stringr::str_replace_all(pattern = "\\.", 
                                                 replacement = "",
                                                 string = rownames(metadata))
rownames(metadata) <- stringr::str_replace_all(pattern = "-", 
                                                 replacement = "\\.",
                                                 string = rownames(metadata))

metadata$ratio <- metadata$human/metadata$bacteria
metadata$ratio[metadata$ratio == Inf] <- 100

metadata <- metadata[order(rownames(metadata)), ]
metadata$pseudo <- rbinom(45, 1, 0.5)
metadata$group <- substr(str_extract(pattern = "[0-9]{1,2}$", 
                                     string = rownames(metadata)),
                         start = 1,
                         stop = 1)
metadata["rel", ] <- NA
metadata$group[is.na(metadata$group)] <- "theoretical"
metadata["P01.E06.82.empty", "group"] <- "empty"
metadata$group[metadata$ratio == 100] <- "pure_human"
#rownames(metadata) <- metadata$Row.names
#metadata$Row.names <- NULL

metadata$ratio <- round(x = metadata$ratio, digits = 2)
metadata$ratio[metadata$ratio == "0.01"] <- "99:1"
metadata$ratio[metadata$ratio == "0.11"] <- "90:10"
metadata$ratio[metadata$ratio == "0.33"] <- "75:25"
metadata$ratio[metadata$ratio == "1"] <- "50:50"
metadata$ratio[metadata$ratio == "3"] <- "25:75"
metadata$ratio[metadata$ratio == "9"] <- "10:90"
metadata$ratio[metadata$ratio == "99"] <- "1:99"
metadata$ratio[metadata$ratio == "100"] <- "pure_human"


metadata["P01.E06.82.empty", c("ratio")] <- "empty"
metadata["rel", c("ratio")] <- "theoretical"


metadata$ConcGroup <- NA
metadata$ConcGroup[metadata$Concentration <= 1.5] <- 1
metadata$ConcGroup[metadata$Concentration > 1.5] <- 2
metadata$ConcGroup[metadata$Concentration > 3] <- 3
metadata$ConcGroup[metadata$Concentration > 6] <- 4
metadata$ConcGroup[metadata$Concentration > 12] <- 5
metadata$ConcGroup[metadata$Concentration > 25] <- 6

metadata$ConcGroup[metadata$group == "empty"] <- "empty"

metadata$ConcGroup[metadata$ConcGroup == 1] <-
  sprintf("%04g", round(mean(metadata$Concentration[metadata$ConcGroup == 1], na.rm = T), digits = 1))
metadata$ConcGroup[metadata$ConcGroup == 2] <- 
  sprintf("%04g", round(mean(metadata$Concentration[metadata$ConcGroup == 2], na.rm = T), digits = 1))
metadata$ConcGroup[metadata$ConcGroup == 3] <- 
  sprintf("%04g", round(mean(metadata$Concentration[metadata$ConcGroup == 3], na.rm = T), digits = 1))
metadata$ConcGroup[metadata$ConcGroup == 4] <- 
  sprintf("%04g", round(mean(metadata$Concentration[metadata$ConcGroup == 4], na.rm = T), digits = 1))
metadata$ConcGroup[metadata$ConcGroup == 5] <- 
  sprintf("%04g", round(mean(metadata$Concentration[metadata$ConcGroup == 5], na.rm = T), digits = 1))
metadata$ConcGroup[metadata$ConcGroup == 6] <- 
  sprintf("%04g", round(mean(metadata$Concentration[metadata$ConcGroup == 6], na.rm = T), digits = 1))

metadata$ConcGroup[metadata$group == "theoretical"] <- "theoretical"
metadata$ConcGroup[metadata$group == "pure_human"] <- "pure_human"

metadata$group[metadata$group == 8] <- "pure_bact"
metadata$ConcGroup[metadata$group == "pure_bact"] <- "pure_bact"
metadata$ratio[metadata$group == "pure_bact"] <- "pure_bact"

metadata$blinds <- metadata$ConcGroup
metadata$blinds[metadata$blinds %in% c("pure_human"#, 
                                   # "empty"
)] <- "blinds"
metadata$blinds[metadata$blinds != "blinds"] <- "smples"

#metadata_red <- metadata[rownames(metadata) %in% rownames(sets$gen_raw_rel), ]
#metadata_red <- metadata_red[, c("pseudo", "Concentration", "human", "bacteria", "ratio", "group", "ConcGroup")]




#####                                                                #####
##########                                                      ##########
###############  read in phyloseq obj OTUs and de decontam  ###############
##########                                                      ##########
#####                                                                #####
# read in phyloseq object from otu.biom
all_otu <- read.table(here::here("input/lotus2_SLV_raw/higherLvl/OTU.txt"),
                      sep = "\t",
                      header = T,
                      row.names = 1)
colnames(all_otu) <- colnames(gen_raw)
otu <- otu_table(all_otu, taxa_are_rows = T)  #check sample names

warning("first create AbundanceAll further down in ths script without the Decontam dataset, then come back up here, create this dataset and rerun AbundanceAll computation!")
sample <- sample_data(metadata) 
warning("changing concentration of empty sample from 0 to 0.01!!!")
sample$Concentration[sample$Concentration == 0] <- 0.01
sample <- sample[-nrow(sample), ]
#check sample names  <- muss die gleichen sein wie in deinen otu object
##Check which samples are not in the out matrix
setdiff(sample_names(sample), sample_names(otu))
biom_file <- import_biom(here::here("input/lotus2_SLV_raw/OTU.biom"))
tax <- tax_table(biom_file) 

# save it all in a merged new phyloseq object
ps <- merge_phyloseq(otu, sample, tax)

df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=blinds)) + geom_point()


sample_data(ps)$is.neg <- sample_data(ps)$blinds == "blinds"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
contamdf.comb <- isContaminant(seqtab = ps, method="combined", neg="is.neg",con = "Concentration")
table(contamdf.prev$contaminant)
which(contamdf.prev$contaminant)
rownames(contamdf.prev)[contamdf.prev$contaminant]

table(contamdf.comb$contaminant)
which(contamdf.comb$contaminant)
rownames(contamdf.comb)[contamdf.comb$contaminant]

# show percentage of reads removed by filtering out Decontam flagged OTUs
sort((colSums(otu_table(ps)[rownames(contamdf.prev)[contamdf.prev$contaminant], ]) / colSums(otu_table(ps))) * 100)

# using only the "empty" sample as contamination control leads to no OTUs being reported as contaminants
# using both "empty" and "pure human" as contaminant control leads flagging of 25 and 4 OTUs for prevalence and combined approach respectively. 
  # With "OTU46" belonging to Cutibctarium being the only somewhat impactfull hit when using "prevalence" Decontam method

# using only "pure human" as control leads to 166 and 8 OTUs being flagged with prevalence and combined approach respectively.



# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$blinds == "blinds", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$blinds == "smples", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(allTaxa, physeq))
}
# remove OTUs identified as contaminants
cleaned_physeq <- pop_taxa(ps, rownames(contamdf.prev)[contamdf.prev$contaminant])

genus <- tax_glom(physeq = cleaned_physeq, taxrank = "Rank6")
data_biom_all <- as.data.frame(tax_table(cleaned_physeq)) # still includes NA's
data_biom_tax <- as.data.frame(tax_table(genus)) #biome dataframe mit taxa -> damit die zuordnung klappt beim mergen
data_biom <-as.data.frame(otu_table(genus)) #biome datafram mit otu

# umbennen der Ranks in phylum etc.
colnames(data_biom_tax)<-c("kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
data_biom_tax$otu <- row.names(data_biom_tax)

# genus_otu_short <- as.data.frame(otu_table(data_biom_tax))
data_biom_tax$hash <- rownames(data_biom_tax)
data_biom$hash <- rownames(data_biom)
joined_genus <- left_join(data_biom_tax, data_biom)
joined_genus$hash <- NULL
joined_genus$Genus[joined_genus$otu == "OTU15"] <- "Staphylococcaceae;?"
joined_genus$Genus[joined_genus$otu == "OTU119"] <- "Paenibacillaceae;?"
joined_genus$Genus[joined_genus$otu == "OTU114"] <- "Bacilli;?;?;?"
joined_genus$Genus[joined_genus$otu == "OTU17"] <- "noHit"
joined_genus <- rbind(joined_genus, 0)
joined_genus$Genus[nrow(joined_genus)] <- "Cutibacterium"

genus <- joined_genus[, - c(1:5,7) ]
genus$Genus <- sub("g__", "",genus$Genus)
#genus$C = paste(genus$otu, genus$Genus, sep="_")
rownames(genus) <- make.names(genus$Genus, unique = T)
genus$Genus <- NULL
genus$otu <- NULL
#genus$C <- NULL
genus <- genus[order(rownames(genus)), ]
#genus <- as.data.frame(t(genus))
genus["Staphylococcus", ] <- genus["Staphylococcus", ] + genus["Staphylococcaceae..", ]
genus <- genus[(rownames(genus) != "Staphylococcaceae.."), ]

#gen_fraction <- gen_decont / gen_raw
#rowSums(gen_decont) /2106626
warning("adding the two staphylococcaceae manually!!!")

preSets <- list("gen_raw" = gen_raw, 
                "gen_decont" = gen_decont, 
                "gen_decontLotus" = gen_decontLotus, 
                "gen_rawDecontam" = genus
                )
gen_abs <- list()
gen_rel <- list()
sets <- list()
for (preSet in names(preSets)) {
  print(preSet)
  temp <- as.data.frame(t(preSets[[preSet]]))
  if (!is.null(temp$`Bacteria;Firmicutes;Bacilli;Staphylococcales;Staphylococcaceae;Staphylococcus`)) {
    temp$`Bacteria;Firmicutes;Bacilli;Staphylococcales;Staphylococcaceae;Staphylococcus` <- 
      temp$`Bacteria;Firmicutes;Bacilli;Staphylococcales;Staphylococcaceae;?` + 
      temp$`Bacteria;Firmicutes;Bacilli;Staphylococcales;Staphylococcaceae;Staphylococcus`
    temp <- temp[, colnames(temp) != "Bacteria;Firmicutes;Bacilli;Staphylococcales;Staphylococcaceae;?"]
  }
  abs <- paste0(preSet, "_abs")
  rel <- paste0(preSet, "_rel")
  # gen_rel[[preSet]] <- as.data.frame(t(apply(temp, 1, function (x) x/sum(x))))
  # gen_abs[[preSet]] <- temp
  sets[[abs]] <- temp
  sets[[rel]] <- as.data.frame(t(apply(temp, 1, function (x) x/sum(x))))
}


# samplesize_raw <- sort(colSums(gen_raw))[2] # take second lowast and discard empty control
# # remove too small samples
# gen_raw_red <- gen_raw[, colSums(gen_raw) >= samplesize_raw]
# rtk_gen_raw <-
#   rtk(
#     gen_raw_red,
#     repeats = 10,
#     depth = samplesize_raw,
#     ReturnMatrix = 1,
#     margin = 2,
#     verbose = FALSE,
#     threads = 1,
#     tmpdir = NULL,
#     seed = 17
#   )
# rtk_use_gen_raw <-
#   as.data.frame(rtk_gen_raw$raremat)
# colSums(rtk_use_gen_raw)
# gen_raw_feat <- as.data.frame(t(rtk_use_gen_raw))



# 
# samplesize_dec <- sort(colSums(gen_decont))[3] # take third lowest and discard empty control and P01-D06-81
# # remove too small samples
# gen_dec_red <- gen_decont[, colSums(gen_decont) >= samplesize_dec]
# 
# rtk_gen_dec <-
#   rtk(
#     gen_dec_red,
#     repeats = 10,
#     depth = samplesize_dec,
#     ReturnMatrix = 1,
#     margin = 2,
#     verbose = FALSE,
#     threads = 1,
#     tmpdir = NULL,
#     seed = 17
#   )
# rtk_use_gen_dec <-
#   as.data.frame(rtk_gen_dec$raremat)
# colSums(rtk_use_gen_dec)
# gen_dec_feat <- as.data.frame(t(rtk_use_gen_dec))

#####                                                                #####
##########                                                      ##########
###############                    metadeconf?                    ###############
##########                                                      ##########
#####                                                                #####

# test <- MetaDeconfound(gen_raw_feat, metadata_red, nnodes = 5, returnLong = T)
# test_plot <- BuildHeatmap(test, d_range = "full")
# ggsave2(plot = test_plot, filename = here::here("figures/gen_raw_metadeconf.svg"), width = 8, height = 4)
# 
# test2 <- MetaDeconfound(gen_dec_feat, metadata_red[rownames(gen_dec_feat), ], nnodes = 5, returnLong = T)
# test_plot2 <- BuildHeatmap(test2, d_range = "full")
# ggsave2(plot = test_plot2, filename = here::here("figures/gen_decon_metadeconf.svg"), width = 8, height = 4)

#this is the real one!!
# test3rel <- MetaDeconfound(gen_raw_rel, metadata_red[rownames(gen_raw_rel), ], nnodes = 5, returnLong = T)
# test_plot3rel <- BuildHeatmap(test3rel, d_range = "full")
# ggsave2(plot = test_plot3rel, filename = here::here("figures/gen_raw_rel_metadeconf.svg"), width = 8, height = 4)



#####                                                                #####
##########                                                      ##########
###############               relative abundance           ###############
##########                                                      ##########
#####                                                                #####

warning("ongoing construction. plug in elements of setsAbundanceAll into script below loop")

setsAbundanceAll <- list()


targets <-
  c("Bacillus",
    "Cutibacterium",
    "Enterococcus",
    "Escherichia-Shigella",
    "Limosilactobacillus",
    "Listeria",
    "noHit",
    "Pseudomonas",
    "Salmonella",
    "Staphylococcus"
  )

for (set in names(sets)) {
  print(set)
  dataType <- set
  gen_raw_rel_renamed <- sets[[set]]
  if (!(set %in% c("gen_rawDecontam_abs", "gen_rawDecontam_rel"))) {
    colnames(gen_raw_rel_renamed) <- sub(pattern = ";", 
                                         replacement = "", 
                                         x = str_extract(string = colnames(gen_raw_rel_renamed), 
                                                         pattern = ";([[:alpha:], [-, ?]]+$)"))
    colnames(gen_raw_rel_renamed)[is.na(colnames(gen_raw_rel_renamed))] <- "noHit"
  }
  tgen <- as.data.frame(t(gen_raw_rel_renamed))
  rownames(tgen)[rownames(tgen) == "Escherichia.Shigella"] <- "Escherichia-Shigella"
  print(tgen[!(rownames(tgen) %in% targets), c("P01.G03.35"), drop = F])
  tgen <- rbind(tgen[(rownames(tgen) %in% targets), ], colSums(tgen[!(rownames(tgen) %in% targets), ]))
  rownames(tgen)[nrow(tgen)] <- "other"
  
  ttheo <- as.data.frame(t(theoComp))
  ttheo$rel <- as.numeric(ttheo$rel)
  gen_raw_merged <- merge(tgen, ttheo, by = 0, all = T)
  rownames(gen_raw_merged) <- gen_raw_merged$Row.names
  gen_raw_merged$Row.names <- NULL
  
  AbundanceAll <- merge(t(gen_raw_merged), metadata, by = 0, all = T)
   rownames(AbundanceAll) <- AbundanceAll$Row.names
   AbundanceAll$Row.names <- NULL

  summary(as.factor(AbundanceAll$group))
  summary(as.factor(AbundanceAll$ratio))
  
  AbundanceDetail <- AbundanceAll
  setsAbundanceAll[[set]] <- AbundanceAll
  
}


#####                                                                #####
##########                                                      ##########
###############             metacal unbiasing              ###############
##########                                                      ##########
#####                                                                #####
# library(tidyverse)
# 
# obsAb <- as.matrix(setsAbundanceAll$gen_decont_abs[setsAbundanceAll$gen_decont_abs$ratio %in% c("pure_bact"), targets[-c(2,7)]])
# actAb <- as.matrix(setsAbundanceAll$gen_decont_abs[setsAbundanceAll$gen_decont_abs$ratio == "theoretical", targets[-c(2,7)]])
# # actAb <- rbind(actAb, actAb)
# # actAb <- rbind(actAb, actAb)
# # actAb <- rbind(actAb, actAb[1:3, ])
# 
# 
# rownames(actAb) <- rownames(obsAb)
# # bias_fit <- estimate_bias(observed = obsAb[7, ,drop = F], actual = actAb[7, ,drop = F], margin = 1, boot = T)
# bias_fit <- estimate_bias(observed = obsAb, actual = actAb, margin = 1, boot = T)
# coef_tb <- summary(bias_fit)$coefficients
# coef_tb %>%
#   mutate(taxon = fct_reorder(taxon, estimate)) %>%
#   ggplot(aes(taxon, estimate, 
#              ymin = estimate / gm_se^2, ymax = estimate * gm_se^2)) +
#   geom_hline(yintercept = 1, color = "grey") +
#   geom_pointrange() +
#   scale_y_log10() +
#   theme_classic() +
#   labs(title = "sequencing efficiency estimates") +
#   coord_flip()
#####                                                                #####
##########                                                      ##########
###############                    rela bund groups        ###############
##########                                                      ##########
#####                                                                #####
warning("sum up other bacteria to one other column")


non_bactNames <- c("pseudo", "Concentration", "human" ,"bacteria", "H2O_for_verdünnung",
                   "Stan_for_verdünnung", "Komment", "ratio", "pseudo" , "group", "ConcGroup", "blinds")

keppGen <- c(colnames(theoComp), "noHit", "Cutibacterium", "other")

cbPalette <-
  c("#999999",
    "#E69F00",
    "#bf0d92",
    "#56B4E9",
    "#082877", 
    "#009E73",
    "#5b4106",
    "#F0E442",
    "#CC79A7",
    "#ffabab",
    "#bfff92")


#####                                                                #####
##########                                                      ##########
###############                 rel abund looped           ###############
##########                                                      ##########
#####                                                                #####

relAbunPlots <- list()
noHitSubsList <- list()
for (k in c("group", "ConcGroup", "ratio")) {
  print(k)
  if (k == "group") {
    next
  }
  
  relAbunPlots[[k]] <- list()
  
  for (AbAll in names(setsAbundanceAll)) {
    if (grepl("_abs", AbAll)) {
      next
    }
    #print(AbAll)
    AbundanceAll <- setsAbundanceAll[[AbAll]]
    abundance_groups <- list()
    for (i in unique(AbundanceAll[[k]])) {
      #print(i)
      abundance_groups[[i]] <-
        colSums(AbundanceAll[AbundanceAll[[k]] == i, !(colnames(AbundanceAll) %in% non_bactNames)])
    }
    
    AbundanceAll_merged <-
      as.data.frame(bind_rows(abundance_groups, .id = "column_label"))
    rownames(AbundanceAll_merged) <- AbundanceAll_merged$column_label
    AbundanceAll_merged$column_label <- NULL
    
    
    rel_PhylumAll <-
      apply(AbundanceAll_merged, 1, function(x)
        x / sum(x, na.rm = T)) # rows 1 als prozente
    # er tranformiert automatisch
    
    rel_PhylumAll <-
      rel_PhylumAll[keppGen, ]
    
    rel_PhylumAll <- as.data.frame(rel_PhylumAll)
    rel_PhylumAll$Phylum <- row.names(rel_PhylumAll)
    
    for (i in 1:(nrow(rel_PhylumAll) - 1)) {
      rel_PhylumAll[, i] <- as.numeric(rel_PhylumAll[, i])
    } # macht aus den abundance für die bacteria numerical zahlen
    
    molten_rel_PhylumAll <-
      melt(data = rel_PhylumAll, id.vars = "Phylum")
    # molten_rel_PhylumAll$labels <-
    #   scales::percent(molten_rel_PhylumAll$value)
    molten_rel_PhylumAll$variable <-
      str_replace((molten_rel_PhylumAll$variable), "Abundance_", "")
    molten_rel_PhylumAll$variable <- sub("pure_bact","pure bact", molten_rel_PhylumAll$variable)
    print(AbAll)
    noHitSubs <- molten_rel_PhylumAll[molten_rel_PhylumAll$Phylum == "noHit", ]
    noHitSubs$value <- scales::percent(noHitSubs$value)
    print(noHitSubs)
    noHitSubsList[[k]][[AbAll]] <- noHitSubs
    if (k == "ratio") {
      molten_rel_PhylumAll$variable <-
        factor(
          molten_rel_PhylumAll$variable,
          levels = c(
            "theoretical",
            "pure bact",
            "99:1",
            "90:10",
            "75:25",
            "50:50",
            "25:75",
            "10:90",
            "1:99",
            "pure_human",
            "empty"
          )
        )
    } else if (k == "ConcGroup") {
      molten_rel_PhylumAll$variable <-
        factor(
          molten_rel_PhylumAll$variable,
          levels = c(
            "theoretical",
            "pure_bact",
            "01.3",
            "02.5",
            "05.2",
            "09.8",
            "19.9",
            "40.7",
            "pure_human",
            "empty"
          )
        )
    }
    
    molten_rel_PhylumAll <- molten_rel_PhylumAll[!(molten_rel_PhylumAll$variable %in% c("theoretical", "pure_human", "empty")), ]
    
    # xlabText <-
    #   switch(k,
    #          "group" = "ratio group as assigned by Theda",
    #          "ConcGroup" = "abs. conc. group",
    #          "ratio" = "ratio group bacterial:human")
    xlabText <- sub(pattern = "gen_", replacement = "", x = AbAll)
    xlabText <- sub(pattern = "_rel", replacement = "", x = xlabText)
    relAbunGroup <-
      ggplot (data = molten_rel_PhylumAll, aes(x = variable, y = value, fill = Phylum)) +
      geom_bar(stat = "identity", color = "black") +
      scale_fill_manual(values = cbPalette) +
      theme_classic() +
      theme(
        #axis.text.x = element_blank(),
        legend.text = element_text(face = "italic"),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1
        ),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank()
      ) +
      xlab(paste0(xlabText)) +
      labs(fill = "genus")
    relAbunPlots[[k]][[AbAll]] <- relAbunGroup
  }
}

percNhoHit <- do.call("rbind", noHitSubsList$ratio)
percNhoHit$rownames <- rownames(percNhoHit)
percNhoHit$rownames <- sub("\\.[0-9]{1,3}", "", percNhoHit$rownames)
write.table(percNhoHit, file = "intermediate/percNoHit.tsv", sep = "\t")
reshape(percNhoHit, timevar = "rownames", direction = "wide")
# extract legend for later plotting in plot_grid
relAbunPlots[["ConcGroup"]][["gen_raw_rel"]] <- relAbunPlots[["ConcGroup"]][["gen_raw_rel"]] +
  theme(legend.background = element_rect(fill = 'transparent'))
legendPlot <- as_ggplot(get_legend(relAbunPlots[["ConcGroup"]][["gen_raw_rel"]], position = "right"))
legendPlot <-
  legendPlot + theme(
    plot.background = element_rect(fill = 'transparent'),
    panel.background = element_rect(fill = 'transparent'),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  )

legendPlotHor <- as_ggplot(get_legend(relAbunPlots[["ConcGroup"]][["gen_raw_rel"]], position = "bottom"))
legendPlotHor <-
  legendPlotHor + theme(
    plot.background = element_rect(fill = 'transparent'),
    panel.background = element_rect(fill = 'transparent'),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  )

# remove legend from all the plots
relAbunPlotsNL <- list()
for (i in names(relAbunPlots)) {
  for (j in names(relAbunPlots$ConcGroup)) {
    relAbunPlotsNL[[i]][[j]] <- relAbunPlots[[i]][[j]] + theme(legend.position = "none")
  }
}
titleConcGroup <- ggdraw() + 
  draw_label(
    "relative genera abundance per\ndilution group (mean absolute DNA concentration in ng/μL)",
    angle = 0,
  ) 
mergedPlots_rel_ConcGroup <- plot_grid(relAbunPlotsNL[["ConcGroup"]][["gen_raw_rel"]], 
                                       relAbunPlotsNL[["ConcGroup"]][["gen_decont_rel"]],
                                       relAbunPlotsNL[["ConcGroup"]][["gen_decontLotus_rel"]], 
                                       relAbunPlotsNL[["ConcGroup"]][["gen_rawDecontam_rel"]],
                                       nrow = 1)
mergedPlots_rel_ConcGroup <- plot_grid(titleConcGroup,
                                       mergedPlots_rel_ConcGroup,
                                       nrow = 2,
                                       rel_heights = c(0.1, 1))


titleRatio <- ggdraw() + 
  draw_label(
    "relative genera abundance per\nbacterial to human DNA ratio group",
    angle = 0,
  )
mergedPlots_rel_ratio <- plot_grid(relAbunPlotsNL[["ratio"]][["gen_raw_rel"]], 
                                   relAbunPlotsNL[["ratio"]][["gen_decont_rel"]],
                                   relAbunPlotsNL[["ratio"]][["gen_decontLotus_rel"]],
                                   relAbunPlotsNL[["ratio"]][["gen_rawDecontam_rel"]],
                                   nrow = 1)
mergedPlots_rel_ratio <- plot_grid(titleRatio, 
                                   mergedPlots_rel_ratio, 
                                   nrow = 2,
                                   rel_heights = c(0.1, 1))




mergedPlots_rel <- plot_grid(mergedPlots_rel_ConcGroup,
                             mergedPlots_rel_ratio, nrow = 2)
mergedPlots_rel <- plot_grid(mergedPlots_rel, legendPlot, nrow = 1, rel_widths = c(1,0.25))
mergedPlots_rel

# 
# ggsaveAll(
#   here::here("figures/"),
#   filename = paste0("relAbund_2by3"),
#   plot = mergedPlots_rel,
#   devices = c("svg", "pdf"),
#   height = 15,
#   width = 15
# )


mergedPlots_rel_ratio <- plot_grid(relAbunPlotsNL[["ratio"]][["gen_raw_rel"]] + xlab("naive"), 
                                   relAbunPlotsNL[["ratio"]][["gen_decont_rel"]] + xlab("filtered pre-clustering"),
                                   relAbunPlotsNL[["ratio"]][["gen_decontLotus_rel"]] + xlab("filtered post-clustering"),
                                   #relAbunPlotsNL[["ratio"]][["gen_rawDecontam_rel"]],
                                   nrow = 1)
mergedPlots_rel_ratio

ggsaveAll(
  here::here("figures/"),
  filename = paste0("fig_4_draft"),
  plot = mergedPlots_rel_ratio,
  devices = c("svg", "pdf"),
  height = 7,
  width = 6
)


#####                                                                #####
##########                                                      ##########
###############                    rel per sample              ###############
##########                                                      ##########
#####                                                                #####

#warning("subset here manually to only ratio 99 etc.")
perSampleAbunPlots <- list()
molten_DF_list <- list()
for (set in names(setsAbundanceAll)) {
  print(set)
  rel_PhylumAllRat <- as.data.frame(t(setsAbundanceAll[[set]]))
  rel_PhylumAllRat <-
    rel_PhylumAllRat[keppGen, ]
  rel_PhylumAllRat$Phylum <- row.names(rel_PhylumAllRat)
  
  
  for (i in 1:(ncol(rel_PhylumAllRat)-1)) {
    rel_PhylumAllRat[, i] <- as.numeric(rel_PhylumAllRat[, i])
  } # macht aus den abundance für die bacteria numerical zahlen
  
  molten_rel_PhylumAllRat <- melt(data = rel_PhylumAllRat, id.vars = "Phylum")
  # molten_rel_PhylumAllRat$labels <- scales::percent(molten_rel_PhylumAllRat$value)
  molten_rel_PhylumAllRat$variable <-
    str_replace((molten_rel_PhylumAllRat$variable), "Abundance_", "")
  molten_rel_PhylumAllRat$ratio <- AbundanceAll[molten_rel_PhylumAllRat$variable, "ratio"]
  molten_rel_PhylumAllRat$ConcGroup <- AbundanceAll[molten_rel_PhylumAllRat$variable, c("ConcGroup")]
  molten_rel_PhylumAllRat$sampleType <- paste0(molten_rel_PhylumAllRat$ratio, " ", molten_rel_PhylumAllRat$ConcGroup)
  molten_rel_PhylumAllRat$sampleType <- sub("pure_human pure_human", "pure human", molten_rel_PhylumAllRat$sampleType)
  molten_rel_PhylumAllRat$sampleType <- sub("pure_bact pure_bact", "pure bact", molten_rel_PhylumAllRat$sampleType)
  molten_rel_PhylumAllRat$sampleType <- sub("empty empty", "empty", molten_rel_PhylumAllRat$sampleType)
  molten_rel_PhylumAllRat$sampleType <- sub("theoretical theoretical", "theoretical", molten_rel_PhylumAllRat$sampleType)
  
  relAbunIndi <-
    ggplot (data = molten_rel_PhylumAllRat, aes(x = sampleType, y = value, fill = Phylum)) +
    geom_bar(stat = "identity", color = "black", show.legend = FALSE) +
    scale_fill_manual(values=cbPalette) +
    theme_classic() +
    theme(axis.title.y = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5),
      axis.text.y = element_blank(),
      axis.ticks = element_blank()) +
    xlab(set) +
    labs(fill = "genus")
  perSampleAbunPlots[[set]] <- relAbunIndi
  molten_DF_list[[set]] <- molten_rel_PhylumAllRat
  if (grepl("abs", set)) {
    x <- molten_rel_PhylumAllRat[molten_rel_PhylumAllRat$variable == "P01.H01.8", ]
    x <- x[!(is.na(x$value)), ]
    p <- molten_rel_PhylumAllRat[molten_rel_PhylumAllRat$variable == "rel", ]
    p <- p[!(is.na(p$value)), ]
    x <- x[x$Phylum %in% p$Phylum, ]
    P <- p[p$Phylum %in% x$Phylum, ]
    xp <- merge(x[, c("Phylum", "value")], p[, c("Phylum", "value")], by = "Phylum")
    #print(paste0("p-value pure-bact vs. rel: ", chisq.test(x = x$value, p = p$value)$p.value))
    print("p-value pure-bact vs. rel:")
    print(chisq.test(x = x$value, p = p$value))
  }
}


# 2023 11 30
# compare each sample to the pure bact sample
# ATTENTION: putting in smpl2 as either y or p argument in chisquare gives extremely different p-values!!!!!
smplNames <- unique(molten_DF_list$gen_decontLotus_abs$variable)
#smplNames <- smplNames[smplNames != "rel"]
chisDF <- as.data.frame(crossing(smplNames, smplNames))
colnames(chisDF) <- c("smpl1", "smpl2")
chisDF$chisquare_p <- NA
chisDF$identifier <- paste0(chisDF$smpl1, chisDF$smpl2)
chisDF <- chisDF[chisDF$smpl2 == "P01.H01.8", ]
#smplNames <- smplNames[smplNames != "P01.H01.8"]


for (smpl1 in smplNames) {
  print(smpl1)
  if (smpl1 %in% c("rel", "P01.H01.8")) {
    next
  }
  x <-
    molten_DF_list$gen_decontLotus_abs[molten_DF_list$gen_decontLotus_abs$variable == smpl1,]
  x <- x[!(is.na(x$value)),]
  x <- x[x$value/sum(x$value) > 0.01, ]
  smpl2 <- "P01.H01.8"
  p <-
    molten_DF_list$gen_decontLotus_abs[molten_DF_list$gen_decontLotus_abs$variable == smpl2,]
  p <- p[!(is.na(p$value)),]
  x <- x[x$Phylum %in% p$Phylum,]
  p <- p[p$Phylum %in% x$Phylum,]
  xp <-
  merge(x[, c("Phylum", "value")], p[, c("Phylum", "value")], by = "Phylum")
  colnames(xp) <- c("Phylum", "value_smpl1", "value_smpl2")
  xp_tab <- as.table(t(xp[, c("value_smpl1", "value_smpl2")]))
  dimnames(xp_tab) <- list(Sample = c("Sample1", "pure_bact"), 
                           Phylum = xp$Phylum)
  print(xp_tab)
  x_tab <- as.table(x$value)
  dimnames(x_tab) <- list(Phylum = xp$Phylum)
  p_tab <- as.table(p$value)
  dimnames(p_tab) <- list(Phylum = xp$Phylum)
  chisDF$chisquare_p[chisDF$identifier == paste0(smpl1, smpl2)] <-
    chisq.test(xp_tab, simulate.p.value = F)$p.value
  
  x_val <- x$value
  names(x_val) <- x$Phylum
  p_val <- p$value
  names(p_val) <- p$Phylum
  
  # chisq.test(x_tab, simulate.p.value = F)
  # chisq.test(x_val, simulate.p.value = F)
  # chisq.test(x = xp_tab, simulate.p.value = F)
  # chisq.test(x = x_tab, y = p_tab, simulate.p.value = F)
  # chisq.test(x = x_tab, p = p_tab/sum(p_tab))
  # chisq.test(x = x_val, y = p_val, simulate.p.value = F)
  # chisq.test(x = x_val, p = p_val/sum(p$value), simulate.p.value = F)
  # chisq.test(x = x$value, y = p$value, simulate.p.value = F)
  
  
  # chisDF$chisquare_p[chisDF$identifier == paste0(smpl1, smpl2)] <-
  #   chisq.test(x = x$value, y = p$value)$p.value
  # chisDF$chisquare_p[chisDF$identifier == paste0(smpl1, smpl2)] <-
  #   chisq.test(x = x$value, p = (p$value/sum(p$value)))$p.value
}

chisDF <- chisDF[!(is.na(chisDF$chisquare_p)), ]
chisDF$FDR <- p.adjust(chisDF$chisquare_p, method = "fdr")

perSampleNiceNames <- molten_rel_PhylumAllRat
unique_niceNames <- perSampleNiceNames[perSampleNiceNames$Phylum == "Bacillus", ]
colnames(unique_niceNames)[2] <- "ID"
#write.table(perSampleNiceNames, file = here::here("intermediate/perSampleNiceNames.tsv"), sep = "\t", row.names = F)

perSampleRelAbunPlot <- plot_grid(perSampleAbunPlots[["gen_raw_rel"]], 
          perSampleAbunPlots[["gen_decont_rel"]],
          perSampleAbunPlots[["gen_decontLotus_rel"]], 
          perSampleAbunPlots[["gen_rawDecontam_rel"]],
          nrow = 1)

perSampleAbsAbunPlot <- plot_grid(perSampleAbunPlots[["gen_raw_abs"]], 
          perSampleAbunPlots[["gen_decont_abs"]],
          perSampleAbunPlots[["gen_decontLotus_abs"]], 
          perSampleAbunPlots[["gen_rawDecontam_abs"]],
          nrow = 1)

# ggsaveAll(
#   here::here("figures/"),
#   filename = paste0("perSampleAbsAbunPlot"),
#   plot = perSampleAbsAbunPlot,
#   devices = c("svg", "pdf", "png"),
#   height = 8,
#   width = 13
# )
# 
# ggsaveAll(
#   here::here("figures/"),
#   filename = paste0("perSampleRelAbunPlot"),
#   plot = perSampleRelAbunPlot,
#   devices = c("svg", "pdf", "png"),
#   height = 8,
#   width = 13
# )
# 

merged_lotus <- plot_grid(perSampleAbunPlots[["gen_decontLotus_abs"]] + 
                            xlab("") + 
                            ylab("absolute") +
                            theme(axis.line.x = element_blank(),
                                  axis.ticks.length.x = unit(0, "mm"),
                                  axis.text.x = element_text(size = 8),
                                  plot.margin = margin(0, 0, -9, 0, "pt")
                                  ), 
                          perSampleAbunPlots[["gen_decontLotus_rel"]] + 
                            xlab("") + 
                            ylab("absolute") +
                            theme(axis.line.x = element_blank(),
                                  axis.text.x = element_blank(),
                                  plot.margin = margin(-5, 0, -20, 0, "pt")
                                  ), 
                          nrow = 2, 
                          rel_heights = c(1.3, 1))
merged_lotus
merged_lotus2 <- plot_grid(merged_lotus, legendPlot, nrow = 1, rel_widths = c(1.3, 0.4))
merged_lotus2
# safe only our recommended-approach plot
ggsaveAll(
  here::here("figures/"),
  filename = paste0("perSampleMergedAbunPlot_lotusOnly"),
  plot = merged_lotus,
  devices = c("svg", "pdf", "png"),
  height = 6,
  width = 7
)

merged_raw <- plot_grid(perSampleAbunPlots[["gen_raw_abs"]] + 
                            xlab("") + 
                            ylab("absolute") +
                            theme(axis.line.x = element_blank(),
                                  axis.ticks.length.x = unit(0, "mm"),
                                  axis.text.x = element_text(size = 8),
                                  plot.margin = margin(0, 0, -9, 0, "pt")
                            ), 
                          perSampleAbunPlots[["gen_raw_rel"]] + 
                            xlab("") + 
                            ylab("absolute") +
                            theme(axis.line.x = element_blank(),
                                  axis.text.x = element_blank(),
                                  plot.margin = margin(-5, 0, -20, 0, "pt")
                            ), 
                          nrow = 2,
                        rel_heights = c(1.3, 1))
merged_raw
merged_raw2 <- plot_grid(merged_raw, legendPlot, nrow = 1, rel_widths = c(1.3, 0.4))
merged_raw2
# safe only our recommended-approach plot
ggsaveAll(
  here::here("figures/"),
  filename = paste0("perSampleMergedAbunPlot_raw"),
  plot = merged_raw,
  devices = c("svg", "pdf", "png"),
  height = 6,
  width = 7
)





# 
# ggsaveAll(
#   here::here("figures/"),
#   filename = paste0("perSampleAbsAbunPlot_lotusOnly"),
#   plot = perSampleAbunPlots[["gen_decontLotus_abs"]] + xlab(""),
#   devices = c("svg", "pdf", "png"),
#   height = 3.5,
#   width = 6
# )
# 
# ggsaveAll(
#   here::here("figures/"),
#   filename = paste0("genusLegendHorizontal"),
#   plot = legendPlotHor,
#   devices = c("svg", "pdf", "png"),
#   height = 1,
#   width = 6
# )


# rel_PhylumAllRat <-
#   AbundanceAll[, keppGen] #nur die bact aus dem standard und den kontrollen
# 
# rel_PhylumAllRat <- as.data.frame(t(rel_PhylumAllRat))
# rel_PhylumAllRat$Phylum <- row.names(rel_PhylumAllRat)
# 
# for (i in 1:(ncol(rel_PhylumAllRat)-1)) {
#   rel_PhylumAllRat[, i] <- as.numeric(rel_PhylumAllRat[, i])
# } # macht aus den abundance für die bacteria numerical zahlen
# 
# molten_rel_PhylumAllRat <- melt(data = rel_PhylumAllRat, id.vars = "Phylum")
# molten_rel_PhylumAllRat$labels <- scales::percent(molten_rel_PhylumAllRat$value)
# molten_rel_PhylumAllRat$variable <-
#   str_replace((molten_rel_PhylumAllRat$variable), "Abundance_", "")
# 
# relAbunIndi <-
#   ggplot (data = molten_rel_PhylumAllRat, aes(x = variable, y = value, fill = Phylum)) +
#   geom_bar(stat = "identity", color = "black") +
#   scale_fill_manual(values=cbPalette) +
#   theme_classic() +
#   theme(#axis.text.x = element_blank(),
#     legend.text = element_text(face = "italic"),
#     axis.title.y = element_blank(),
#     axis.line.y = element_blank(),
#     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#     axis.text.y = element_blank(),
#     axis.ticks = element_blank()) +
#   xlab("sample") +
#   labs(fill = "genus")
# relAbunIndi
# ggsaveAll(
#   here::here("figures/"),
#   filename = paste0("relAbundundancePerSampleRatio99_", dataType),
#   plot = relAbunIndi,
#   devices = c("svg", "pdf"),
#   height = 6,
#   width = 4
# )


#####                                                                #####
##########                                                      ##########
###############                 abs abund looped           ###############
##########                                                      ##########
#####                                                                #####

absAbunPlots <- list()
for (k in c("group", "ConcGroup", "ratio")) {
  print(k)
  if (k == "group") {
    next
  }
  absAbunPlots[[k]] <- list()
  
  for (AbAll in names(setsAbundanceAll)) {
     #AbAll <- "gen_rawDecontam_abs"
    if (grepl("_rel", AbAll)) {
      next
    }
    print(AbAll)
    AbundanceAll <- setsAbundanceAll[[AbAll]]
    AbundanceAll <- AbundanceAll[rownames(AbundanceAll) != "rel", ]
    if (AbAll == "gen_decont_abs") {
      print("association of cutibacterium readCount to sample dilution.")
      plot(AbundanceAll$Concentration[AbundanceAll$ratio == "1:99"], 
           AbundanceAll$Cutibacterium[AbundanceAll$ratio == "1:99"])
    }
    abundance_groups <- list()
    for (i in unique(AbundanceAll[[k]])) {
      print(i)
      abundance_groups[[i]] <-
        colSums(AbundanceAll[AbundanceAll[[k]] == i, !(colnames(AbundanceAll) %in% non_bactNames)])
    }
    
    AbundanceAll_merged <-
      as.data.frame(bind_rows(abundance_groups, .id = "column_label"))
    rownames(AbundanceAll_merged) <- AbundanceAll_merged$column_label
    AbundanceAll_merged$column_label <- NULL
    
    rel_PhylumAll <-
      t(AbundanceAll_merged[, keppGen])
    
    rel_PhylumAll <- as.data.frame(rel_PhylumAll)
    rel_PhylumAll$Phylum <- row.names(rel_PhylumAll)
    
    for (i in 1:(ncol(rel_PhylumAll) - 1)) {
      print(i)
      rel_PhylumAll[, i] <- as.numeric(rel_PhylumAll[, i])
    } # macht aus den abundance für die bacteria numerical zahlen
    
    molten_rel_PhylumAll <-
      melt(data = rel_PhylumAll, id.vars = "Phylum")
    # molten_rel_PhylumAll$labels <-
    #   scales::percent(molten_rel_PhylumAll$value)
    molten_rel_PhylumAll$variable <-
      str_replace((molten_rel_PhylumAll$variable), "Abundance_", "")
    if (k == "ratio") {
      molten_rel_PhylumAll$variable <-
        factor(
          molten_rel_PhylumAll$variable,
          levels = c(
            "pure_bact",
            "99:1",
            "90:10",
            "75:25",
            "50:50",
            "25:75",
            "10:90",
            "1:99",
            "pure_human",
            "empty"
          )
        )
    } else if (k == "ConcGroup") {
      molten_rel_PhylumAll$variable <-
        factor(
          molten_rel_PhylumAll$variable,
          levels = c(
            "pure_bact",
            "01.3",
            "02.5",
            "05.2",
            "09.8",
            "19.9",
            "40.7",
            "pure_human",
            "empty"
          )
        )
    }
    
    # xlabText <-
    #   switch(k,
    #          "group" = "ratio group as assigned by Theda",
    #          "ConcGroup" = "abs. concentration group, 1 is lowest",
    #          "ratio" = "ratio human:bacterial")
    xlabText <- sub(pattern = "gen_", replacement = "", x = AbAll)
    xlabText <- sub(pattern = "_abs", replacement = "", x = xlabText)
    
    relAbunGroup <-
      ggplot (data = molten_rel_PhylumAll, aes(x = variable, y = value, fill = Phylum)) +
      geom_bar(stat = "identity", color = "black") +
      scale_fill_manual(values = cbPalette) +
      theme_classic() +
      theme(
        #axis.text.x = element_blank(),
        legend.text = element_text(face = "italic"),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1
        ),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank()
      ) +
      xlab(paste0(xlabText)) +
      labs(fill = "genus")
    absAbunPlots[[k]][[AbAll]] <- relAbunGroup
  }
}
# 
# mergedPlots_abs_group <- plot_grid(absAbunPlots[["group"]][["gen_raw_abs"]], 
#                                    absAbunPlots[["group"]][["gen_decont_abs"]],
#                                    absAbunPlots[["group"]][["gen_decontLotus_abs"]], ncol = 1)
# mergedPlots_abs_ConcGroup <- plot_grid(absAbunPlots[["ConcGroup"]][["gen_raw_abs"]], 
#                                        absAbunPlots[["ConcGroup"]][["gen_decont_abs"]],
#                                        absAbunPlots[["ConcGroup"]][["gen_decontLotus_abs"]], ncol = 1)
# mergedPlots_abs_ratio <- plot_grid(absAbunPlots[["ratio"]][["gen_raw_abs"]], 
#                                    absAbunPlots[["ratio"]][["gen_decont_abs"]],
#                                    absAbunPlots[["ratio"]][["gen_decontLotus_abs"]], ncol = 1)
# 
# mergedPlots_abs <- plot_grid(mergedPlots_abs_group, 
#                              mergedPlots_abs_ConcGroup,
#                              mergedPlots_abs_ratio, ncol = 3)

# extract legend for later plotting in plot_grid
absAbunPlots[["ConcGroup"]][["gen_raw_abs"]] <- absAbunPlots[["ConcGroup"]][["gen_raw_abs"]] +
  theme(legend.background = element_rect(fill = 'transparent'))
legendPlot <- as_ggplot(get_legend(absAbunPlots[["ConcGroup"]][["gen_raw_abs"]], position = "right"))
legendPlot <-
  legendPlot + theme(
    plot.background = element_rect(fill = 'transparent'),
    panel.background = element_rect(fill = 'transparent'),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  )

# remove legend from all the plots
absAbunPlotsNL <- list()
for (i in names(absAbunPlots)) {
  for (j in names(absAbunPlots$ConcGroup)) {
    absAbunPlotsNL[[i]][[j]] <- absAbunPlots[[i]][[j]] + theme(legend.position = "none")
  }
}
titleConcGroup <- ggdraw() + 
  draw_label(
    "absolute genera abundance per\ndilution group (mean absolute DNA concentration in ng/μL)",
    angle = 0,
  ) 
mergedPlots_abs_ConcGroup <- plot_grid(absAbunPlotsNL[["ConcGroup"]][["gen_raw_abs"]], 
                                       absAbunPlotsNL[["ConcGroup"]][["gen_decont_abs"]],
                                       absAbunPlotsNL[["ConcGroup"]][["gen_decontLotus_abs"]], 
                                       absAbunPlotsNL[["ConcGroup"]][["gen_rawDecontam_abs"]],
                                       nrow = 1)
mergedPlots_abs_ConcGroup <- plot_grid(titleConcGroup, mergedPlots_abs_ConcGroup,
                                       nrow = 2,
                                       rel_heights = c(0.1, 1))


titleRatio <- ggdraw() + 
  draw_label(
    "absolute genera abundance per\nbacterial to human DNA ratio group",
    angle = 0,
  )
mergedPlots_abs_ratio <- plot_grid(absAbunPlotsNL[["ratio"]][["gen_raw_abs"]], 
                                   absAbunPlotsNL[["ratio"]][["gen_decont_abs"]],
                                   absAbunPlotsNL[["ratio"]][["gen_decontLotus_abs"]],
                                   absAbunPlotsNL[["ratio"]][["gen_rawDecontam_abs"]],
                                   nrow = 1)
mergedPlots_abs_ratio <- plot_grid(titleRatio, 
                                   mergedPlots_abs_ratio, 
                                   nrow = 2,
                                   rel_heights = c(0.1, 1))




mergedPlots_abs <- plot_grid(mergedPlots_abs_ConcGroup,
                             mergedPlots_abs_ratio, nrow = 2)
mergedPlots_abs <- plot_grid(mergedPlots_abs, legendPlot, nrow = 1, rel_widths = c(1,0.25))
mergedPlots_abs

# 
# ggsaveAll(
#   here::here("figures/"),
#   filename = paste0("absAbund_2by3"),
#   plot = mergedPlots_abs,
#   devices = c("svg", "pdf"),
#   height = 15,
#   width = 15
# )


# 2024_04_02 final rel abund figures

plot_grid(relAbunPlotsNL[["ratio"]][["gen_raw_rel"]], 
          relAbunPlotsNL[["ratio"]][["gen_decont_rel"]],
          relAbunPlotsNL[["ratio"]][["gen_decontLotus_rel"]], nrow = 1)

 plot_grid(absAbunPlotsNL[["ratio"]][["gen_raw_abs"]], 
                                   absAbunPlotsNL[["ratio"]][["gen_decont_abs"]],
                                   absAbunPlotsNL[["ratio"]][["gen_decontLotus_abs"]], nrow = 1)

merged_per_ratio_abs_rel <- list()
 for (cleaningMeth in c("raw", "decont", "decontLotus", "rawDecontam")) {
   merged_per_ratio_abs_rel[[cleaningMeth]] <- plot_grid(absAbunPlotsNL[["ratio"]][[paste0("gen_", cleaningMeth, "_abs")]] + 
               xlab("") + 
               ylab("absolute") +
               theme(axis.line.x = element_blank(),
                     axis.ticks.length.x = unit(0, "mm"),
                     axis.text.x = element_text(size = 12),
                     plot.margin = margin(0, 0, -9, 0, "pt")
               ), 
             relAbunPlotsNL[["ratio"]][[paste0("gen_", cleaningMeth, "_rel")]] + 
               xlab("") + 
               ylab("absolute") +
               theme(axis.line.x = element_blank(),
                     axis.text.x = element_blank(),
                     plot.margin = margin(-5, 0, -20, 0, "pt")
               ), 
             nrow = 2, 
             rel_heights = c(1.3, 1))
 }

merged_merged_per_ratio_abs_rel <-
  plot_grid(
    merged_per_ratio_abs_rel[[1]],
    merged_per_ratio_abs_rel[[2]],
    merged_per_ratio_abs_rel[[3]],
    nrow = 1,
    labels = c("naive", "filtered pre-clustering", "filtered post-clustering")
  )
merged_merged_per_ratio_abs_rel

ggsaveAll(
  here::here("figures/"),
  filename = "abs_rel_merged_3ways",
  plot = merged_merged_per_ratio_abs_rel,
  devices = c("svg", "pdf"),
  height = 7,
  width = 15
)



ggsaveAll(
  here::here("figures/"),
  filename = "suppl_abs_rel_decontam",
  plot =  plot_grid(relAbunPlotsNL[["ratio"]][["gen_raw_rel"]] + xlab("naive"), 
                    relAbunPlotsNL[["ratio"]][["gen_decont_rel"]] + xlab("filtered pre-clustering"), 
                    relAbunPlotsNL[["ratio"]][["gen_decontLotus_rel"]] + xlab("filtered post-clustering"), 
                    relAbunPlotsNL[["ratio"]][["gen_rawDecontam_rel"]] + xlab("Decontam filtered"), 
                    nrow = 1
  ),
  devices = c("svg", "pdf"),
  height = 7,
  width = 10
)




plot_grid(relAbunPlotsNL[["ConcGroup"]][["gen_raw_rel"]], 
          relAbunPlotsNL[["ConcGroup"]][["gen_decont_rel"]],
          relAbunPlotsNL[["ConcGroup"]][["gen_decontLotus_rel"]], nrow = 1)

plot_grid(absAbunPlotsNL[["ConcGroup"]][["gen_raw_abs"]], 
          absAbunPlotsNL[["ConcGroup"]][["gen_decont_abs"]],
          absAbunPlotsNL[["ConcGroup"]][["gen_decontLotus_abs"]], nrow = 1)

merged_per_Conc_abs_rel <- list()
for (cleaningMeth in c("raw", "decont", "decontLotus", "rawDecontam")) {
  merged_per_Conc_abs_rel[[cleaningMeth]] <- plot_grid(absAbunPlotsNL[["ConcGroup"]][[paste0("gen_", cleaningMeth, "_abs")]] + 
                                                          xlab("") + 
                                                          ylab("absolute") +
                                                          theme(axis.line.x = element_blank(),
                                                                axis.ticks.length.x = unit(0, "mm"),
                                                                axis.text.x = element_text(size = 12),
                                                                plot.margin = margin(0, 0, -9, 0, "pt")
                                                          ), 
                                                        relAbunPlotsNL[["ConcGroup"]][[paste0("gen_", cleaningMeth, "_rel")]] + 
                                                          xlab("") + 
                                                          ylab("absolute") +
                                                          theme(axis.line.x = element_blank(),
                                                                axis.text.x = element_blank(),
                                                                plot.margin = margin(-5, 0, -20, 0, "pt")
                                                          ), 
                                                        nrow = 2, 
                                                        rel_heights = c(1.3, 1))
}

merged_merged_per_Conc_abs_rel <-
  plot_grid(
    merged_per_ratio_abs_rel[[1]],
    merged_per_ratio_abs_rel[[2]],
    merged_per_ratio_abs_rel[[3]],
    nrow = 1,
    labels = c("naive", "filtered pre-clustering", "filtered post-clustering")
  )
merged_merged_per_Conc_abs_rel




ggsaveAll(
  here::here("figures/"),
  filename = "abs_rel__concGroup_merged_3ways",
  plot = merged_merged_per_Conc_abs_rel,
  devices = c("svg", "pdf"),
  height = 7,
  width = 15
)


#####                                                                #####
##########                                                      ##########
###############                   abs per sample              ###############
##########                                                      ##########
#####                                                                #####

warning("subset here manually to only ratio 99 etc.")
rel_PhylumAllRat <-
  AbundanceAll[AbundanceAll$ratio == "99", keppGen] #nur die bact aus dem standard und den kontrollen
rel_PhylumAllRat <-
  AbundanceAll[, keppGen]

rel_PhylumAllRat <- as.data.frame(t(rel_PhylumAllRat))
rel_PhylumAllRat$Phylum <- row.names(rel_PhylumAllRat)

for (i in 1:(ncol(rel_PhylumAllRat)-1)) {
  rel_PhylumAllRat[, i] <- as.numeric(rel_PhylumAllRat[, i])
} # macht aus den abundance für die bacteria numerical zahlen

molten_rel_PhylumAllRat <- melt(data = rel_PhylumAllRat, id.vars = "Phylum")
# molten_rel_PhylumAllRat$labels <- scales::percent(molten_rel_PhylumAllRat$value)
molten_rel_PhylumAllRat$variable <-
  str_replace((molten_rel_PhylumAllRat$variable), "Abundance_", "")

relAbunIndi <-
  ggplot (data = molten_rel_PhylumAllRat, aes(x = variable, y = value, fill = Phylum)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values=cbPalette) +
  theme_classic() +
  theme(#axis.text.x = element_blank(),
    legend.text = element_text(face = "italic"),
    axis.title.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()) +
  xlab("sample") +
  labs(fill = "genus")
relAbunIndi
# ggsaveAll(
#   here::here("figures/"),
#   filename = paste0("relAbundundancePerSampleRatio99_", dataType),
#   plot = relAbunIndi,
#   devices = c("svg", "pdf"), 
#   height = 6, 
#   width = 4
# )



#####                                                                #####
##########                                                      ##########
###############                    betadiv                 ###############
##########                                                      ##########
#####                                                                #####

beta_input <- setsAbundanceAll$gen_raw_rel
beta_input$ID <- rownames(beta_input)
rownames(beta_input) <- paste0(rownames(beta_input), "_raw") 
beta_input$processing <- "raw"
for (i in names(setsAbundanceAll)) {
  if (grepl("_abs", i) | i == "gen_raw_rel") {
    next
  }
  tempor <- setsAbundanceAll[[i]]
  processing <- stringr::str_match(pattern = "_[a-z].+_", string = i)
  processing <- str_remove_all(string = processing, pattern = "_")
  tempor$ID <- rownames(tempor)
  rownames(tempor) <- paste0(rownames(tempor), "_", processing) 
  tempor$processing <- processing
  beta_input <- rbind(beta_input, tempor)
}

beta_meta <- beta_input

warning("setting NAs to zeros (important for theoretic sample!)")
#beta_input <- beta_input[!(beta_input$ConcGroup %in% c("empty", "pure_human", "theoretical")), ]
#beta_input <- beta_input[!(beta_input$ConcGroup %in% c("empty", "pure_human")), ]
beta_input <- beta_input[!(beta_input$ConcGroup %in% c("empty", "theoretical")), ]
beta_input <- beta_input[, !(colnames(beta_input) %in% c(non_bactNames, "processing", "ID"))]
beta_input[is.na(beta_input)] <- 0

# collect complete set for later analysis
bray_distAll <- vegdist(beta_input, na.rm = T)
# dist of theoretical to pure_bact
as.matrix(bray_distAll)["rel_raw", "P01.H01.8_raw"]
# subset to only _raw samples
bray_dist_allMatrix <- as.matrix(bray_distAll)[1:44, 1:44]
# standard deviation of pure-bact vs. everything else
sd(bray_dist_allMatrix["P01.H01.8_raw", -c(39, 44)])
# is theoretical to pure bact > als pure_bact to everythig else?
t.test(x = bray_dist_allMatrix["P01.H01.8_raw", -c(39, 44)], as.matrix(bray_distAll)["rel_raw", "P01.H01.8_raw"])
# subset for raw or lotus-cleaned samples
beta_inputRaw <- beta_input[grepl("_raw$", rownames(beta_input)), ]
bray_distRaw <- vegdist(beta_inputRaw, na.rm = T)
pcoaEHRaw <- as.data.frame(cmdscale (bray_distRaw, k = 2))
pcoaEHRaw <- merge(pcoaEHRaw, beta_meta[, c(non_bactNames, "processing", "ID")], by = 0)
rownames(pcoaEHRaw) <- pcoaEHRaw$Row.names
pcoaEHRaw$Row.names <- NULL
pcoaEHRaw <- pcoaEHRaw[!(pcoaEHRaw$ConcGroup %in% c("empty", "pure_human")), ]
pcoaEHRaw <- merge(pcoaEHRaw, unique_niceNames[, c("ID", "sampleType")])
colnames(pcoaEHRaw)[colnames(pcoaEHRaw) == "sampleType"] <- "label"
centroidsRaw <- aggregate(cbind(V1,V2)~ratio,pcoaEHRaw,mean)
centroidsRaw$label <- centroidsRaw$ratio
centroidsRawConc <- aggregate(cbind(V1,V2)~ConcGroup,pcoaEHRaw,mean)[1:6, ]
centroidsRawConc$label <- centroidsRawConc$ConcGroup
colnames(centroidsRawConc)[1] <- "ratio"

#pcoaEHRaw$label <- NULL

beta_inputLotus <- beta_input[grepl("_decontLotus$", rownames(beta_input)), ]
bray_distLotus <- vegdist(beta_inputLotus, na.rm = T)
pcoaEHLotus <- as.data.frame(cmdscale (bray_distLotus, k = 2))
pcoaEHLotus <- merge(pcoaEHLotus, beta_meta[, c(non_bactNames, "processing", "ID")], by = 0)
rownames(pcoaEHLotus) <- pcoaEHLotus$Row.names
pcoaEHLotus$Row.names <- NULL
pcoaEHLotus <- pcoaEHLotus[!(pcoaEHLotus$ConcGroup %in% c("empty", "pure_human")), ]
pcoaEHLotus <- merge(pcoaEHLotus, unique_niceNames[, c("ID", "sampleType")])
colnames(pcoaEHLotus)[colnames(pcoaEHLotus) == "sampleType"] <- "label"
centroidsLotus <- aggregate(cbind(V1,V2)~ratio,pcoaEHLotus,mean)
centroidsLotus$label <- centroidsLotus$ratio


pcoaGenRaw <-
  ggplot (
    data = pcoaEHRaw,
    aes (
      x = V1,
      y = V2,
      fill = as.factor(ratio)#,
     # label = as.factor(label)
    )
  ) +
  theme_classic () +
  geom_point (
    size = 3,
    stroke = 0.8,
    pch = 21,
    #alpha = 0.6,
    aes(alpha = as.factor(ConcGroup))
  ) +
  # geom_point(data = centroidsRaw, size = 5, shape = 16, color = "black") + # centroides hinzufügen
  # geom_point(data = centroidsRaw, size = 4, shape = 16) +
  # geom_label_repel(data = centroidsRaw, force = 5, max.overlaps = 45, aes(label = as.factor(label))) +
  # geom_line(col = "gray") +
  xlab ("PCo 1") +
  ylab ("PCo 2") +
  #stat_ellipse(level = 0.8) +

   scale_fill_manual(values = cbPalette) +
  # geom_point(data = centroidsRawConc, size = 5, shape = 17, color = "black") +
  # geom_label_repel(data = centroidsRawConc, force = 5, max.overlaps = 45, aes(label = as.factor(label)), color = "black") +
  
  #geom_point(data = centroidsRawConc, size = 4, shape = 17) +
  theme (
    axis.title.x = element_text (size = 13),
    axis.text.x = element_text (size = 13),
    axis.text.y = element_text (size = 13),
    axis.title.y = element_text (size = 13),
    legend.position = "bottom"
  ) +
  #labs(color = "concentration group")
  labs(color = "ratio") +
  guides(alpha = "none")
print(pcoaGenRaw)

pcoaGenLotus <-
  ggplot (
    data = pcoaEHLotus,
    aes (
      x = V1,
      y = V2,
      color = as.factor(ratio)#,
      # label = as.factor(label)
    )
  ) +
  theme_classic () +
  geom_point (
    size = 3,
    stroke = 1.3,
    alpha = 0.8
  ) +
  geom_point(data = centroidsLotus, size = 5, shape = 16, color = "black") + # centroides hinzufügen
  geom_point(data = centroidsLotus, size = 4, shape = 16) +
  geom_label_repel(data = centroidsLotus, force = 5, max.overlaps = 45, aes(label = as.factor(label))) +
  # geom_line(col = "gray") +
  xlab ("PCo 1") +
  ylab ("PCo 2") +
  stat_ellipse(level = 0.8) +
  
  scale_color_manual(values = cbPalette) +
  theme (
    axis.title.x = element_text (size = 13),
    axis.text.x = element_text (size = 13),
    axis.text.y = element_text (size = 13),
    axis.title.y = element_text (size = 13),
    legend.position = "bottom"
  ) +
  #labs(color = "concentration group")
  labs(color = "ratio")
print(pcoaGenLotus)

# repeat for both sets combined, try to link sample pairs
beta_inputBoth <- beta_input[rownames(beta_input) %in% c(rownames(beta_inputLotus), rownames(beta_inputRaw)), ]
bray_distBoth <- vegdist(beta_inputBoth, na.rm = T)
pcoaEHBoth <- as.data.frame(cmdscale (bray_distBoth, k = 2))
pcoaEHBoth <- merge(pcoaEHBoth, beta_meta[, c(non_bactNames, "processing", "ID")], by = 0)
rownames(pcoaEHBoth) <- pcoaEHBoth$Row.names
pcoaEHBoth$Row.names <- NULL
pcoaEHBoth <- pcoaEHBoth[!(pcoaEHBoth$ConcGroup %in% c("empty", "pure_human")), ]
pcoaEHBoth <- merge(pcoaEHBoth, unique_niceNames[, c("ID", "sampleType")])
colnames(pcoaEHBoth)[colnames(pcoaEHBoth) == "sampleType"] <- "label"
centroidsBoth <- aggregate(cbind(V1,V2)~ratio + processing,pcoaEHBoth,mean)
centroidsBoth$label <- paste0(centroidsBoth$ratio, "_", centroidsBoth$processing)

pcoaGenBoth <-
  ggplot (
    data = pcoaEHBoth,
    aes (
      x = V1,
      y = V2,
      color = as.factor(ratio)#,
      # label = as.factor(label)
    )
  ) +
  theme_classic () +
  # geom_point (
  #   size = 3,
  #   stroke = 1.3,
  #   alpha = 0.1
  # ) +
  #geom_line(aes(group = label)) +
  geom_point(data = centroidsBoth, size = 5, shape = 16, color = "black") + # centroides hinzufügen
  geom_point(data = centroidsBoth, size = 4, shape = 16) +
  geom_label_repel(data = centroidsBoth, force = 5, max.overlaps = 45, aes(label = as.factor(label)), alpha = 0.8, size = 3) +
  geom_line(data = centroidsBoth, aes(group = ratio)) +
  # geom_line(col = "gray") +
  xlab ("PCo 1") +
  ylab ("PCo 2") +
  #stat_ellipse(level = 0.8) +
  scale_color_manual(values = cbPalette) +
  theme (
    axis.title.x = element_text (size = 13),
    axis.text.x = element_text (size = 13),
    axis.text.y = element_text (size = 13),
    axis.title.y = element_text (size = 13),
    legend.position = "bottom"
  ) +
  #labs(color = "concentration group")
  labs(color = "ratio")
print(pcoaGenBoth)





dummyMetadec <- metadeconfoundR::MetaDeconfound(reduced_feature, metaMatMetformin, nnodes = 3)
dummyMetadecLegend <- BuildHeatmap(dummyMetadec, d_col = c("white", "grey", "black"), d_range = "full")



warning("check which file to save to!")
ggsaveAll(
  here::here("figures/"),
  filename = "betaDiv_RelGenRaw_Conc_vs_ratio",
  plot = plot_grid(pcoaGenRaw, dummyMetadecLegend),
  devices = c("svg", "pdf"),
  height = 5,
  width = 10
)
# ggsaveAll(
#   here::here("figures/"),
#   filename = "betaDiv_RelGenLotusEllipseCentroids",
#   plot = pcoaGenLotus,
#   devices = c("svg", "pdf"), 
#   height = 10, 
#   width = 12
# )
# 
ggsaveAll(
  here::here("figures/"),
  filename = "betaDiv_RelGenBothLinkedCentroids",
  plot = pcoaGenBoth,
  devices = c("svg", "pdf"),
  height = 6,
  width = 6
)



# ggsaveAll(
#   here::here("figures/"),
#   filename = paste0("betadivZoomed_", dataType),
#   plot = pcoaGen,
#   devices = c("svg", "pdf"), 
#   height = 6, 
#   width = 6
# )
# 
# ggsaveAll(
#   here::here("figures/"),
#   filename = paste0("betadivRat9_99_", dataType),
#   plot = pcoaGen,
#   devices = c("svg", "pdf"), 
#   height = 6, 
#   width = 6
# )
# 
# ggsaveAll(
#   here::here("figures/"),
#   filename = paste0("betadivAll_threePre"),
#   plot = pcoaGen,
#   devices = c("svg", "pdf"), 
#   height = 12, 
#   width = 12
# )
# 
# ggsaveAll(
#   here::here("figures/"),
#   filename = "betadivGenLotus_AllSamples",
#   plot = pcoaGen,
#   devices = c("svg", "pdf"), 
#   height = 6, 
#   width = 6
# )
# 
# ggsaveAll(
#   here::here("figures/"),
#   filename = "betadivGenLotus_humanRemoved",
#   plot = pcoaGen,
#   devices = c("svg", "pdf"), 
#   height = 10, 
#   width = 12
# )


# auch noch alles einmal mit absolute abundances
# 


#####                                                                #####
##########                                                      ##########
###############                    betadiv on OTU               ###############
##########                                                      ##########
#####                                                                #####
 # this needs to be reworked completely
# # read in OTUs and redo analysis
# beta_input <- setsAbundanceAll$gen_raw_rel
# rownames(beta_input) <- paste0(rownames(beta_input), "_raw") 
# beta_input$processing <- "raw"
# for (i in names(setsAbundanceAll)) {
#   if (grepl("_abs", i) | i == "gen_raw_rel") {
#     next
#   }
#   tempor <- setsAbundanceAll[[i]]
#   processing <- stringr::str_match(pattern = "_[a-z].+_", string = i)
#   processing <- str_remove_all(string = processing, pattern = "_")
#   rownames(tempor) <- paste0(rownames(tempor), "_", processing) 
#   tempor$processing <- processing
#   beta_input <- rbind(beta_input, tempor)
# }
# 
# beta_meta <- beta_input
# 
# warning("setting NAs to zeros (important for theoretic sample!)")
# beta_input <- beta_input[, !(colnames(beta_input) %in% c(non_bactNames, "processing"))]
# beta_input[is.na(beta_input)] <- 0
# bray_dist <- vegdist(beta_input, na.rm = T)
# bray_distAll <- bray_dist
# pcoaEH <- as.data.frame(cmdscale (bray_dist, k = 2))
# 
# # apply congruence analysis
# 
# 
# pcoaEH <- merge(pcoaEH, beta_meta[, c(non_bactNames, "processing", "ID")], by = 0)
# 
# 
# rownames(pcoaEH) <- pcoaEH$Row.names
# pcoaEH$Row.names <- NULL
# pcoaEH <- pcoaEH[!(pcoaEH$ConcGroup %in% c("empty", "pure_human")), ]
# pcoaEH$label <- as.character(paste0("C:", pcoaEH$ConcGroup, " ", "R:", pcoaEH$ratio, " ", "P:", pcoaEH$processing))
# 
# pcoaeAll <- pcoaEH
# pcoaEH <- pcoaEH[(pcoaEH$ratio %in% c("99:1", "90:10")), ]
# 
# 
# 
# pcoaGen <-
#   ggplot (
#     pcoaEH[pcoaEH$processing == "decont", ],
#     aes (
#       x = V1,
#       y = V2,
#       color = as.character(ConcGroup),
#       label = label,
#       #group = ConcGroup
#     )
#   ) +
#   theme_classic () +
#   geom_point (
#     aes (color = as.character(ConcGroup)),
#     size = 3,
#     stroke = 1.3,
#     alpha = 0.8
#   ) +
#   geom_label_repel(force = 5, max.overlaps = 45) +
#   # geom_line(col = "gray") +
#   xlab ("PCo 1") +
#   ylab ("PCo 2") +
#   scale_color_manual(values = cbPalette) +
#   theme (
#     axis.title.x = element_text (size = 13),
#     axis.text.x = element_text (size = 13),
#     axis.text.y = element_text (size = 13),
#     axis.title.y = element_text (size = 13),
#     legend.position = "bottom"
#   ) +
#   labs(color = "concentration group")
# print(pcoaGen)
# 
# warning("check wich file to save to!")
# ggsaveAll(
#   here::here("figures/"),
#   filename = paste0("betadivAll_", dataType),
#   plot = pcoaGen,
#   devices = c("svg", "pdf"), 
#   height = 6, 
#   width = 6
# )
# 
# ggsaveAll(
#   here::here("figures/"),
#   filename = paste0("betadivZoomed_", dataType),
#   plot = pcoaGen,
#   devices = c("svg", "pdf"), 
#   height = 6, 
#   width = 6
# )
# 
# ggsaveAll(
#   here::here("figures/"),
#   filename = paste0("betadivRat9_99_", dataType),
#   plot = pcoaGen,
#   devices = c("svg", "pdf"), 
#   height = 6, 
#   width = 6
# )
# 
# ggsaveAll(
#   here::here("figures/"),
#   filename = paste0("betadivAll_threePre"),
#   plot = pcoaGen,
#   devices = c("svg", "pdf"), 
#   height = 12, 
#   width = 12
# )


# auch noch alles einmal mit absolute abundances
# 
# permanova tsting 

#bray_dist <- vegdist(beta_input[rownames(pcoaEH)[pcoaEH$processing %in% c("raw", "decontLotus")], ], na.rm = T)
#bray_dist <- bray_distBoth
bray_distBothDF <- as.data.frame(as.matrix(bray_distBoth))
rownames(pcoaEHBoth) <- paste0(pcoaEHBoth$ID, "_", pcoaEHBoth$processing)
#permanova <- adonis2(bray_dist ~ as.factor(human) + as.factor(ConcGroup) + as.factor(processing), strata = pcoaEH$ID,data = pcoaEH[pcoaEH$processing %in% c("raw", "decontLotus"), ])
permanova <- adonis2(bray_distBoth ~ human + ratio + Concentration + processing + ID, data = pcoaEHBoth)
permanova
#permanova <- adonis2(brayDistDF[, -(ncol(brayDistDF))] ~ as.factor(human) + as.factor(ConcGroup) + as.factor(processing) + as.factor(ID),data = pcoaEH)


#bray_dist <- vegdist(beta_input[rownames(pcoaEH)[pcoaEH$processing == "decontLotus"], ], na.rm = T)
rownames(pcoaEHLotus) <- paste0(pcoaEHLotus$ID, "_", pcoaEHLotus$processing)
adonis2(bray_distLotus ~ as.factor(human) + as.factor(ratio) + ConcGroup + Concentration, data = pcoaEHLotus)

adonis2(bray_distLotus ~ human + Concentration, data = pcoaEHLotus)
adonis2(bray_distLotus ~ ratio , data = pcoaEHLotus)
adonis2(bray_distLotus ~ as.factor(ratio) , data = pcoaEHLotus)
adonis2(bray_distLotus ~ human , data = pcoaEHLotus)



#bray_dist <- vegdist(beta_input[rownames(pcoaEH)[pcoaEH$processing == "raw"], ], na.rm = T)
adonis2(bray_distRaw ~ human + ratio + as.factor(ConcGroup) + as.factor(Concentration), data = pcoaEHRaw)
adonis2(bray_distRaw ~ human + ratio + as.factor(Concentration), data = pcoaEHRaw)
adonis2(bray_distRaw ~ human + ratio + as.factor(ConcGroup), data = pcoaEHRaw)
adonis2(bray_distRaw ~ human + ratio + Concentration, data = pcoaEHRaw)

permanova <- adonis2(bray_distRaw ~ as.factor(bacteria) , data = pcoaEHRaw)

permanova <-  adonis2(bray_distRaw ~ human + ratio + as.factor(ConcGroup), data = pcoaEHRaw)
permanova



#mm <- lmer (data = dist_long, rank (brayDist) ~ same_ID + same_processing + same_ratio + same_ConcGroup + (1 | Sample1) + (1 | Sample2), REML = F)
varianceTable <- as.data.frame(permanova)
# e.g. ID explains 8%, extraction 6%, collection <<1%
varianceTable$VarExplained <- varianceTable$R2
varianceTable$Variable <- rownames(varianceTable)
#varianceTable[nrow(varianceTable) +1, ] <- c(rep(1, 4), (1 - sum(varianceTable$VarExplained)), "Residual")
varianceTable$VarExplained <- as.numeric(varianceTable$VarExplained)
varianceTable$VarLabels <- scales::percent(varianceTable$VarExplained)
print("Variance table (ANOVA output:)")
print(varianceTable)

#remove insignificant variable
varianceTable <- varianceTable[!rownames(varianceTable) == "Total", ]

varianceTable2 <- varianceTable %>% 
  arrange(desc(Variable)) %>%
  mutate(prop = VarExplained / sum(varianceTable$VarExplained) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )
# varianceTable %>%
#   mutate(VarLabels = fct_reorder(VarLabels, VarExplained)) %>%
ggplot (varianceTable2) +
  geom_bar(aes (x = "", y = prop, fill = Variable), width = 1, stat = "identity", color = "black") +
  coord_polar("y", start = 0) +
  theme_void() +
  geom_text(aes(x=1.7, y = ypos, label=VarLabels)) ->
  pVarExpl
pVarExpl



#####                                                                #####
##########                                                      ##########
###############              congruence analysis            ###############
##########                                                      ##########
#####                                                                #####


# warning("this is not finished yet!!!")
# brayDist <- as.matrix(vegdist (t(rtk_use), method = "bray", na.rm = T))
# 
# pcoaE$dog_sampled <- 0
# hum_w_dog <- unique(subset(pcoaE, dog > 0 & species == "human")$study_id)
# all_dogs <- rownames(pcoaE)[grep(pattern = "[1-9]_D", x = rownames(pcoaE))]
# 
# hum_dog_congru <- as.data.frame(hum_w_dog)
# rownames(hum_dog_congru) <- hum_dog_congru$hum_w_dog
# hum_dog_congru$hum_w_dog <- NULL
# hum_dog_congru$own <- NA
# hum_dog_congru$other <- NA
# for (human in hum_w_dog) {
#   print(human)
#   intermed <- subset(pcoaE, study_id == human)
#   if (length(grep(pattern = "[1-9]_D", x = rownames(intermed))) > 0) {
#     pcoaE$dog_sampled[pcoaE$study_id == human] <- 1
#     print(paste0("dog sampled for ", human))
#     humSampls <- rownames(intermed)[grep(pattern = "[1-9]_D", x = rownames(intermed), invert = T)]
#     dogSampls <- rownames(intermed)[grep(pattern = "[1-9]_D", x = rownames(intermed))]
#     other_dogs <- all_dogs[!(all_dogs %in% dogSampls)]
#     # get mean of all human--dog distances for own dog vs. all other dogs
#     own_dist <- c()
#     other_dist <- c()
#     for (sample in humSampls) {
#       own_dist <- c(own_dist, brayDist[sample, dogSampls])
#       other_dist <- c(other_dist, brayDist[sample, other_dogs])
#     }
#     own_dist_mean <- mean(own_dist)
#     other_dist_mean <- mean(other_dist)
#     print(paste0("dist to own dog is ", own_dist_mean))
#     print(paste0("dist to other dogs is ", other_dist_mean))
#     hum_dog_congru[human, "own"] <- own_dist_mean
#     hum_dog_congru[human, "other"] <- other_dist_mean
#   }
# }
# colSums(hum_dog_congru, na.rm = T)

warning("update this, as a pcoaEH object of all sampes is no longer generated up above")
bray_dist <- vegdist(beta_input[rownames(pcoaEH)[pcoaEH$processing %in% c("raw", "decontLotus")], ], na.rm = T)
brayDistDF <- as.data.frame(as.matrix(bray_dist))
brayDistDF <- brayDistDF[rownames(brayDistDF) %in% rownames(pcoaEH[pcoaEH$processing %in% c("raw", "decontLotus"), ]), 
                         colnames(brayDistDF) %in% rownames(pcoaEH[pcoaEH$processing %in% c("raw", "decontLotus"), ])]
brayDistDF$Sample1 <- rownames(brayDistDF)

  dist_long <- pivot_longer(data = brayDistDF, 
                            cols = !Sample1, names_to = "Sample2", values_to = "brayDist")
  #View(dist_long)
  dist_long <- dist_long[dist_long$brayDist != 0, ]
  dist_long$same_ConcGroup <- "NO"
  dist_long$same_ratio <- "NO"
  dist_long$same_processing <- "NO"
  dist_long$same_ID <- "NO"
  
  # dist_long$ageDiff <- NA
  # dist_long$BMIDiff <- NA
  # dist_long$IgEDiff <- NA
  for (i in seq_along(dist_long$Sample1)) {
    s1 <- dist_long$Sample1[i]
    s2 <- dist_long$Sample2[i]
    if (pcoaEH[s1, "ConcGroup"] == pcoaEH[s2, "ConcGroup"]) {
      dist_long$same_ConcGroup[i] <- "YES"
    }
    if (pcoaEH[s1, "ratio"] == pcoaEH[s2, "ratio"]) {
      dist_long$same_ratio[i] <- "YES"
    }
    if (pcoaEH[s1, "processing"] == pcoaEH[s2, "processing"]) {
      dist_long$same_processing[i] <- "YES"
    }
    if (pcoaEH[s1, "ID"] == pcoaEH[s2, "ID"]) {
      dist_long$same_ID[i] <- "YES"
    }
    # if (pcoaEH[s1, "enterotype"] == pcoaEH[s2, "enterotype"]) {
    #   dist_long$same_Enterotype[i] <- "YES"
    # }
    # dist_long$ageDiff[i] <- abs(pcoaEH[s1, "age_in_years"] - pcoaEH[s2, "age_in_years"])
    # dist_long$BMIDiff[i] <- abs(pcoaEH[s1, "bmi"] - pcoaEH[s2, "bmi"])
    # dist_long$IgEDiff[i] <- abs(pcoaEH[s1, "mr_total_ige_1"] - pcoaEH[s2, "mr_total_ige_1"])
    # dist_long$mean_dog[i] <- mean(c(pcoaEH[s1, "dog"], pcoaEH[s2, "dog"]))
  }
 

  
  p_processing <- lrtest(lmer (data = dist_long, rank (brayDist) ~ same_ID + same_processing + same_ratio + same_ConcGroup + (1 | Sample1) + (1 | Sample2), REML = F), 
                         lmer (data = dist_long, rank (brayDist) ~ same_ID + same_ratio + same_ConcGroup + (1 | Sample1) + (1 | Sample2), REML = F))$'Pr(>Chisq)' [2]
  if (p_processing < 2.2e-16) {
    p_processing <- "< 2.2e-16"
  }
  
  p_ratio <- lrtest(lmer (data = dist_long, 
                          rank (brayDist) ~ same_ID + same_processing + same_ratio + same_ConcGroup + (1 | Sample1) + (1 | Sample2), 
                          REML = F), 
                    lmer (data = dist_long, 
                          rank (brayDist) ~ same_ID + same_processing + same_ConcGroup + (1 | Sample1) + (1 | Sample2), 
                          REML = F))$'Pr(>Chisq)' [2]
  
  
  if (p_ratio < 2.2e-16) {
    p_ratio <- "< 2.2e-16"
  }
  p_ConcGroup <- lrtest(lmer (data = dist_long, rank (brayDist) ~ same_ID + same_processing + same_ratio + same_ConcGroup + (1 | Sample1) + (1 | Sample2), REML = F), 
                        lmer (data = dist_long, rank (brayDist) ~ same_ID + same_processing + same_ratio + (1 | Sample1) + (1 | Sample2), REML = F))$'Pr(>Chisq)' [2]
  if (p_ConcGroup < 2.2e-16) {
    p_ConcGroup <- "< 2.2e-16"
  }
  
  p_ID <- lrtest(lmer (data = dist_long, rank (brayDist) ~ same_ID + same_processing + same_ratio + same_ConcGroup + (1 | Sample1) + (1 | Sample2), REML = F), 
                        lmer (data = dist_long, rank (brayDist) ~ same_processing + same_ratio + same_ConcGroup + (1 | Sample1) + (1 | Sample2), REML = F))$'Pr(>Chisq)' [2]
  if (p_ID < 2.2e-16) {
    p_ID <- "< 2.2e-16"
  }

vioPlot_same_ratio <- ggplot(data = dist_long, aes(x = brayDist, y = same_ratio)) + 
  geom_violin(aes(color = same_ratio, fill = same_ratio, alpha = 0.5)) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA, color = "grey60") + 
  labs(y = "same human to bacterial DNA ratio", x = "Bray-Curstis distance between two samples", subtitle = paste0("lmer lrtest p-value: ", formatC(p_ratio, format = "e"))) +
  coord_flip() +
  theme_classic() +
  theme(legend.position = "none")
vioPlot_same_ConcGroup <- ggplot(data = dist_long, aes(x = brayDist, y = same_ConcGroup)) + 
  geom_violin(aes(color = same_ConcGroup, fill = same_ConcGroup, alpha = 0.5)) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA, color = "grey60") + 
  labs(y = "same absolut bacterial DNA concentration", x = "Bray-Curstis distance between two samples", subtitle = paste0("lmer lrtest p-value: ", formatC(p_ConcGroup, format = "e"))) +
  coord_flip() +
  theme_classic() +
  theme(legend.position = "none")
vioPlot_same_processing <- ggplot(data = dist_long, aes(x = brayDist, y = same_processing)) + 
  geom_violin(aes(color = same_processing, fill = same_processing, alpha = 0.5)) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA, color = "grey60") + 
  labs(y = "same human DNA filtering procedure", x = "Bray-Curstis distance between two samples", subtitle = paste0("lmer lrtest p-value: ", formatC(p_processing, format = "e"))) +
  coord_flip() +
  theme_classic() +
  theme(legend.position = "none")
vioPlot_same_ID <- ggplot(data = dist_long, aes(x = brayDist, y = same_ID)) + 
  geom_violin(aes(color = same_ID, fill = same_ID, alpha = 0.5)) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA, color = "grey60") + 
  labs(y = "same ID", x = "Bray-Curstis distance between two samples", subtitle = paste0("lmer lrtest p-value: ", formatC(p_ID, format = "e"))) +
  coord_flip() +
  theme_classic() +
  theme(legend.position = "none")

congruenceViolins <- plot_grid(vioPlot_same_ConcGroup, 
          vioPlot_same_ID, 
          vioPlot_same_processing, 
          vioPlot_same_ratio, ncol = 4)
# 
# ggsaveAll(
#   here::here("figures/"),
#   filename = "betaCongruenceViolins_rawVSlotus",
#   plot = congruenceViolins,
#   devices = c("svg", "pdf", "png"), 
#   height = 5, 
#   width = 13
# )


mm <- lmer (data = dist_long, rank (brayDist) ~ same_ID + same_processing + same_ratio + same_ConcGroup + (1 | Sample1) + (1 | Sample2), REML = F)
varianceTable <- as.data.frame(anova (mm))
# e.g. ID explains 8%, extraction 6%, collection <<1%
varianceTable$VarExplained <- varianceTable$`Sum Sq` / sum (resid (mm)^2)
varianceTable$Variable <- rownames(varianceTable)
varianceTable[nrow(varianceTable) +1, ] <- c(rep(1, 4), (1 - sum(varianceTable$VarExplained)), "Residuals")
varianceTable$VarExplained <- as.numeric(varianceTable$VarExplained)
varianceTable$VarLabels <- scales::percent(varianceTable$VarExplained)
print("Variance table (ANOVA output:)")
print(varianceTable)
 
 #remove insignificant variable
varianceTable <- varianceTable[!rownames(varianceTable) == "same_processing", ]
  
  varianceTable2 <- varianceTable %>% 
    arrange(desc(Variable)) %>%
    mutate(prop = VarExplained / sum(varianceTable$VarExplained) *100) %>%
    mutate(ypos = cumsum(prop)- 0.5*prop )
  # varianceTable %>%
  #   mutate(VarLabels = fct_reorder(VarLabels, VarExplained)) %>%
  ggplot (varianceTable2) +
    geom_bar(aes (x = "", y = prop, fill = Variable), width = 1, stat = "identity", color = "black") +
    coord_polar("y", start = 0) +
    theme_void() +
    geom_text(aes(x=1.7, y = ypos, label=VarLabels)) ->
    pVarExpl
  pVarExpl
  
  
  # ggsaveAll(
  #   here::here("figures/"),
  #   filename = paste0("ExplainedVarPieCongruence", dataType),
  #   plot = pVarExpl,
  #   devices = c("svg", "pdf", "png"), 
  #   height = 5, 
  #   width = 5
  # )
  # 
  


#####                                                                #####
##########                                                      ##########
###############              plot wells as pie charts      ###############
##########                                                      ##########
#####                                                                #####
wellMeta <- metadata[, c("human", "bacteria", "ConcGroup")]
#wellMeta$ratio <- metadata$human / metadata$bacteria
#wellMeta$partHuman <- metadata$human / 100
#wellMeta$partBact <- metadata$bacteria / 100
wellMeta$sampleID <- rownames(wellMeta)
wellMetaLong <- tidyr::gather(data = wellMeta, key = "type", value = "percentage", human:bacteria)
wellMetaLong <- na.exclude(wellMetaLong)
wellMetaLong <- wellMetaLong[!(wellMetaLong$ConcGroup %in% c("pure_bact", "pure_human")), ]
wellMetaLong$ConcGroup <- as.numeric(wellMetaLong$ConcGroup)
wellMetaLong$ConcGroup <- round((wellMetaLong$ConcGroup * 0.01 ) / 0.407, digits = 3)

wellMetaLong <- wellMetaLong[wellMetaLong$ConcGroup == 1.0, ]

rowPies <- list()
# for (alpha in seq(1, 0.05, length.out = 6)) {
for (alpha in c(1, 0.62, 0.43, 0.24, 0.125, 0.07)) {
  pies <- list()
  for (sample in unique(wellMetaLong$sampleID)) {
    varianceTable2 <- wellMetaLong[wellMetaLong$sampleID == sample, ] %>% 
      mutate(ypos = cumsum(percentage)- 0.5*percentage )
    #print(varianceTable2)
    
    ggplot (varianceTable2) +
      geom_bar(aes (x = "", y = percentage, fill = type), alpha = alpha, width = 0.8, stat = "identity", color = "grey15", show.legend = FALSE) +
      scale_fill_manual(values = c("#964B00", "#FF6A5E")) + 
      coord_polar("y", start = 0) +
      theme(plot.margin = margin(l = -10, r = -10, unit = "cm")) +
      theme_void() ->
      pies[[sample]]
  }
  rowPies[[as.character(alpha)]] <- plot_grid(pies[["P01.A01.1"]],
                                NULL,
                                pies[["P01.B01.2"]], 
                                NULL,
                                pies[["P01.G01.7"]],
                                NULL,
                                pies[["P01.C01.3"]],
                                NULL,
                                pies[["P01.F01.6"]], 
                                NULL,
                                pies[["P01.E01.5"]], 
                                NULL,
                                pies[["P01.D01.4"]], 
                                nrow = 1, rel_widths = c(1, -0.05, 1, -0.05, 1, -0.05, 1, -0.05, 1, -0.05, 1, -0.05,  1))
}


plot_grid(rowPies[[1]],
          NULL,
          rowPies[[2]],
          NULL,
          rowPies[[3]],
          NULL,
          rowPies[[4]],
          NULL,
          rowPies[[5]],
          NULL,
          rowPies[[6]],
          ncol = 1, rel_heights = c(1, -0.05, 1, -0.05, 1, -0.05, 1, -0.05, 1, -0.05, 1)) -> allPies
allPies

# ggsaveAll(
#   here::here("figures/"),
#   filename = "wellPiesFadingNEW",
#   plot = allPies,
#   devices = c("svg", "pdf", "png"), 
#   height = 8, 
#   width = 8
# )

#####                                                                #####
##########                                                      ##########
###############         deviation from ground truth        ###############
##########                                                      ##########
#####                                                                #####

deviate_rel_raw <- molten_DF_list$gen_raw_rel
deviate_rel_raw <- deviate_rel_raw[!(deviate_rel_raw$Phylum %in% c("other", "noHit", "Cutibacterium")), ]
#deviate_rel_raw_pureBact <- deviate_rel_raw[deviate_rel_raw$ratio == "pure_bact", ]
deviate_rel_raw_pureBact <- deviate_rel_raw[deviate_rel_raw$variable == "rel", ]
#deviate_rel_raw <- deviate_rel_raw[!(deviate_rel_raw$ratio %in% c("pure_bact", "theoretical", "empty", "pure_human")), ]
deviate_rel_raw <- deviate_rel_raw[!(deviate_rel_raw$ratio %in% c("empty", "pure_human")), ]

deviate_rel_raw2 <- deviate_rel_raw
for (bact in unique(deviate_rel_raw2$Phylum)) {
  deviate_rel_raw2[deviate_rel_raw2$Phylum == bact, "value"] <- deviate_rel_raw2[deviate_rel_raw2$Phylum == bact, "value"]/deviate_rel_raw_pureBact[deviate_rel_raw_pureBact$Phylum == bact, "value"]
}

deviate_rel_raw2 <- deviate_rel_raw2[!is.na(deviate_rel_raw2$ratio), ]
ratio_compare_to_theoretical <- ggplot (deviate_rel_raw2, aes(x = ratio, y = value)) +  #subset(deviate_rel_raw2, ratio == "pure_bact")
  #geom_text(aes(label = scales::number(value, accuracy = 0.01)), nudge_x = 0.11) +
  theme_classic () +
  geom_hline(yintercept = 1, color = "grey") +
  geom_point(aes (x = ratio, y = value, fill = Phylum), pch = 21, size = 5) +
  scale_fill_manual(values = cbPalette[-c(2, 7, 8)]) +
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 15)) +
  scale_y_continuous(trans='log10') + 
  ylab("ratio sample : theoretical")

ggsaveAll(
  here::here("figures/"),
  filename = "ratio_compare_to_theoretical",
  plot = ratio_compare_to_theoretical,
  devices = c("svg", "pdf", "png"),
  height = 8,
  width = 8
)

  

# mean over concGroup
deviate_rel_raw3 <- deviate_rel_raw_pureBact[, c("Phylum", "ratio", "value")]

for (bact in unique(deviate_rel_raw$Phylum)) {
  bactTemp <- deviate_rel_raw[deviate_rel_raw$Phylum == bact, ]
  for (ratio in unique(deviate_rel_raw$ratio)) {
    ratioTemp <- bactTemp[bactTemp$ratio == ratio, ]
    ratioTempMean <- mean(ratioTemp$value)
    deviate_rel_raw3 <- rbind(deviate_rel_raw3, c(ratioTemp$Phylum[1], ratioTemp$ratio[1], ratioTempMean))
  }
  
}
deviate_rel_raw3 <- deviate_rel_raw3[deviate_rel_raw3$ratio != "pure_bact", ]
deviate_rel_raw3$value <- as.numeric(deviate_rel_raw3$value)
deviate_rel_raw4 <- deviate_rel_raw3

for (bact in unique(deviate_rel_raw2$Phylum)) {
  deviate_rel_raw4[deviate_rel_raw4$Phylum == bact, "value"] <- deviate_rel_raw3[deviate_rel_raw3$Phylum == bact, "value"]/deviate_rel_raw_pureBact[deviate_rel_raw_pureBact$Phylum == bact, "value"]
}
deviate_rel_raw4$processing <- "raw"


deviate_rel_raw_meanedPlot <- ggplot (deviate_rel_raw4, aes(x = ratio, y = value)) +
  geom_boxplot(aes (x = ratio, y = value, fill = ratio)) +
  geom_point(aes (x = ratio, y = value, color = Phylum, size = 3)) +
  geom_hline(yintercept=1) +
  theme_classic () + 
  labs(title = "raw meaned") +
  ylim(0.75, 1.5)






deviate_rel_decontLotus <- molten_DF_list$gen_decontLotus_rel
deviate_rel_decontLotus <- deviate_rel_decontLotus[!(deviate_rel_decontLotus$Phylum %in% c("other", "noHit", "Cutibacterium")), ]
deviate_rel_decontLotus_pureBact <- deviate_rel_decontLotus[deviate_rel_decontLotus$ratio == "pure_bact", ]
deviate_rel_decontLotus <- deviate_rel_decontLotus[!(deviate_rel_decontLotus$ratio %in% c("pure_bact", "theoretical", "empty", "pure_human")), ]
deviate_rel_decontLotus2 <- deviate_rel_decontLotus
for (bact in unique(deviate_rel_decontLotus2$Phylum)) {
  deviate_rel_decontLotus2[deviate_rel_decontLotus2$Phylum == bact, "value"] <- 
    deviate_rel_decontLotus2[deviate_rel_decontLotus2$Phylum == bact, "value"]/
    deviate_rel_decontLotus_pureBact[deviate_rel_decontLotus_pureBact$Phylum == bact, "value"]
}


ggplot (deviate_rel_decontLotus2, aes(x = ratio, y = value)) +
  geom_boxplot(aes (x = ratio, y = log(value), fill = Phylum)) +
  theme_classic ()

# mean over concGroup
deviate_rel_decontLotus3 <- deviate_rel_decontLotus_pureBact[, c("Phylum", "ratio", "value")]

for (bact in unique(deviate_rel_decontLotus$Phylum)) {
  bactTemp <- deviate_rel_decontLotus[deviate_rel_decontLotus$Phylum == bact, ]
  for (ratio in unique(deviate_rel_decontLotus$ratio)) {
    ratioTemp <- bactTemp[bactTemp$ratio == ratio, ]
    ratioTempMean <- mean(ratioTemp$value)
    deviate_rel_decontLotus3 <- rbind(deviate_rel_decontLotus3, c(ratioTemp$Phylum[1], ratioTemp$ratio[1], ratioTempMean))
  }
  
}
deviate_rel_decontLotus3 <- deviate_rel_decontLotus3[deviate_rel_decontLotus3$ratio != "pure_bact", ]
deviate_rel_decontLotus3$value <- as.numeric(deviate_rel_decontLotus3$value)
deviate_rel_decontLotus4 <- deviate_rel_decontLotus3

for (bact in unique(deviate_rel_decontLotus2$Phylum)) {
  deviate_rel_decontLotus4[deviate_rel_decontLotus4$Phylum == bact, "value"] <- deviate_rel_decontLotus3[deviate_rel_decontLotus3$Phylum == bact, "value"]/deviate_rel_decontLotus_pureBact[deviate_rel_decontLotus_pureBact$Phylum == bact, "value"]
}
deviate_rel_decontLotus4$processing <- "decontLotus"

deviate_rel_decontLotus_meanedPlot <- ggplot (deviate_rel_decontLotus4, aes(x = ratio, y = value)) +
  geom_boxplot(aes (x = ratio, y = value, fill = ratio)) +
  geom_point(aes (x = ratio, y = value, color = Phylum, size = 3)) +
  geom_hline(yintercept=1) +
  theme_classic () + 
  labs(title = "decontLotus meaned") +
  ylim(0.75, 1.5)

#plot_grid(deviate_rel_raw_meanedPlot, deviate_rel_decontLotus_meanedPlot, align = "hv")


deviate_rel_both <- rbind(deviate_rel_decontLotus4, deviate_rel_raw4)
deviate_rel_both$label <- paste0(deviate_rel_both$Phylum, "_", deviate_rel_both$ratio)

deviate_rel_both_meanedPlot <- ggplot (deviate_rel_both, aes(x = ratio, y = value)) +
  geom_boxplot(aes (fill = processing), alpha = 0.5) +
  #geom_dotplot(dotsize = 0.5, binaxis = "y", binpositions = "all", aes(fill = Phylum), stackgroups = T) +
  geom_point(aes (color = Phylum, group = processing), alpha = 0.9, size = 3, position = position_dodge(0.75)) +
  #geom_line(aes(group = label, color = Phylum)) +
  geom_hline(yintercept=1, color = "grey") +
  theme_classic () + 
  ylab("relative abundance difference compared to pure_bact") +
  xlab("bacterial:human DNA ratio group") +
  ylim(0.75, 1.5) +
  scale_color_manual(values = cbPalette[-c(2, 7, 8)]) +
  scale_fill_discrete(name = "off-target read removal", labels = c("LotuS", "none")) +
  labs(color = "genus")
deviate_rel_both_meanedPlot

ggsaveAll(
  here::here("figures/"),
  filename = "rel_dif_compared_to_pure_bact_per_ratio_group",
  plot = deviate_rel_both_meanedPlot,
  devices = c("svg", "pdf", "png"),
  height = 5,
  width = 6
)

for (rat in unique(deviate_rel_both$ratio)) {
  print(rat)
  raw <- deviate_rel_both$value[deviate_rel_both$ratio == rat & deviate_rel_both$processing == "raw"]
  lotus <- deviate_rel_both$value[deviate_rel_both$ratio == rat & deviate_rel_both$processing == "decontLotus"]
  #print(wilcox.test(x = abs(raw), y = abs(lotus), paired = F, conf.int = T))
  print(t.test(x = abs(raw), y = abs(lotus), paired = T, conf.int = T, var.equal = F))
  
  if (rat == "1:99") {
    stop()
  }
}
x <- raw
y <- lotus
nx <- length(x)
ny <- length(y)

# 20231204 TB
backingfileCompl <- tempfile()
backingfileCompl <- path.expand(backingfileCompl)
splitTemp <- strsplit(backingfileCompl, "[\\/\\\\]")
backingfile <- splitTemp[[1]][length(splitTemp[[1]])]
dom <- bigmemory::big.matrix(nrow = nx,
                             ncol = ny,
                             backingfile = backingfile,
                             descriptorfile = paste0(backingfile, ".desc", collapse = ""),
                             backingpath = tempdir()
)

# computation of cliff's delta
nxny <- length(dom)
preSumPSc <- vector(length = ncol(dom))
preSumBiggers <- vector(length = ncol(dom))
for (i in 1:ncol(dom)) {
  preSumPSc[i] <- sum(dom[, i] < 0)
  preSumBiggers[i] <- sum(dom[, i] > 0)
}
rm(dom)
unlink(paste0(backingfileCompl, "*", collapse = "")) # 20231204 TB
PSc <- sum(preSumPSc)/nxny
PSbigger <- sum(preSumBiggers)/nxny
dc <- PSc - PSbigger




#####                                                                #####
##########                                                      ##########
###############                    metadeconf raw vs decont      ###############
##########                                                      ##########
#####                                                                #####

features <- beta_meta[beta_meta$processing %in% c("raw"), c(1:11)]
metadata <- beta_meta[beta_meta$processing %in% c("raw"), c(12:ncol(beta_meta))]
metadata <- metadata[, c("processing", "Concentration", "human", "ratio", "ID", "ConcGroup")]
metadata$processing[metadata$processing == "raw"] <- 0
metadata$processing[metadata$processing == "decontLotus"] <- 1
metadata$processing <- as.numeric(metadata$processing)
metadata <- na.exclude(metadata)
features <- features[rownames(metadata), ]
metadata$processing[metadata$ratio == "pure_bact"] <- 1
stop("not running metadeconounfdR automatically!")
testMetadecon <- MetaDeconfound(features, metadata, nnodes = 7, returnLong = T, robustCutoff = 1)
testMetadecon$Ds[testMetadecon$metaVariable == "human"] <- testMetadecon$Ds[testMetadecon$metaVariable == "human"]*(-1)

prettyMeta <- as.data.frame(colnames(metadata))
prettyMeta$niceNames <- prettyMeta$`colnames(metadata)`
prettyMeta$niceNames[prettyMeta$niceNames == "Concentration"] <- "total DNA concentration"
prettyMeta$niceNames[prettyMeta$niceNames == "human"] <- "percent human DNA"
metadeconfPlot <- BuildHeatmap(testMetadecon, d_range = "full", metaVariableNames = prettyMeta) + coord_flip()
metadeconfPlot <-
  metadeconfPlot + theme(
    axis.title = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12, face = "italic", angle = 45, vjust = 1)
  )
metadeconfPlot

ggsaveAll(
  here::here("figures/"),
  filename = "metadeconf_concVSratio",
  plot = metadeconfPlot,
  devices = c("svg", "pdf", "png"),
  height = 2.9,
  width = 7
)

# only run metadeconfoundR on decontLotus set, so we see pure associations after filtering
  # (the associations to  "human" in the above combined raw+lotus set are based on a mixture of both sets)
features <- beta_meta[beta_meta$processing %in% c("decontLotus"), c(1:11)]
metadata <- beta_meta[beta_meta$processing %in% c("decontLotus"), c(12:ncol(beta_meta))]
metadata <- metadata[, c("pseudo", "Concentration", "human", "ratio", "ConcGroup")]
metadata <- na.exclude(metadata)
features <- features[rownames(metadata), ]

testMetadeconLotus <- MetaDeconfound(features, metadata, nnodes = 7, returnLong = T)
BuildHeatmap(testMetadeconLotus, d_range = "full", keepMeta = colnames(metadata))

mergedFeMe <- merge(metadata, features, by = 0)

ggplot (mergedFeMe) +
  geom_violin(aes (x = as.factor(pseudo), y = rank(other)), draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_point(aes (x = as.factor(pseudo), y = rank(other), color = ratio)) +
  scale_color_manual(values = cbPalette)

View(mergedFeMe[, c("Row.names", "pseudo", "other", "Cutibacterium")])

table(mergedFeMe$other, mergedFeMe$pseudo)

table(mergedFeMe$Concentration, mergedFeMe$pseudo)

ggplot (mergedFeMe) +
  geom_violin(aes (x = as.factor(pseudo), y = rank(Concentration)), draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_point(aes (x = as.factor(pseudo), y = rank(Concentration)))


#####                                                                #####
##########                                                      ##########
###############                    digital qpcr           ###############
##########                                                      ##########
#####                                                                #####

#qpcr <- read.table(here::here("ddqpcr/all_ddqpcr_results.csv"), sep = ",", header = T)
qpcr <- read.table(here::here("ddqpcr/19124_all_ddqpcr_results.csv"), sep = ",", header = T)

# normalize 16S and 18S values to the 50:50 samples --> 16S and 18S should have dsame level for the 50:50 sample
qpcr$all_18S <- qpcr$all_18S*(mean(c(25520000, 21700000))/mean(c(106020, 77680)))
qpcr_long <- melt(data = qpcr)
qpcr_long$variable <- as.character(qpcr_long$variable)
qpcr_long$variable[qpcr_long$variable == "all_16S"] <- "16S"
qpcr_long$variable[qpcr_long$variable == "all_18S"] <- "18S"
qpcr_long$verdünnung[qpcr_long$verdünnung == "1 zu 99"] <- "99:1"
qpcr_long$verdünnung[qpcr_long$verdünnung == "1 zu 90"] <- "90:10"
qpcr_long$verdünnung[qpcr_long$verdünnung == "25 zu 75"] <- "75:25"
qpcr_long$verdünnung[qpcr_long$verdünnung == "50 zu 50 "] <- "50:50"
qpcr_long$verdünnung[qpcr_long$verdünnung == "75 zu 25"] <- "25:75"
qpcr_long$verdünnung[qpcr_long$verdünnung == "90 zu 1"] <- "10:90"
qpcr_long$verdünnung[qpcr_long$verdünnung == "99 zu 1"] <- "1:99"
qpcr_long$verdünnung <-
  factor(
    qpcr_long$verdünnung,
    levels = c(
      "99:1",
      "90:10",
      "75:25",
      "50:50",
      "25:75",
      "10:90",
      "1:99"
    ))

digi_qpcr_centered_to_50_50 <-
  ggplot(data = qpcr_long, 
         aes(x = verdünnung,
             y = log10(value), 
             fill = variable)) +
  geom_point(aes(size = 4), pch = 21) +
  #geom_line() +
  theme_classic() +
  xlab("ratio") +
  ylab(expression("log10 transformed normalized copy number")) +
  guides(size = "none") +
  theme(legend.title = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.position = c(0.3,0.2))

digi_qpcr_centered_to_50_50
mean(qpcr_long[7:8, "value"])/mean(qpcr_long[21:22, "value"])


ggsaveAll(
  here::here("figures/"),
  filename = paste0("digi_qpcr_centered_to_50_50"),
  plot = digi_qpcr_centered_to_50_50,
  devices = c("svg", "pdf", "png"),
  height = 3.9,
  width = 3.9
)
