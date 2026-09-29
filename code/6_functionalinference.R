## Functional inference from amplicon data
# We'll use picrust2 for 16s and funguild for ITS
library(ggpicrust2)
library(FUNGuildR)

## PICRUSt2  ====
# Analysing outputs after running standalone picrust2
metacyc = pathway_annotation(file = "./data/picrust2_outputs/pathways_out/path_abun_unstrat.tsv.gz",
                          pathway = "MetaCyc")
metacyc = as.data.frame(metacyc)

# Checking NSTI
sample_nsti = read.delim("./data/picrust2_outputs/EC_metagenome_out/weighted_nsti.tsv.gz", 
                         header = TRUE, sep = "\t")
asv_nsti = read.delim("./data/picrust2_outputs/combined_marker_predicted_and_nsti.tsv.gz", 
                   header = TRUE, sep = "\t")
summary(sample_nsti)
plot(sample_nsti$weighted_NSTI)

asv_reads = apply(clean_bac,2,sum)
nsti_df = data.frame(asv = names(asv_reads), 
                     counts = unname(asv_reads),
                     rel_abun = unname(asv_reads/sum(asv_reads)))
nsti_df = left_join(nsti_df, asv_nsti,
                    by = join_by(asv == sequence))
sum(nsti_df$rel_abun[nsti_df$metadata_NSTI < 2]) *100 # Proportion of reads with NSTI <2
sum(nsti_df$rel_abun[nsti_df$metadata_NSTI < 0.2]) *100 # with NSTI <0.2


## NMDS
metacyc_spp = t(as.matrix(metacyc[-c(1:2)]))

set.seed(1630)
func_nmds = metaMDS(metacyc_spp, try = 100, trymax = 1000)
func_scores = scores(func_nmds)$sites

nmds_df = data.frame(nmds1 = func_scores[,1],
                     nmds2 = func_scores[,2],
                     timepoint = main_df$timepoint,
                     treatment = main_df$treatment,
                     swd = main_df$swd)

figS13 = ggplot(nmds_df[1:110,])+
  # Dry-down points
  geom_point(data = nmds_df[nmds_df$treatment == "field",], 
             aes(x = nmds1, y = nmds2, col = swd), size = 4)+
  scale_colour_gradient(name = "dry-down SWD",
                        low = "#025d96", high = "#c8e7fa")+
  geom_path(data = nmds_df[nmds_df$treatment == "field",],
            aes(x = nmds1, y = nmds2, col = swd),
            linewidth = 2)+
  new_scale_colour()+
  # Rewet-up points
  geom_point(data = nmds_df[nmds_df$treatment == "drought",],
             aes(x = nmds1, y = nmds2, col = swd), size = 4)+
  geom_path(data = nmds_df[nmds_df$treatment == "drought",],
            aes(x = nmds1, y = nmds2, col = swd),
            linewidth = 2)+
  scale_colour_gradient(name = "rewet-up SWD", 
                        low = "#ad0303", high = "#f5d7d7")+
  facet_wrap(.~timepoint)+
  labs(x = "NMDS1", y = "NMDS2", title = "Inferred MetaCyc pathway abundances (Prokaryotes)")+
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 12))+
  coord_fixed()

#ggsave("./figures/figS13.svg", figS14, width=8, height=5)

## Funguild ====
my_taxonomy = paste(fun_taxonomy$Kingdom, fun_taxonomy$Phylum,
                    fun_taxonomy$Class, fun_taxonomy$Order,
                    fun_taxonomy$Family, fun_taxonomy$Genus,
                    fun_taxonomy$Species,
                    sep = ";")
my_taxonomy = data.frame(asv = rownames(fun_taxonomy),
                         Taxonomy = my_taxonomy)
my_guilds = funguild_assign(my_taxonomy)

## Checking assignment coverage
funconf_df = apply(fun_rare, 2, sum)
funconf_df = data.frame(asv = names(funconf_df),
                       counts = unname(funconf_df),
                       rel_abun = unname(funconf_df)/sum(funconf_df))
funconf_df = left_join(funconf_df, my_guilds,
                       by = join_by(asv))

# Probable and higher
sum(funconf_df$rel_abun[funconf_df$confidenceRanking %in% c("Probable", "Highly Probable")]) * 100

## Getting counts
fun_counts = t(fun_rare)
fun_counts = apply(fun_counts, 2, function(x) x/sum(x))
fun_counts = as.data.frame(fun_counts)
fun_counts$asv = rownames(fun_counts)
fun_counts = left_join(fun_counts, my_guilds, by = join_by(asv))

# Keep only probable assignments
fun_counts = fun_counts[fun_counts$confidenceRanking %in% c("Probable", "Highly Probable"),]
fun_counts

## Assigning to major trophic modes
fun_counts$pathotroph = grepl("Pathotroph", fun_counts$trophicMode)
fun_counts$symbiotroph = grepl("Symbiotroph", fun_counts$trophicMode)
fun_counts$saprotroph = grepl("Saprotroph", fun_counts$trophicMode)

## Pathotrophs
pathotroph = fun_counts[fun_counts$pathotroph == TRUE, 1:119]
pathotroph = apply(pathotroph, 2, sum)

## Symbiotrophs
symbiotroph = fun_counts[fun_counts$symbiotroph == TRUE, 1:119]
symbiotroph = apply(symbiotroph, 2, sum)

## Symbiotrophs
saprotroph = fun_counts[fun_counts$saprotroph == TRUE, 1:119]
saprotroph = apply(saprotroph, 2, sum)

guild_df = data.frame(pathotroph,
                      symbiotroph,
                      saprotroph,
                      timepoint = main_df$timepoint,
                      treatment = main_df$treatment,
                      swd = main_df$swd)

## Fitting GAMs for figure
pathotroph_gam = gam(pathotroph ~ s(swd, k = 5) +
                       s(swd, treatment, bs = "sz") +
                       s(swd, timepoint, bs = "sz") +
                       s(swd, timepoint, treatment, bs = "sz"),
                     method = "REML",
                     select = TRUE,
                     data = guild_df[1:110,])
summary(pathotroph_gam)
gam.check(pathotroph_gam)

symbiotroph_gam = gam(symbiotroph ~ s(swd, k = 5) +
                       s(swd, treatment, bs = "sz") +
                       s(swd, timepoint, bs = "sz") +
                       s(swd, timepoint, treatment, bs = "sz"),
                     method = "REML",
                     select = TRUE,
                     data = guild_df[1:110,])
summary(symbiotroph_gam)
gam.check(symbiotroph_gam)

saprotroph_gam = gam(saprotroph ~ s(swd, k = 5) +
                        s(swd, treatment, bs = "sz") +
                        s(swd, timepoint, bs = "sz") +
                        s(swd, timepoint, treatment, bs = "sz"),
                      method = "REML",
                      select = TRUE,
                      data = guild_df[1:110,])
summary(saprotroph_gam)
gam.check(saprotroph_gam)

## Getting predictions
guild_df = droplevels(guild_df[1:110,])
guild_preds = data.frame(treatment = rep(rep(levels(guild_df$treatment), each = 1000),5),
                         swd = rep(seq(min(guild_df$swd), max(guild_df$swd), length.out = 1000), 2*5),
                         timepoint = rep(levels(guild_df$timepoint), each = 2*1000),
                         stringsAsFactors = TRUE)
guild_preds$timepoint = factor(guild_preds$timepoint, 
                               levels = c("3 days",
                                          "7 days",
                                          "14 days",
                                          "35 days",
                                          "70 days"))
guild_preds$pathotroph = get_gam_ci(pathotroph_gam, guild_preds)
guild_preds$symbiotroph = get_gam_ci(symbiotroph_gam, guild_preds)
guild_preds$saprotroph = get_gam_ci(saprotroph_gam, guild_preds)

path_plot = ggplot(data = guild_preds)+
  geom_ribbon(aes(x = swd, 
                  ymin = pathotroph$sim_lo, 
                  ymax = pathotroph$sim_hi,
                  fill = treatment), alpha = 0.3)+
  geom_point(data = guild_df[1:110,],
             aes(x = swd, y = pathotroph, colour = treatment))+
  geom_line(aes(x = swd, y = pathotroph$fit,
                colour = treatment), linewidth = 1.25)+
  facet_wrap(.~timepoint, nrow = 1)+
  scale_colour_manual(values = mycols)+
  scale_fill_manual(values = mycols)+
  lims(y = c(0,1))+
  guides(colour = "none", fill = "none")

sym_plot = ggplot(data = guild_preds)+
  geom_ribbon(aes(x = swd, 
                  ymin = symbiotroph$sim_lo, 
                  ymax = symbiotroph$sim_hi,
                  fill = treatment), alpha = 0.3)+
  geom_point(data = guild_df[1:110,],
             aes(x = swd, y = symbiotroph, colour = treatment))+
  geom_line(aes(x = swd, y = symbiotroph$fit,
                colour = treatment), linewidth = 1.25)+
  facet_wrap(.~timepoint, nrow = 1)+
  scale_colour_manual(values = mycols)+
  scale_fill_manual(values = mycols)+
  lims(y = c(0,1))+
  guides(colour = "none", fill = "none")

sap_plot = ggplot(data = guild_preds)+
  geom_ribbon(aes(x = swd, 
                  ymin = saprotroph$sim_lo, 
                  ymax = saprotroph$sim_hi,
                  fill = treatment), alpha = 0.3)+
  geom_point(data = guild_df[1:110,],
             aes(x = swd, y = saprotroph, colour = treatment))+
  geom_line(aes(x = swd, y = saprotroph$fit,
                colour = treatment), linewidth = 1.25)+
  facet_wrap(.~timepoint, nrow = 1)+
  scale_colour_manual(values = mycols)+
  scale_fill_manual(values = mycols)+
  lims(y = c(0,1))+
  guides(colour = "none", fill = "none")

figS14 = grid.arrange(path_plot,
                      sym_plot,
                      sap_plot, nrow = 3)
#ggsave("./figures/figS14.svg", figS14, width=8, height=6)
