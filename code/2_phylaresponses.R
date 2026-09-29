# 2: Dominant phyla GAMs ====

source("./code/initial_setup.R")

## Prokaryotic phyla ====
## Grouping by phyla (warning is due to some NAs in taxonomy)
bacphy_df = rowsum(t(bac_rare), clean_bac_taxonomy$Phylum)

## Getting relative abundances
bac_raresize = apply(bac_rare, 1, sum)
bacphy_rel = apply(bacphy_df, 1, function(x) {
  x/bac_raresize})

## Only keep taxa that make up at least 10% of community in at least 1 sample
bacphy_keep = apply(bacphy_rel, 2, function(x) sum(x>0.1))
bacphy_keep = bacphy_keep[bacphy_keep>0]
bacphy_hi = bacphy_rel[,colnames(bacphy_rel) %in% names(bacphy_keep) ,drop = FALSE]

bacphy_hi = as.data.frame(bacphy_hi)
bacphy_hi = cbind(bacphy_hi, 
                  swd = main_df$swd,
                  timepoint = main_df$timepoint,
                  treatment = main_df$treatment)
bacphy_hi = melt(bacphy_hi[1:110,], id.vars = c("swd", "timepoint", "treatment"))
bacphy_hi = droplevels(bacphy_hi)
bacphy_hi$int_term = interaction(bacphy_hi$timepoint, 
                                 bacphy_hi$treatment,
                                 bacphy_hi$variable)

# Running the GAM (slow)
bacphy_mod = gam(value ~ swd + int_term +
                   ti(swd, by = int_term, k = 6, bs = "tp"),
                 method = "REML",
                 select = TRUE,
                 data = bacphy_hi)

## Creating dataframe to predict values for
bacphy_pred = data.frame(swd = rep(seq(min(bacphy_hi$swd), 
                                       max(bacphy_hi$swd), 
                                       length.out = 100), 2*5),
                         int_term = rep(levels(bacphy_hi$int_term), each = 100),
                         stringsAsFactors = TRUE)
bacphy_pred2 = get_gam_ci(bacphy_mod, bacphy_pred)

# Breaking apart the int_term
temp = unlist(strsplit(as.character(bacphy_pred$int_term), "\\."))
temp = matrix(temp, ncol = 3, byrow = TRUE)
temp = as.data.frame(temp); names(temp) = c("timepoint", "treatment", "phylum")
bacphy_pred2 = cbind(temp, swd = bacphy_pred$swd, bacphy_pred2)
head(bacphy_pred2)

# timepoints
bacphy_pred2$timepoint = factor(bacphy_pred2$timepoint, 
                                levels = c("3 days", "7 days", "14 days",
                                           "35 days", "70 days"))
# treatments
bacphy_pred2$treatment = factor(bacphy_pred2$treatment,
                                levels = c("field",
                                           "drought"))
bacphy_pred2$phylum = factor(bacphy_pred2$phylum,
                             levels = levels(bacphy_hi$variable))

## Fig 2e:
facet_names = c(
  `field` = "dry-down",
  `drought` = "rewet-up",
  `3 days` = "3 days",
  `7 days` = "7 days",
  `14 days` = "14 days",
  `35 days` = "35 days",
  `70 days` = "70 days"
)

bac_cols = c("#648FFF",
             "#785EF0",
             "#DC267F",
             "#FE6100",
             "#FFB000",
             "#4CE8A2")

bacphy_legend2 = get_legend(
  ggplot(data = bacphy_pred2)+
    geom_line(aes(x = swd, y = fit, col = phylum),
              linewidth = 2)+
    scale_color_manual(values = bac_cols,
                       breaks = c("Acidobacteriota",
                                  "Actinomycetota",
                                  "Bacillota",
                                  "Chloroflexota",
                                  "Pseudomonadota",
                                  "Verrucomicrobiota"))+
    theme(legend.text = element_text(size = 14),
          legend.title = element_text(size = 16)))

fig2e = ggplot(data = bacphy_pred2[bacphy_pred2$timepoint %in% c("70 days"),])+
  geom_ribbon(aes(x = swd, 
                  ymin = sim_lo, 
                  ymax = sim_hi,
                  fill = phylum), alpha = 0.3)+
  geom_line(data = bacphy_pred2[bacphy_pred2$timepoint %in% c("70 days") &
                                  bacphy_pred2$treatment == "field",],
            aes(x = swd, y = fit, col = phylum),
            linewidth = 1.5, 
            arrow = arrow(ends = "last", type = "open", 
                          length = unit(0.3, "cm")))+
  geom_line(data = bacphy_pred2[bacphy_pred2$timepoint %in% c("70 days") &
                                  bacphy_pred2$treatment == "drought",],
            aes(x = swd, y = fit, col = phylum),
            linewidth = 1.5, 
            arrow = arrow(ends = "first", type = "open", 
                          length = unit(0.3, "cm")))+
  geom_point(data = bacphy_hi[bacphy_hi$timepoint %in% c("70 days"),],
             aes(x = swd, y = value, col = variable))+
  facet_wrap(.~treatment, nrow = 2,
             labeller = as_labeller(facet_names))+
  scale_fill_manual(values = bac_cols,
                    breaks = c("Acidobacteriota",
                               "Actinomycetota",
                               "Bacillota",
                               "Chloroflexota",
                               "Pseudomonadota",
                               "Verrucomicrobiota"))+
  scale_color_manual(values = bac_cols,
                     breaks = c("Acidobacteriota",
                                "Actinomycetota",
                                "Bacillota",
                                "Chloroflexota",
                                "Pseudomonadota",
                                "Verrucomicrobiota"))+
  scale_y_continuous(labels = function(x) sprintf("%.2f",x))+ # 2dp to be align with fungal plot
  labs(x = "SWD", y = "Proportion", title = "Prokaryotes")+
  guides(col = "none", fill = "none")+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 22,
                                  vjust = 2))

fig2e_final = grid.arrange(fig2e, bacphy_legend2, nrow = 1,
                     widths = c(5,3))
plot(fig2e_final)

## Fungal phyla ====
funphy_df = rowsum(t(fun_rare), clean_fun_taxonomy$Phylum)
fun_raresize = apply(fun_rare, 1, sum)
funphy_rel = apply(funphy_df, 1, function(x) {
  x/fun_raresize})

## Only keep taxa that make up at least 10% of community in at least 1 sample
funphy_keep = apply(funphy_rel, 2, function(x) sum(x>0.1))
funphy_keep = funphy_keep[funphy_keep>0]
funphy_hi = funphy_rel[,colnames(funphy_rel) %in% names(funphy_keep) ,drop = FALSE]

funphy_hi = as.data.frame(funphy_hi)
funphy_hi = cbind(funphy_hi, 
                  swd = main_df$swd,
                  timepoint = main_df$timepoint,
                  treatment = main_df$treatment)
funphy_hi = melt(funphy_hi[1:110,],
                 id.vars = c("swd", "timepoint", "treatment"))
funphy_hi = droplevels(funphy_hi)
funphy_hi$int_term = interaction(funphy_hi$timepoint, 
                                 funphy_hi$treatment,
                                 funphy_hi$variable)

# Running the GAM
funphy_mod = gam(value ~ swd + int_term +
                   ti(swd, by = int_term, k = 5, bs = "tp"),
                 method = "REML",
                 select = TRUE,
                 family = tw,
                 data = funphy_hi)

## Creating dataframe to predict values for
funphy_pred = data.frame(swd = rep(seq(min(funphy_hi$swd), 
                                       max(funphy_hi$swd), 
                                       length.out = 100), 2*5),
                         int_term = rep(levels(funphy_hi$int_term), each = 100),
                         stringsAsFactors = TRUE)
funphy_pred2 = get_gam_ci(funphy_mod, funphy_pred)

# Breaking apart the int_term
temp = unlist(strsplit(as.character(funphy_pred$int_term), "\\."))
temp = matrix(temp, ncol = 3, byrow = TRUE)
temp = as.data.frame(temp); names(temp) = c("timepoint", "treatment", "phylum")
funphy_pred2 = cbind(temp, swd = funphy_pred$swd, funphy_pred2)
head(funphy_pred2)

# timepoints
funphy_pred2$timepoint = factor(funphy_pred2$timepoint, 
                                levels = c("3 days", "7 days", "14 days",
                                           "35 days", "70 days"))
# treatments
funphy_pred2$treatment = factor(funphy_pred2$treatment,
                                levels = c("field",
                                           "drought"))
funphy_pred2$phylum = factor(funphy_pred2$phylum,
                             levels = levels(funphy_hi$variable))

## Backtransforming from Tweedie distribution
funphy_pred2$fit = funphy_mod$family$linkinv(funphy_pred2$fit)
funphy_pred2$point_hi = funphy_mod$family$linkinv(funphy_pred2$point_hi)
funphy_pred2$point_lo = funphy_mod$family$linkinv(funphy_pred2$point_lo)
funphy_pred2$sim_hi = funphy_mod$family$linkinv(funphy_pred2$sim_hi)
funphy_pred2$sim_lo = funphy_mod$family$linkinv(funphy_pred2$sim_lo)

## Fig 2f:
fun_cols = c("#648FFF",
             "#DC267F",
             "#FFB000")

funphy_legend2 = get_legend(
  ggplot(data = funphy_pred2)+
    geom_line(aes(x = swd, y = fit, colour = phylum),
              linewidth = 2) + 
    scale_colour_manual(values = fun_cols,
                        labels = sub("p__", "", levels(funphy_pred2$phylum)))+
    theme(legend.text = element_text(size = 14),
          legend.title = element_text(size = 16)))

fig2f = ggplot(data = funphy_pred2[funphy_pred2$timepoint %in% c("70 days"),])+
  geom_ribbon(aes(x = swd, 
                  ymin = sim_lo, 
                  ymax = sim_hi,
                  fill = phylum), alpha = 0.3)+
  geom_line(data = funphy_pred2[funphy_pred2$timepoint %in% c("70 days") & 
                                  funphy_pred2$treatment == "field",],
            aes(x = swd, y = fit, col = phylum),
            linewidth = 1.5, 
            arrow = arrow(ends = "last", type = "open", 
                          length = unit(0.3, "cm")))+
  geom_line(data = funphy_pred2[funphy_pred2$timepoint %in% c("70 days") &
                                  funphy_pred2$treatment == "drought",],
            aes(x = swd, y = fit, col = phylum),
            linewidth = 1.5, 
            arrow = arrow(ends = "first", type = "open", 
                          length = unit(0.3, "cm")))+
  geom_point(data = funphy_hi[funphy_hi$timepoint %in% c("70 days"),],
             aes(x = swd, y = value, col = variable))+
  facet_wrap(.~treatment, nrow = 2,
             labeller = as_labeller(facet_names))+
  labs(x = "SWD", y = "Proportion", title = "Fungi")+
  scale_fill_manual(values = fun_cols,
                    labels = sub("p__", "", levels(funphy_pred2$phylum)))+
  scale_colour_manual(values = fun_cols,
                      labels = sub("p__", "", levels(funphy_pred2$phylum)))+
  guides(col = "none", fill = "none")+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 22,
                                  vjust = 2))

fig2f_final = grid.arrange(fig2f, funphy_legend2, nrow = 1,
                     widths = c(5,3))
plot(fig2f)

# ggsave("./figures/fig2e.svg", fig2e_final, height = 6, width = 6)
# ggsave("./figures/fig2f.svg", fig2f_final, height = 6, width = 6)

# Additional SI figures (Figs S4-5) ====

###(a) Fig S4: Prokaryotic phyla all timepoints ====
bacphy_legendSI = get_legend(
  ggplot(data = bacphy_pred2)+
    geom_line(aes(x = swd, y = fit, col = phylum),
              linewidth = 2)+
    scale_color_manual(values = bac_cols,
                       breaks = c("Acidobacteriota",
                                  "Actinomycetota",
                                  "Bacillota",
                                  "Chloroflexota",
                                  "Pseudomonadota",
                                  "Verrucomicrobiota"))+
    theme(legend.text = element_text(size = 12),
          legend.title = element_text(size = 14)))

figS4 = ggplot(data = bacphy_pred2)+
  geom_ribbon(aes(x = swd, 
                  ymin = sim_lo, 
                  ymax = sim_hi,
                  fill = phylum), alpha = 0.3)+
  geom_line(data = bacphy_pred2[bacphy_pred2$treatment == "field",],
            aes(x = swd, y = fit, col = phylum),
            linewidth = 1.5, 
            arrow = arrow(ends = "last", type = "open", 
                          length = unit(0.3, "cm")))+
  geom_line(data = bacphy_pred2[bacphy_pred2$treatment == "drought",],
            aes(x = swd, y = fit, col = phylum),
            linewidth = 1.5, 
            arrow = arrow(ends = "first", type = "open", 
                          length = unit(0.3, "cm")))+
  geom_point(data = bacphy_hi,
             aes(x = swd, y = value, col = variable))+
  facet_wrap(.~treatment*timepoint, nrow = 2,
             labeller = as_labeller(facet_names))+
  scale_fill_manual(values = bac_cols,
                    breaks = c("Acidobacteriota",
                               "Actinomycetota",
                               "Bacillota",
                               "Chloroflexota",
                               "Pseudomonadota",
                               "Verrucomicrobiota"))+
  scale_color_manual(values = bac_cols,
                     breaks = c("Acidobacteriota",
                                "Actinomycetota",
                                "Bacillota",
                                "Chloroflexota",
                                "Pseudomonadota",
                                "Verrucomicrobiota"))+
  labs(x = "SWD", y = "Proportion")+
  guides(col = "none", fill = "none")+
  theme(axis.text.x = element_text(size = 8),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 14))

figS4_final = grid.arrange(figS4, bacphy_legendSI, 
                           nrow = 1,
                     widths = c(5,1))
plot(figS4_final)

#ggsave("./figures/figS4.svg", figS4_final, width=11, height=6)

### (b) Fig S5: Fungal phyla all timepoints ====
funphy_legendSI = get_legend(
  ggplot(data = funphy_pred2)+
    geom_line(aes(x = swd, y = fit, colour = phylum),
              linewidth = 2) + 
    scale_colour_manual(values = fun_cols,
                        labels = sub("p__", "", levels(funphy_pred2$phylum)))+
    theme(legend.text = element_text(size = 12),
          legend.title = element_text(size = 14)))

figS5 = ggplot(data = funphy_pred2)+
  geom_ribbon(aes(x = swd, 
                  ymin = sim_lo, 
                  ymax = sim_hi,
                  fill = phylum), alpha = 0.3)+
  geom_line(data = funphy_pred2[funphy_pred2$treatment == "field",],
            aes(x = swd, y = fit, col = phylum),
            linewidth = 1.5, 
            arrow = arrow(ends = "last", type = "open", 
                          length = unit(0.3, "cm")))+
  geom_line(data = funphy_pred2[funphy_pred2$treatment == "drought",],
            aes(x = swd, y = fit, col = phylum),
            linewidth = 1.5, 
            arrow = arrow(ends = "first", type = "open", 
                          length = unit(0.3, "cm")))+
  geom_point(data = funphy_hi,
             aes(x = swd, y = value, col = variable))+
  facet_wrap(.~treatment*timepoint, nrow = 2,
             labeller = as_labeller(facet_names))+
  labs(x = "SWD", y = "Proportion")+
  scale_fill_manual(values = fun_cols,
                    labels = sub("p__", "", levels(funphy_pred2$phylum)))+
  scale_colour_manual(values = fun_cols,
                      labels = sub("p__", "", levels(funphy_pred2$phylum)))+
  guides(col = "none", fill = "none")+
  theme(axis.text.x = element_text(size = 8),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 14))

figS5_final = grid.arrange(figS5, funphy_legendSI, nrow = 1,
                           widths = c(5,1))
plot(figS5)

#ggsave("./figures/figS5.svg", figS5_final, width=11, height=6)

