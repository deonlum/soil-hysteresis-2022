# 1. NMDS & PERMANOVA ====

source("./code/initial_setup.R")

## Repeated rarefaction for diversity indices
## To replicate same data in paper exactly (slow)
# set.seed(4123)
# bac_avgdiv = repeat_rarefy_div(clean_bac, 1000)
# 
# set.seed(2490)
# fun_avgdiv = repeat_rarefy_div(clean_fun, 1000)

# Quick run (similar results)
bac_avgdiv = repeat_rarefy_div(clean_bac, 10)
fun_avgdiv = repeat_rarefy_div(clean_fun, 10)

## Saving out values
bac_bcavg = as.dist(bac_avgdiv$avg_matrix)
fun_bcavg = as.dist(fun_avgdiv$avg_matrix)

## Prokaryotes ====
set.seed(956)
bac_nmds = metaMDS(bac_bcavg,
                   k = 2, try = 100, 
                   trymax = 10000, maxit = 999)
bac_scores = as.data.frame(scores(bac_nmds))
main_df$bac_nmds1 = bac_scores$NMDS1
main_df$bac_nmds2 = bac_scores$NMDS2

## Generating isoclines for NMDS plots (similar to ordisurf)
# Using only the final timepoint
bac_nmds_gam = gam(swd ~ s(bac_nmds1, bac_nmds2,
                           k = 5, bs = "tp"),
                   method = "REML",
                   data = main_df[89:110,])
bac_nmds_predict = data.frame(bac_nmds1 = seq(min(main_df$bac_nmds1)*1.1, 
                                              max(main_df$bac_nmds1)*1.1, 
                                              length.out = 100),
                              bac_nmds2 = seq(min(main_df$bac_nmds2)*1.1, 
                                              max(main_df$bac_nmds2)*1.1,
                                              length.out = 100))
bac_nmds_predict = expand.grid(bac_nmds_predict)
bac_predicted_swd = predict.gam(bac_nmds_gam, newdata = bac_nmds_predict)
bac_nmds_gam_df = cbind(bac_nmds_predict, bac_predicted_swd)

## Plotting
fig1b = ggplot() +
  geom_contour(data = bac_nmds_gam_df, 
               aes(x = bac_nmds1, y = bac_nmds2, 
                   z = bac_predicted_swd),
               breaks = seq(0.3,0.9,0.05),
               col = "grey", alpha = 0.5)+
  # Control points
  geom_point(data = main_df[main_df$treatment == "control",],
             aes(x = bac_nmds1, y = bac_nmds2), 
             shape = "square", col = "black",
             size = 4, alpha = 0.5)+
  geom_point(data = main_df[main_df$treatment == "control-drought",],
             aes(x = bac_nmds1, y = bac_nmds2), 
             shape = 25, col = "black", fill = "black",
             size = 4, alpha = 0.5)+
  geom_point(data = main_df[main_df$treatment == "control-field",],
             aes(x = bac_nmds1, y = bac_nmds2), 
             shape = 17, col = "black",
             size = 4, alpha = 0.5)+
  # Dry-down points
  geom_point(data = main_df[main_df$treatment == "field",], 
             aes(x = bac_nmds1, y = bac_nmds2, col = swd), size = 4)+
  scale_colour_gradient(name = "dry-down SWD",
                        low = "#025d96", high = "#c8e7fa")+
  new_scale_colour()+
  # Rewet-up points
  geom_point(data = main_df[main_df$treatment == "drought",],
             aes(x = bac_nmds1, y = bac_nmds2, col = swd), size = 4)+
  scale_colour_gradient(name = "rewet-up SWD", 
                        low = "#ad0303", high = "#f5d7d7")+
  geom_text_contour(data = bac_nmds_gam_df, 
                    aes(x = bac_nmds1, y = bac_nmds2, 
                        z = bac_predicted_swd),
                    breaks = seq(0.3,0.9,0.05),
                    skip = 1, 
                    label.placer = label_placer_flattest())+
  labs(x = "NMDS1", y = "NMDS2", title = "Prokaryotes")+
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 22))+
  coord_fixed()

fig1b

## Fungi ====
set.seed(8213)
fun_nmds = metaMDS(fun_bcavg,
                   k = 2, try = 100, 
                   trymax = 10000, maxit = 999)
fun_scores = as.data.frame(scores(fun_nmds))
main_df$fun_nmds1 = fun_scores$NMDS1
main_df$fun_nmds2 = fun_scores$NMDS2

## Isoclines
fun_nmds_gam = gam(swd ~ s(fun_nmds1, fun_nmds2,
                           k = 5, bs = "tp"),
                   method = "REML",
                   data = main_df[89:110,])
fun_nmds_predict = data.frame(fun_nmds1 = seq(min(main_df$fun_nmds1)*1.1, 
                                              max(main_df$fun_nmds1)*1.1, 
                                              length.out = 100),
                              fun_nmds2 = seq(min(main_df$fun_nmds2)*1.1, 
                                              max(main_df$fun_nmds2)*1.1,
                                              length.out = 100))
fun_nmds_predict = expand.grid(fun_nmds_predict)
fun_predicted_swd = predict.gam(fun_nmds_gam, newdata = fun_nmds_predict)
fun_nmds_gam_df = cbind(fun_nmds_predict, fun_predicted_swd)

## Plotting
fig1c = ggplot() + 
  geom_contour(data = fun_nmds_gam_df, 
               aes(x = fun_nmds1, y = fun_nmds2, 
                   z = fun_predicted_swd),
               breaks = seq(0.3,0.9,0.05),
               col = "grey", alpha = 0.5)+
  # Dry-down points
  geom_point(data = main_df[main_df$treatment == "field",], 
             aes(x = fun_nmds1, y = fun_nmds2, col = swd), size = 4)+
  scale_colour_gradient(name = "dry-down SWD",
                        low = "#025d96", high = "#c8e7fa")+
  new_scale_colour()+
  # Rewet-up points
  geom_point(data = main_df[main_df$treatment == "drought",],
             aes(x = fun_nmds1, y = fun_nmds2, col = swd), size = 4)+
  scale_colour_gradient(name = "rewet-up SWD", 
                        low = "#ad0303", high = "#f5d7d7")+
  geom_text_contour(data = fun_nmds_gam_df, 
                    aes(x = fun_nmds1, y = fun_nmds2, 
                        z = fun_predicted_swd),
                    breaks = seq(0.3,0.9,0.05),
                    skip = 1, 
                    label.placer = label_placer_fraction(frac = 1))+
  # Control points
  geom_point(data = main_df[main_df$treatment == "control",],
             aes(x = fun_nmds1, y = fun_nmds2), 
             shape = "square", col = "black",
             size = 4, alpha = 0.5)+
  geom_point(data = main_df[main_df$treatment == "control-drought",],
             aes(x = fun_nmds1, y = fun_nmds2), 
             shape = 25, col = "black", fill = "black",
             size = 4, alpha = 0.5)+
  geom_point(data = main_df[main_df$treatment == "control-field",],
             aes(x = fun_nmds1, y = fun_nmds2), 
             shape = 17, col = "black",
             size = 4, alpha = 0.5)+
  labs(x = "NMDS1", y = "NMDS2", title = "Fungi")+
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 22))+
  coord_fixed()

fig1c
nmds_legend = get_legend(fig1b)
fig1 = arrangeGrob(fig1b + theme(legend.position = "none"),
                   nmds_legend,
                   fig1c + theme(legend.position = "none"),
                   layout_matrix = rbind(c(1,2,3)),
                   widths = c(4,1.2,4),
                   heights = c(4))
plot(fig1)

#ggsave("./figures/fig1.svg", fig1, width=14, height=8)

# Additional SI plots (S1-3) ====

## (a) Figure S1: NMDS by time ====
figS1a = ggplot()+
  # Dry-down points
  geom_point(data = main_df[main_df$treatment == "field",], 
             aes(x = bac_nmds1, y = bac_nmds2, col = swd), size = 4)+
  geom_path(data = main_df[main_df$treatment == "field",],
            aes(x = bac_nmds1, y = bac_nmds2, col = swd),
            linewidth = 2)+
  scale_colour_gradient(name = "dry-down SWD",
                        low = "#025d96", high = "#c8e7fa")+
  new_scale_colour()+
  # Rewet-up points
  geom_point(data = main_df[main_df$treatment == "drought",],
             aes(x = bac_nmds1, y = bac_nmds2, col = swd), size = 4)+
  geom_path(data = main_df[main_df$treatment == "drought",],
            aes(x = bac_nmds1, y = bac_nmds2, col = swd),
            linewidth = 2)+
  scale_colour_gradient(name = "rewet-up SWD", 
                        low = "#ad0303", high = "#f5d7d7")+
  scale_linetype_manual(values = c(1,2))+
  labs(x = "NMDS1", y = "NMDS2", title = "Prokaryotes")+
  theme(legend.title = element_text(hjust = 0.5),
        plot.title = element_text(size = 20,
                                  hjust = 0.5,
                                  face = "bold"),
        strip.text = element_text(size = 14))+
  coord_fixed()+
  facet_wrap(.~timepoint, ncol = 1)

figS1a

figS1b = ggplot()+
  # Dry-down points
  geom_point(data = main_df[main_df$treatment == "field",], 
             aes(x = fun_nmds1, y = fun_nmds2, col = swd), size = 4)+
  geom_path(data = main_df[main_df$treatment == "field",],
            aes(x = fun_nmds1, y = fun_nmds2, col = swd),
            linewidth = 2)+
  scale_colour_gradient(name = "dry-down SWD",
                        low = "#025d96", high = "#c8e7fa")+
  new_scale_colour()+
  # Rewet-up points
  geom_point(data = main_df[main_df$treatment == "drought",],
             aes(x = fun_nmds1, y = fun_nmds2, col = swd), size = 4)+
  geom_path(data = main_df[main_df$treatment == "drought",],
            aes(x = fun_nmds1, y = fun_nmds2, col = swd),
            linewidth = 2)+
  scale_colour_gradient(name = "rewet-up SWD", 
                        low = "#ad0303", high = "#f5d7d7")+
  scale_linetype_manual(values = c(1,2))+
  labs(x = "NMDS1", y = "NMDS2", title = "Fungi")+
  theme(legend.title = element_text(hjust = 0.5),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 20),
        strip.text = element_text(size = 14))+
  coord_fixed()+
  facet_wrap(.~timepoint, ncol = 1)

figS1 = grid.arrange(figS1a + theme(legend.position = "none"),
                     nmds_legend,
                     figS1b + theme(legend.position = "none"),
                     layout_matrix = rbind(c(1,2,3)),
                     widths = c(10,0.8,10),
                     heights = 12)
#ggsave("./figures/figS1.svg", figS1, width=10, height=10)

## (b) Figure S2/S3: Alternative plots showing response over time with binned SWD levels ====
main_df$swd_level = factor(paste("SWD level", c(rep(c(1:11, 11:1),5), rep(c(11,1,11), each = 3))),
                           levels = paste("SWD level", 1:11))

main_df$timepoint_integer = ifelse(main_df$timepoint == "3 days", 1, 
                                   ifelse(main_df$timepoint ==  "7 days", 2,
                                          ifelse(main_df$timepoint == "14 days", 3,
                                                 ifelse(main_df$timepoint == "35 days", 4,
                                                        ifelse(main_df$timepoint == "70 days", 5, 0)))))

# These plots are imprecise as the SWD levels do not correspond exactly to SWD:
ggplot(main_df[main_df$timepoint_integer != 0,])+
  geom_text(aes(x = timepoint, y = swd, label = swd_level, colour = treatment))
# However, the SWD levels are similar enough across timepoints to  
# provide some useful information about how communities are responding
# within each level

## For legend
figS2_legend = get_legend(ggplot(data = main_df[main_df$treatment %in% c("field", "drought"),], 
                                 aes(x = bac_nmds1, y = bac_nmds2, 
                                     colour = treatment))+
                            geom_point()+
                            scale_colour_manual(name = "initial state",
                                                labels = c("well-watered",
                                                           "severely droughted"), 
                                                values = c("blue", "red"))+
                            theme(legend.text = element_text(size = 12),
                                  legend.title = element_text(size = 16)))

figS2 = ggplot()+
  # Initial points
  geom_point(data = main_df[!main_df$treatment %in% c("field", "drought"),-24],
             aes(x = bac_nmds1, y = bac_nmds2,
                 shape = treatment), colour = "grey")+
  # Dry-down points
  geom_path(data = main_df[main_df$treatment == "field",],
            aes(x = bac_nmds1, y = bac_nmds2), colour = "blue",
            linewidth = 1, alpha = 0.2)+
  geom_text(data = main_df[main_df$treatment == "field",], 
            aes(x = bac_nmds1, y = bac_nmds2,
                label = timepoint_integer), colour = "blue", size = 4)+
  # Rewet-up points
  geom_path(data = main_df[main_df$treatment == "drought",],
            aes(x = bac_nmds1, y = bac_nmds2), colour = "red",
            linewidth = 1, alpha = 0.2)+
  geom_text(data = main_df[main_df$treatment == "drought",], 
            aes(x = bac_nmds1, y = bac_nmds2,
                label = timepoint_integer), colour = "red", size = 4)+
  scale_linetype_manual(values = c(1,2))+
  labs(x = "NMDS1", y = "NMDS2", title = "Prokaryotes")+
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 20),
        strip.text = element_text(size = 14),
        axis.text = element_text(size = 8))+
  guides(shape = "none")+
  coord_fixed()+
  facet_wrap(.~swd_level, nrow = 4)

figS3 = ggplot()+
  # Initial points
  geom_point(data = main_df[!main_df$treatment %in% c("field", "drought"),-24],
             aes(x = fun_nmds1, y = fun_nmds2,
                 shape = treatment), colour = "grey")+
  # Dry-down points
  geom_path(data = main_df[main_df$treatment == "field",],
            aes(x = fun_nmds1, y = fun_nmds2), colour = "blue",
            linewidth = 1, alpha = 0.2)+
  geom_text(data = main_df[main_df$treatment == "field",], 
            aes(x = fun_nmds1, y = fun_nmds2,
                label = timepoint_integer), colour = "blue", size = 4)+
  # Rewet-up points
  geom_path(data = main_df[main_df$treatment == "drought",],
            aes(x = fun_nmds1, y = fun_nmds2), colour = "red",
            linewidth = 1, alpha = 0.2)+
  geom_text(data = main_df[main_df$treatment == "drought",], 
            aes(x = fun_nmds1, y = fun_nmds2,
                label = timepoint_integer), colour = "red", size = 4)+
  scale_linetype_manual(values = c(1,2))+
  labs(x = "NMDS1", y = "NMDS2", title = "Fungi")+
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 20),
        strip.text = element_text(size = 14),
        axis.text = element_text(size = 8))+
  guides(shape = "none")+
  coord_fixed()+
  facet_wrap(.~swd_level, nrow = 4)

figS2_final = grid.arrange(figS2,figS2_legend,nrow = 1,
                           widths = c(4,1))
figS3_final = grid.arrange(figS3, figS2_legend, nrow = 1,
                           widths = c(4,1))

#ggsave("./figures/figS2.svg", figS2_final, width=9, height=9)
#ggsave("./figures/figS3.svg", figS3_final, width=9, height=9)

## (c) Table S2/3: PERMANOVA ====
get_nmds_df = function(my_dist, my_rows){
  as.dist(as.matrix(my_dist)[my_rows, my_rows])
}

# For all five timepoints in gradient phase
adonis2(get_nmds_df(bac_bcavg, 1:110)~swd*treatment*timepoint,
        data = main_df[1:110,], by = "terms")
adonis2(get_nmds_df(fun_bcavg, 1:110)~swd*treatment*timepoint,
        data = main_df[1:110,], by = "terms")

# By timepoints, prokaryotes
adonis2(get_nmds_df(bac_bcavg, 1:22)~swd*treatment,
        data = main_df[1:22,], by = "terms")
adonis2(get_nmds_df(bac_bcavg, 23:44)~swd*treatment,
        data = main_df[23:44,], by = "terms")
adonis2(get_nmds_df(bac_bcavg, 45:66)~swd*treatment,
        data = main_df[45:66,], by = "terms")
adonis2(get_nmds_df(bac_bcavg, 67:88)~swd*treatment,
        data = main_df[67:88,], by = "terms")
adonis2(get_nmds_df(bac_bcavg, 89:110)~swd*treatment,
        data = main_df[89:110,], by = "terms")

# By timepoints, fungi
adonis2(get_nmds_df(fun_bcavg, 1:22)~swd*treatment,
        data = main_df[1:22,], by = "terms")
adonis2(get_nmds_df(fun_bcavg, 23:44)~swd*treatment,
        data = main_df[23:44,], by = "terms")
adonis2(get_nmds_df(fun_bcavg, 45:66)~swd*treatment,
        data = main_df[45:66,], by = "terms")
adonis2(get_nmds_df(fun_bcavg, 67:88)~swd*treatment,
        data = main_df[67:88,], by = "terms")
adonis2(get_nmds_df(fun_bcavg, 89:110)~swd*treatment,
        data = main_df[89:110,], by = "terms")
