# 4. GAMs for all other measures ====
# gam.check and other diagnostics have been run previously to check for k and distrbutions

source("./code/initial_setup.R")

library(marginaleffects)

## Creating dataframe to store predicted values
temp_df = droplevels(main_df[1:110,])
predicted_data = data.frame(treatment = rep(rep(levels(temp_df$treatment), each = 1000),5),
                            swd = rep(seq(min(temp_df$swd), max(temp_df$swd), length.out = 1000), 2*5),
                            timepoint = rep(levels(temp_df$timepoint), each = 2*1000),
                            stringsAsFactors = TRUE)

## (a) Prokaryotes Shannon ====
bac_shannon_gam = gam(bac_shannon ~ s(swd, k = 5) +
                        s(swd, treatment, bs = "sz") +
                        s(swd, timepoint, bs = "sz") +
                        s(swd, timepoint, treatment, bs = "sz"),
                      method = "REML",
                      select = TRUE,
                      data = main_df[1:110,])

## Quadratic model
bac_shannon_poly = lm(bac_shannon~poly(swd,2)*treatment*timepoint, 
                      data = main_df[1:110,])

predicted_data$bac_shannon_gam = get_gam_ci(bac_shannon_gam, predicted_data)
predicted_data$bac_shannon_poly = as.data.frame(predict(bac_shannon_poly, predicted_data, interval = "confidence"))

## (b) Fungal Shannon ====
fun_shannon_gam = gam(fun_shannon ~ s(swd, k = 5) +
                        s(swd, treatment, bs = "sz") +
                        s(swd, timepoint, bs = "sz") +
                        s(swd, timepoint, treatment, bs = "sz"),
                      method = "REML",
                      select = TRUE,
                      data = main_df[1:110,])

fun_shannon_poly = lm(fun_shannon~poly(swd,2)*treatment*timepoint, 
                      data = main_df[1:110,])
predicted_data$fun_shannon_gam = get_gam_ci(fun_shannon_gam, predicted_data)
predicted_data$fun_shannon_poly = as.data.frame(predict(fun_shannon_poly, predicted_data, 
                                                        interval = "confidence"))

## (c) MBC ====
mbc_gam = gam(mbc ~ s(swd, k = 5) +
                s(swd, treatment, bs = "sz") +
                s(swd, timepoint, bs = "sz") +
                s(swd, timepoint, treatment, bs = "sz"),
              method = "REML",
              select = TRUE,
              data = main_df[1:110,])

# Polynomial
mbc_poly = lm(mbc~poly(swd,2)*treatment*timepoint, 
              data = main_df[1:110,])
predicted_data$mbc_gam = get_gam_ci(mbc_gam, predicted_data)
predicted_data$mbc_poly = as.data.frame(predict(mbc_poly, predicted_data, 
                                                interval = "confidence"))

## (d) MBN ====
mbn_gam = gam(mbn ~ s(swd, k = 4) + 
                s(swd, treatment, bs = "sz") +
                s(swd, timepoint, bs = "sz") +
                s(swd, treatment, timepoint, bs = "sz"),
              method = "REML",
              select = TRUE,
              data = main_df[1:110,])

mbn_poly = lm(mbn~poly(swd,2)*treatment*timepoint, data = main_df[1:110,])
predicted_data$mbn_gam = get_gam_ci(mbn_gam, predicted_data)
predicted_data$mbn_poly = as.data.frame(predict(mbn_poly, predicted_data, 
                                                interval = "confidence"))

## (e) DOC ====
doc_gam = gam(doc ~ s(swd, k = 6) + 
                s(swd, treatment, bs = "sz") +
                s(swd, timepoint, bs = "sz") +
                s(swd, treatment, timepoint, bs = "sz"),
              method = "REML",
              select = TRUE,
              data = main_df[1:110,])

doc_poly = lm(doc~poly(swd,2)*treatment*timepoint, data = main_df[1:110,])
predicted_data$doc_gam = get_gam_ci(doc_gam, predicted_data)
predicted_data$doc_poly = as.data.frame(predict(doc_poly, predicted_data, 
                                                interval = "confidence"))

## (f) DON ====
don_gam = gam(don ~ s(swd, k = 6) + 
                s(swd, treatment, bs = "sz") +
                s(swd, timepoint, bs = "sz") +
                s(swd, treatment, timepoint, bs = "sz"),
              family = tw,
              method = "REML",
              select = TRUE,
              data = main_df[1:110,])

don_poly = lm(don~poly(swd,2)*treatment*timepoint, data = main_df[1:110,])
don_predicted = get_gam_ci(don_gam, predicted_data)

# Backtransforming from Tweedie distribution
don_predicted$fit = don_gam$family$linkinv(don_predicted$fit)
don_predicted$point_hi = don_gam$family$linkinv(don_predicted$point_hi)
don_predicted$point_lo = don_gam$family$linkinv(don_predicted$point_lo)
don_predicted$sim_hi = don_gam$family$linkinv(don_predicted$sim_hi)
don_predicted$sim_lo = don_gam$family$linkinv(don_predicted$sim_lo)

predicted_data$don_gam = don_predicted
predicted_data$don_poly = as.data.frame(predict(don_poly, predicted_data, 
                                                interval = "confidence"))

## (g) Prokaryotes richness (SI) ====
bac_richness_gam = gam(bac_richness ~ s(swd, k = 5) +
                         s(swd, treatment, bs = "sz") +
                         s(swd, timepoint, bs = "sz") +
                         s(swd, timepoint, treatment, bs = "sz"),
                       method = "REML",
                       select = TRUE,
                       data = main_df[1:110,])

bac_richness_poly = lm(bac_richness~poly(swd,2)*treatment*timepoint, 
                       data = main_df[1:110,])
predicted_data$bac_richness_gam = get_gam_ci(bac_richness_gam, predicted_data)
predicted_data$bac_richness_poly = as.data.frame(predict(bac_richness_poly, predicted_data, 
                                                         interval = "confidence"))

## (h) Fungal richness (SI) ====
fun_richness_gam = gam(fun_richness ~ s(swd, k = 5) +
                         s(swd, treatment, bs = "sz") +
                         s(swd, timepoint, bs = "sz") +
                         s(swd, timepoint, treatment, bs = "sz"),
                       method = "REML",
                       select = TRUE,
                       data = main_df[1:110,])

# Polynomial model
fun_richness_poly = lm(fun_richness~poly(swd,2)*treatment*timepoint, 
                       data = main_df[1:110,])
predicted_data$fun_richness_gam = get_gam_ci(fun_richness_gam, predicted_data)
predicted_data$fun_richness_poly = as.data.frame(predict(fun_richness_poly, predicted_data, 
                                                         interval = "confidence"))

## (i) BC to controls prokaryotes (SI) ====
bac_bccontrol_gam = gam(bac_bccontrol ~ s(swd, k = 6) +
                          s(swd, treatment, bs = "sz") +
                          s(swd, timepoint, bs = "sz") +
                          s(swd, timepoint, treatment, bs = "sz"),
                        method = "REML",
                        select = TRUE,
                        data = main_df[1:110,])

# Polynomial model
bac_bccontrol_poly = lm(bac_bccontrol~poly(swd,2)*treatment*timepoint, 
                        data = main_df[1:110,])
predicted_data$bac_bccontrol_gam = get_gam_ci(bac_bccontrol_gam, predicted_data)
predicted_data$bac_bccontrol_poly = as.data.frame(predict(bac_bccontrol_poly, predicted_data, 
                                                          interval = "confidence"))

## (j) BC to controls fungi (SI) ====
fun_bccontrol_gam = gam(fun_bccontrol ~ s(swd, k = 5) +
                          s(swd, treatment, bs = "sz") +
                          s(swd, timepoint, bs = "sz") +
                          s(swd, timepoint, treatment, bs = "sz"),
                        method = "REML",
                        select = TRUE,
                        data = main_df[1:110,])

# Polynomial model
fun_bccontrol_poly = lm(fun_bccontrol~poly(swd,2)*treatment*timepoint, 
                        data = main_df[1:110,])
predicted_data$fun_bccontrol_gam = get_gam_ci(fun_bccontrol_gam, predicted_data)
predicted_data$fun_bccontrol_poly = as.data.frame(predict(fun_bccontrol_poly, predicted_data, 
                                                          interval = "confidence"))

## (k) BC to initial prokaryotes (SI) ====
bac_bcinitial_gam = gam(bac_bcinitial ~ s(swd, k = 6) +
                          s(swd, treatment, bs = "sz") +
                          s(swd, timepoint, bs = "sz") +
                          s(swd, timepoint, treatment, bs = "sz"),
                        method = "REML",
                        select = TRUE,
                        data = main_df[1:110,])

# Polynomial model
bac_bcinitial_poly = lm(bac_bcinitial~poly(swd,2)*treatment*timepoint, 
                        data = main_df[1:110,])
predicted_data$bac_bcinitial_gam = get_gam_ci(bac_bcinitial_gam, predicted_data)
predicted_data$bac_bcinitial_poly = as.data.frame(predict(bac_bcinitial_poly, predicted_data, 
                                                          interval = "confidence"))

## (l) BC to initial fungi (SI) ====
fun_bcinitial_gam = gam(fun_bcinitial ~ s(swd, k = 5) +
                          s(swd, treatment, bs = "sz") +
                          s(swd, timepoint, bs = "sz") +
                          s(swd, timepoint, treatment, bs = "sz"),
                        method = "REML",
                        select = TRUE,
                        data = main_df[1:110,])

# Polynomial model
fun_bcinitial_poly = lm(fun_bcinitial~poly(swd,2)*treatment*timepoint, 
                        data = main_df[1:110,])
predicted_data$fun_bcinitial_gam = get_gam_ci(fun_bcinitial_gam, predicted_data)
predicted_data$fun_bcinitial_poly = as.data.frame(predict(fun_bcinitial_poly, predicted_data, 
                                                          interval = "confidence"))

# Fig. 4: GAM plots (T1 and T5) ====
str(predicted_data)
predicted_data$timepoint = factor(predicted_data$timepoint,
                                  levels = c("3 days",
                                             "7 days",
                                             "14 days",
                                             "35 days",
                                             "70 days"))

predicted_data$treatment = factor(predicted_data$treatment,
                                  levels = c("field",
                                             "drought"))

## For generating plots
plot_GAM = function(my_model, my_measure, my_ylab, my_title, my_timepoints){
  
  # my_model = fitted GAM or polynomial model
  # my_measure = selected measure (pulls from main_df)
  # my_ylab = y label
  # my_title = alphabet title for subplots
  # my_timepoints = selected timepoints to show
  
  curr_data = data.frame(predicted_data[,colnames(predicted_data) == my_model])
  curr_data = cbind(curr_data, predicted_data[c(1:3)])
  
  curr_plot = ggplot(data = curr_data[curr_data$timepoint %in% my_timepoints,])+
    # 95% CI
    geom_ribbon(aes(x = swd, ymin = sim_lo, 
                    ymax = sim_hi,
                    fill = treatment), alpha = 0.3)+
    # Points
    geom_point(data = main_df[main_df$timepoint %in% my_timepoints,],
               aes(x = swd, y = eval(parse(text = my_measure)), colour = treatment),
               size = 4, alpha = 0.5)+
    # Dry-down arrow
    geom_line(data = curr_data[curr_data$timepoint %in% my_timepoints &
                                 curr_data$treatment == "drought",],
              aes(x = swd, y = fit), colour = mycols[2],
              linewidth = 2, 
              arrow = arrow(ends = "first", type = "open", 
                            length = unit(0.3, "cm")))+
    # Rewet-up arrow
    geom_line(data = curr_data[curr_data$timepoint %in% my_timepoints &
                                 curr_data$treatment == "field",],
              aes(x = swd, y = fit), colour = mycols[1],
              linewidth = 2, 
              arrow = arrow(ends = "last", type = "closed", 
                            length = unit(0.3, "cm")))+
    facet_wrap(.~timepoint, nrow = 1)+
    labs(x = "SWD", 
         y = my_ylab,
         title = my_title)+ 
    scale_colour_manual(name = "",
                        values = mycols,
                        labels = c("dry-down",
                                   "rewet-up"))+
    scale_fill_manual(name = "",
                      values = mycols,
                      labels = c("dry-down",
                                 "rewet-up"))+
    theme(legend.position = "none",
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 11),
          axis.text.x = element_text(size = 10),
          strip.text = element_text(size = 13),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = -0.1,
                                    size = 24,
                                    face = "bold"))
  return(curr_plot)
}

fig4a = plot_GAM("bac_shannon_gam", "bac_shannon", "Prokaryotic Shannon", "a", 
         levels(main_df$timepoint)[c(1,5)])
fig4b = plot_GAM("fun_shannon_gam", "fun_shannon", "Fungal Shannon", "b", 
         levels(main_df$timepoint)[c(1,5)])
fig4c = plot_GAM("mbc_gam", "mbc",
         expression(paste("Microbial C ", 
                          "(\U00B5g C ", 
                          g^-1, 
                          "soil)")),
         "c",
         levels(main_df$timepoint)[c(1,5)])+
  scale_y_continuous(labels=function(x)x*1000)

fig4d = plot_GAM("mbn_gam", "mbn",
         expression(paste("Microbial N ", 
                          "(\U00B5g N ", 
                          g^-1, 
                          "soil)")),
         "d",
         levels(main_df$timepoint)[c(1,5)])+
  scale_y_continuous(labels=function(x)x*1000)

fig4e = plot_GAM("doc_gam", "doc",
         expression(paste("DOC ", 
                          "(\U00B5g C ", 
                          g^-1, 
                          "soil)")),
         "e",
         levels(main_df$timepoint)[c(1,5)])+
  scale_y_continuous(labels=function(x)x*1000)

fig4f = plot_GAM("don_gam", "don",
         expression(paste("DON ", 
                          "(\U00B5g N ", 
                          g^-1, 
                          "soil)")),
         "f",
         levels(main_df$timepoint)[c(1,5)])+
  scale_y_continuous(labels=function(x)x*1000)+
  facet_wrap(.~timepoint, scales = "free")

gam_legend = get_legend(fig4a + 
                          theme(legend.position = "top",
                                legend.text = element_text(size = 20)))
text1 <- data.frame(label = c("wetter", ""),
                    timepoint   = c("3 days", "70 days"))
text2 <- data.frame(label = c("drier", ""),
                    timepoint   = c("3 days", "70 days"))

x1 = range(predicted_data$swd)[1] + 
  (range(predicted_data$swd)[2]-range(predicted_data$swd)[1])*0.35
x2 = range(predicted_data$swd)[1] + 
  (range(predicted_data$swd)[2]-range(predicted_data$swd)[1])*0.65

fig4a_annotated = fig4a +
  geom_text(data = text1,
            aes(x = x1, y = 5.45, label = label),
            size = 5, hjust = 1.2, vjust = 0.3)+
  geom_text(data = text2,
            aes(x = x2, y = 5.45, label = label),
            size = 5, hjust = -0.35, vjust = 0.3)+
  geom_segment(data = data.frame(timepoint = c("3 days")), 
               x = x1, xend = x2, y = 5.45, yend=  5.45,
               arrow = arrow(length = unit(6, "pt"),
                             ends = "both",
                             type = "closed"),
               linewidth = 1.1)

fig4 = grid.arrange(gam_legend, 
                    fig4a_annotated, fig4b, fig4c, fig4d, fig4e, fig4f,
                    layout_matrix = rbind(c(1,NA,1),
                                          c(2,NA,3),
                                          c(4,NA,5),
                                          c(6,NA,7)),
                    heights = c(0.3,1,1,1),
                    widths = c(1, 0.05, 1))
#ggsave("./figures/fig4.svg", fig4, width = 10, height = 12)

# Additional SI plots ====

## (a) Fig. S9: All timepoint GAMs ====

# For deviance explained labels
bac_bc_text = data.frame(label = paste0("D[exp] == '", 
                                        round(summary(bac_bccontrol_gam)$dev.expl*100, 1), "%'"),
                         timepoint = factor("70 days"))
fun_bc_text = data.frame(label = paste0("D[exp] == '", 
                                        round(summary(fun_bccontrol_gam)$dev.expl*100, 1), "%'"),
                         timepoint = factor("70 days"))
bac_bcinit_text = data.frame(label = paste0("D[exp] == '", 
                                        round(summary(bac_bcinitial_gam)$dev.expl*100, 1), "%'"),
                         timepoint = factor("70 days"))
fun_bcinit_text = data.frame(label = paste0("D[exp] == '", 
                                        round(summary(fun_bcinitial_gam)$dev.expl*100, 1), "%'"),
                         timepoint = factor("70 days"))
bac_shan_text = data.frame(label = paste0("D[exp] == '", 
                                          round(summary(bac_shannon_gam)$dev.expl*100, 1), "%'"),
                           timepoint = factor("70 days"))
bac_rich_text = data.frame(label = paste0("D[exp] == '", 
                                          round(summary(bac_richness_gam)$dev.expl*100, 1), "%'"),
                           timepoint = factor("70 days"))
fun_shan_text = data.frame(label = paste0("D[exp] == '", 
                                          round(summary(fun_shannon_gam)$dev.expl*100, 1), "%'"),
                           timepoint = factor("70 days"))
fun_rich_text = data.frame(label = paste0("D[exp] == '", 
                                          round(summary(fun_richness_gam)$dev.expl*100, 1), "%'"),
                           timepoint = factor("70 days"))
mbc_text = data.frame(label = paste0("D[exp] == '", 
                                     round(summary(mbc_gam)$dev.expl*100, 1), "%'"),
                      timepoint = factor("70 days"))
mbn_text = data.frame(label = paste0("D[exp] == '", 
                                     round(summary(mbn_gam)$dev.expl*100, 1), "%'"),
                      timepoint = factor("70 days"))
doc_text = data.frame(label = paste0("D[exp] == '", 
                                     round(summary(doc_gam)$dev.expl*100, 1), "%'"),
                      timepoint = factor("70 days"))
don_text = data.frame(label = paste0("D[exp] == '", 
                                     round(summary(don_gam)$dev.expl*100, 1), "%'"),
                      timepoint = factor("70 days"))

figS9a = plot_GAM("bac_bccontrol_gam", "bac_bccontrol", 
                  "Prokaryotic similarity\n(to original community)", "", 
                  levels(main_df$timepoint)[c(1:5)])+
  geom_text(data = bac_bc_text, aes(label = label),
            x = Inf, y = -Inf, parse = TRUE,
            vjust = -0.5, hjust = 1.1, size = 4)

figS9b = plot_GAM("fun_bccontrol_gam", "fun_bccontrol", 
                  "Fungal similarity\n(to original community)", "", 
                  levels(main_df$timepoint)[c(1:5)])+
  geom_text(data = fun_bc_text, aes(label = label),
            x = Inf, y = -Inf, parse = TRUE,
            vjust = -0.5, hjust = 1.1, size = 4)

figS9c = plot_GAM("bac_bcinitial_gam", "bac_bcinitial", 
                  "Prokaryotic similarity\n(to initial state)", "", 
                  levels(main_df$timepoint)[c(1:5)])+
  geom_text(data = bac_bcinit_text, aes(label = label),
            x = Inf, y = -Inf, parse = TRUE,
            vjust = -0.5, hjust = 1.1, size = 4)

figS9d = plot_GAM("fun_bcinitial_gam", "fun_bcinitial", 
                  "Fungal similarity\n(to initial state)", "", 
                  levels(main_df$timepoint)[c(1:5)])+
  geom_text(data = fun_bcinit_text, aes(label = label),
            x = Inf, y = -Inf, parse = TRUE,
            vjust = -0.5, hjust = 1.1, size = 4)

figS9e = plot_GAM("bac_shannon_gam", "bac_shannon", "Prokaryotic Shannon", "", 
                 levels(main_df$timepoint)[c(1:5)])+
  geom_text(data = bac_shan_text, aes(label = label),
            x = Inf, y = Inf, parse = TRUE,
            vjust = 1.5, hjust = 1.1, size = 4)

figS9f = plot_GAM("fun_shannon_gam", "fun_shannon", "Fungal Shannon", "", 
                  levels(main_df$timepoint)[c(1:5)])+
  geom_text(data = fun_shan_text, aes(label = label),
            x = Inf, y = -Inf, parse = TRUE,
            vjust = -0.5, hjust = 1.1, size = 4)

figS9g = plot_GAM("bac_richness_gam", "bac_richness", "Prokaryotic richness", "", 
                  levels(main_df$timepoint)[c(1:5)])+
  geom_text(data = bac_rich_text, aes(label = label),
            x = Inf, y = Inf, parse = TRUE,
            vjust = 1.5, hjust = 1.1, size = 4)

figS9h = plot_GAM("fun_richness_gam", "fun_richness", "Fungal richness", "", 
                  levels(main_df$timepoint)[c(1:5)])+
  geom_text(data = fun_rich_text, aes(label = label),
            x = Inf, y = -Inf, parse = TRUE,
            vjust = -0.5, hjust = 1.1, size = 4)

figS9i = plot_GAM("mbc_gam", "mbc",
                  expression(paste("Microbial C ", 
                                   "(\U00B5g C ", 
                                   g^-1, 
                                   "soil)")),
                  "",
                  levels(main_df$timepoint)[c(1:5)])+
  scale_y_continuous(labels=function(x)x*1000)+
  geom_text(data = mbc_text, aes(label = label),
            x = Inf, y = Inf, parse = TRUE,
            vjust = 1.5, hjust = 1.1, size = 4)

figS9j = plot_GAM("mbn_gam", "mbn",
                 expression(paste("Microbial N ", 
                                  "(\U00B5g N ", 
                                  g^-1, 
                                  "soil)")),
                 "",
                 levels(main_df$timepoint)[c(1:5)])+
  scale_y_continuous(labels=function(x)x*1000)+
  geom_text(data = mbn_text, aes(label = label),
            x = Inf, y = Inf, parse = TRUE,
            vjust = 1.5, hjust = 1.1, size = 4)

figS9k = plot_GAM("doc_gam", "doc",
                 expression(paste("DOC ", 
                                  "(\U00B5g C ", 
                                  g^-1, 
                                  "soil)")),
                 "",
                 levels(main_df$timepoint)[c(1:5)])+
  scale_y_continuous(labels=function(x)x*1000)+
  geom_text(data = doc_text, aes(label = label),
            x = Inf, y = Inf, parse = TRUE,
            vjust = 1.5, hjust = 1.1, size = 4)


figS9l = plot_GAM("don_gam", "don",
                 expression(paste("DON ", 
                                  "(\U00B5g N ", 
                                  g^-1, 
                                  "soil)")),
                 "",
                 levels(main_df$timepoint)[c(1:5)])+
  scale_y_continuous(labels=function(x)x*1000)+
  facet_wrap(.~timepoint, scales = "free", nrow = 1)+
  geom_text(data = don_text, aes(label = label),
            x = Inf, y = Inf, parse = TRUE,
            vjust = 1.5, hjust = 1.1, size = 4)

figS9 = grid.arrange(gam_legend,
                     figS9a, figS9b,
                     figS9c, figS9d,
                     figS9e, figS9f,
                     figS9g, figS9h,
                     figS9i, figS9j,
                     figS9k, figS9l,
                     layout_matrix = rbind(c(1,NA,1),
                                           c(2,NA,3),
                                           c(4,NA,5),
                                           c(6,NA,7),
                                           c(8,NA,9),
                                           c(10,NA,11),
                                           c(12,NA,13)),
                     heights = c(0.3,1,1,1,1,1,1),
                     widths = c(1, 0.05, 1))
#ggsave("figures/figS9.svg", figS9, width=20, height=18)

## (b) Fig. S10: GAM slopes ====

## Getting slopes
get_slopes = function(x){
  plot_slopes(eval(parse(text = x)), variables = c('swd'),
              condition = c('swd', 'treatment', 'timepoint'),
              type = 'link',
              draw = FALSE)
}

bac_bc_slopes = get_slopes("bac_bccontrol_gam")
fun_bc_slopes = get_slopes("fun_bccontrol_gam")
bac_bcinit_slopes = get_slopes("bac_bcinitial_gam")
fun_bcinit_slopes = get_slopes("fun_bcinitial_gam")
bac_shan_slopes = get_slopes("bac_shannon_gam")
fun_shan_slopes = get_slopes("fun_shannon_gam")
bac_rich_slopes = get_slopes("bac_richness_gam")
fun_rich_slopes = get_slopes("fun_richness_gam")
mbc_slopes = get_slopes("mbc_gam")
mbn_slopes = get_slopes("mbn_gam")
doc_slopes = get_slopes("doc_gam")
don_slopes = get_slopes("don_gam")

## Plotting
plotting_slopes = function(x, my_ylab){
  curr_plot = ggplot(data = x)+
    geom_ribbon(aes(x = swd, ymin = conf.low, ymax = conf.high, fill = treatment), alpha = 0.3)+
    geom_line(aes(x = swd, y = estimate, col = treatment), linewidth = 2)+
    geom_hline(yintercept = 0, linetype = "dashed", 
               linewidth = 0.5, col = "black")+
    facet_wrap(.~timepoint, nrow = 1)+
    labs(x = "SWD", y = my_ylab)+
    scale_colour_manual(labels = c("dry-down", "rewet-up"), 
                        values = mycols)+
    scale_fill_manual(labels = c("dry-down", "rewet-up"),
                      values = mycols)+
    theme(legend.position = "none",
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 11),
          axis.text.x = element_text(size = 10),
          strip.text = element_text(size = 13),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = -0.1,
                                    size = 24,
                                    face = "bold"))
  return(curr_plot)
}

bac_bc_slopeplot = plotting_slopes(bac_bc_slopes, "Prokaryotic similarity\n(to original community)")
fun_bc_slopeplot = plotting_slopes(fun_bc_slopes, "Fungal similarity\n(to original community)")
bac_bcinit_slopeplot = plotting_slopes(bac_bcinit_slopes, "Prokaryotic similarity\n(to initial state)")
fun_bcinit_slopeplot = plotting_slopes(fun_bcinit_slopes, "Fungal similarity\n(to initial state)")
bac_shan_slopeplot = plotting_slopes(bac_shan_slopes, "Prokaryotic Shannon")
fun_shan_slopeplot = plotting_slopes(fun_shan_slopes, "Fungal Shannon")
bac_rich_slopeplot = plotting_slopes(bac_rich_slopes, "Prokaryotic richness")
fun_rich_slopeplot = plotting_slopes(fun_rich_slopes, "Fungal richness")
mbc_slopeplot = plotting_slopes(mbc_slopes, "Microbial C")
mbn_slopeplot = plotting_slopes(mbn_slopes, "Microbial N")
doc_slopeplot = plotting_slopes(doc_slopes, "DOC")
don_slopeplot = plotting_slopes(don_slopes, "DON")

figS10 = grid.arrange(gam_legend, 
                      bac_bc_slopeplot,
                      fun_bc_slopeplot, 
                      bac_bcinit_slopeplot,
                      fun_bcinit_slopeplot, 
                      bac_shan_slopeplot,
                      fun_shan_slopeplot, 
                      bac_rich_slopeplot,
                      fun_rich_slopeplot,
                      mbc_slopeplot, 
                      mbn_slopeplot,
                      doc_slopeplot, 
                      don_slopeplot,
                      layout_matrix = rbind(c(1,NA,1),
                                            c(2,NA,3),
                                            c(4,NA,5),
                                            c(6,NA,7),
                                            c(8,NA,9),
                                            c(10,NA,11),
                                            c(12,NA,13)),
                      heights = c(0.3,1,1,1,1,1,1),
                      widths = c(1, 0.05, 1))

#ggsave("./figures/figS10.svg", figS10, width=20, height=18)

## (c) Fig. S11: Polynomial plots ====

# Getting R2 values for plots
bac_bc_polytext = data.frame(label = paste0('R["adj"]^2 ==', 
                                            round(summary(bac_bccontrol_poly)$r.sq,3)),
                             timepoint = factor("70 days"))
fun_bc_polytext = data.frame(label = paste0('R["adj"]^2 ==', 
                                            round(summary(fun_bccontrol_poly)$r.sq,3)),
                             timepoint = factor("70 days"))

bac_bcinit_polytext = data.frame(label = paste0('R["adj"]^2 ==', 
                                            round(summary(bac_bcinitial_poly)$r.sq,3)),
                             timepoint = factor("70 days"))
fun_bcinit_polytext = data.frame(label = paste0('R["adj"]^2 ==', 
                                            round(summary(fun_bcinitial_poly)$r.sq,3)),
                             timepoint = factor("70 days"))

bac_shan_polytext = data.frame(label = paste0('R["adj"]^2 ==', 
                                              round(summary(bac_shannon_poly)$r.sq,3)),
                               timepoint = factor("70 days"))
fun_shan_polytext = data.frame(label = paste0('R["adj"]^2 ==', 
                                              round(summary(fun_shannon_poly)$r.sq,3)),
                               timepoint = factor("70 days"))

bac_rich_polytext = data.frame(label = paste0('R["adj"]^2 ==', 
                                              round(summary(bac_richness_poly)$r.sq,3)),
                               timepoint = factor("70 days"))
fun_rich_polytext = data.frame(label = paste0('R["adj"]^2 ==', 
                                              round(summary(fun_richness_poly)$r.sq,3)),
                               timepoint = factor("70 days"))

mbc_polytext = data.frame(label = paste0('R["adj"]^2 ==', 
                                         round(summary(mbc_poly)$r.sq,3)),
                          timepoint = factor("70 days"))
mbn_polytext = data.frame(label = paste0('R["adj"]^2 ==', 
                                         round(summary(mbn_poly)$r.sq,3)),
                          timepoint = factor("70 days"))
doc_polytext = data.frame(label = paste0('R["adj"]^2 ==', 
                                         round(summary(doc_poly)$r.sq,3)),
                          timepoint = factor("70 days"))
don_polytext = data.frame(label = paste0('R["adj"]^2 ==', 
                                         round(summary(don_poly)$r.sq,3)),
                          timepoint = factor("70 days"))

## For generating plots
plot_poly = function(my_model, my_measure, my_ylab, my_title, my_timepoints){
  
  # my_model = fitted GAM or polynomial model
  # my_measure = selected measure (pulls from main_df)
  # my_ylab = y label
  # my_title = alphabet title for subplots
  # my_timepoints = selected timepoints to show
  
  curr_data = data.frame(predicted_data[,colnames(predicted_data) == my_model])
  curr_data = cbind(curr_data, predicted_data[c(1:3)])
  
  curr_plot = ggplot(data = curr_data[curr_data$timepoint %in% my_timepoints,])+
    # 95% CI
    geom_ribbon(aes(x = swd, ymin = lwr, 
                    ymax = upr,
                    fill = treatment), alpha = 0.3)+
    # Points
    geom_point(data = main_df[main_df$timepoint %in% my_timepoints,],
               aes(x = swd, y = eval(parse(text = my_measure)), colour = treatment),
               size = 4, alpha = 0.5)+
    # Dry-down arrow
    geom_line(data = curr_data[curr_data$timepoint %in% my_timepoints &
                                 curr_data$treatment == "drought",],
              aes(x = swd, y = fit), colour = mycols[2],
              linewidth = 2, 
              arrow = arrow(ends = "first", type = "open", 
                            length = unit(0.3, "cm")))+
    # Rewet-up arrow
    geom_line(data = curr_data[curr_data$timepoint %in% my_timepoints &
                                 curr_data$treatment == "field",],
              aes(x = swd, y = fit), colour = mycols[1],
              linewidth = 2, 
              arrow = arrow(ends = "last", type = "closed", 
                            length = unit(0.3, "cm")))+
    facet_wrap(.~timepoint, nrow = 1)+
    labs(x = "SWD", 
         y = my_ylab,
         title = my_title)+ 
    scale_colour_manual(name = "",
                        values = mycols,
                        labels = c("dry-down",
                                   "rewet-up"))+
    scale_fill_manual(name = "",
                      values = mycols,
                      labels = c("dry-down",
                                 "rewet-up"))+
    theme(legend.position = "none",
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 11),
          axis.text.x = element_text(size = 10),
          strip.text = element_text(size = 13),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = -0.1,
                                    size = 24,
                                    face = "bold"))
  return(curr_plot)
}


figS11a = plot_poly("bac_bccontrol_poly", "bac_bccontrol", 
                  "Prokaryotic similarity\n(to original community)", "", 
                  levels(main_df$timepoint)[c(1:5)])+
  geom_text(data = bac_bc_polytext, aes(label = label),
            x = Inf, y = -Inf, parse = TRUE,
            vjust = -0.5, hjust = 1.1, size = 4)

figS11b = plot_poly("fun_bccontrol_poly", "fun_bccontrol", 
                  "Fungal similarity\n(to original community)", "", 
                  levels(main_df$timepoint)[c(1:5)])+
  geom_text(data = fun_bc_polytext, aes(label = label),
            x = Inf, y = Inf, parse = TRUE,
            vjust = 1.5, hjust = 1.1, size = 4)

figS11c = plot_poly("bac_bcinitial_poly", "bac_bcinitial", 
                    "Prokaryotic similarity\n(to initial state)", "", 
                    levels(main_df$timepoint)[c(1:5)])+
  geom_text(data = bac_bcinit_polytext, aes(label = label),
            x = Inf, y = -Inf, parse = TRUE,
            vjust = -0.5, hjust = 1.1, size = 4)

figS11d = plot_poly("fun_bcinitial_poly", "fun_bcinitial", 
                    "Fungal similarity\n(to initial state)", "", 
                    levels(main_df$timepoint)[c(1:5)])+
  geom_text(data = fun_bcinit_polytext, aes(label = label),
            x = Inf, y = Inf, parse = TRUE,
            vjust = 1.5, hjust = 1.1, size = 4)

figS11e = plot_poly("bac_shannon_poly", "bac_shannon", "Prokaryotic Shannon", "", 
                  levels(main_df$timepoint)[c(1:5)])+
  geom_text(data = bac_shan_polytext, aes(label = label),
            x = Inf, y = Inf, parse = TRUE,
            vjust = 1.5, hjust = 1.1, size = 4)

figS11f = plot_poly("fun_shannon_poly", "fun_shannon", "Fungal Shannon", "", 
                  levels(main_df$timepoint)[c(1:5)])+
  geom_text(data = fun_shan_polytext, aes(label = label),
            x = Inf, y = -Inf, parse = TRUE,
            vjust = -0.5, hjust = 1.1, size = 4)

figS11g = plot_poly("bac_richness_poly", "bac_richness", "Prokaryotic richness", "", 
                  levels(main_df$timepoint)[c(1:5)])+
  geom_text(data = bac_rich_polytext, aes(label = label),
            x = Inf, y = Inf, parse = TRUE,
            vjust = 1.5, hjust = 1.1, size = 4)

figS11h = plot_poly("fun_richness_poly", "fun_richness", "Fungal richness", "", 
                  levels(main_df$timepoint)[c(1:5)])+
  geom_text(data = fun_rich_polytext, aes(label = label),
            x = Inf, y = -Inf, parse = TRUE,
            vjust = -0.5, hjust = 1.1, size = 4)

figS11i = plot_poly("mbc_poly", "mbc",
                  expression(paste("Microbial C ", 
                                   "(\U00B5g C ", 
                                   g^-1, 
                                   "soil)")),
                  "",
                  levels(main_df$timepoint)[c(1:5)])+
  scale_y_continuous(labels=function(x)x*1000)+
  geom_text(data = mbc_polytext, aes(label = label),
            x = Inf, y = Inf, parse = TRUE,
            vjust = 1.5, hjust = 1.1, size = 4)

figS11j = plot_poly("mbn_poly", "mbn",
                  expression(paste("Microbial N ", 
                                   "(\U00B5g N ", 
                                   g^-1, 
                                   "soil)")),
                  "",
                  levels(main_df$timepoint)[c(1:5)])+
  scale_y_continuous(labels=function(x)x*1000)+
  geom_text(data = mbn_polytext, aes(label = label),
            x = Inf, y = Inf, parse = TRUE,
            vjust = 1.5, hjust = 1.1, size = 4)

figS11k = plot_poly("doc_poly", "doc",
                  expression(paste("DOC ", 
                                   "(\U00B5g C ", 
                                   g^-1, 
                                   "soil)")),
                  "",
                  levels(main_df$timepoint)[c(1:5)])+
  scale_y_continuous(labels=function(x)x*1000)+
  geom_text(data = doc_polytext, aes(label = label),
            x = Inf, y = Inf, parse = TRUE,
            vjust = 1.5, hjust = 1.1, size = 4)


figS11l = plot_poly("don_poly", "don",
                  expression(paste("DON ", 
                                   "(\U00B5g N ", 
                                   g^-1, 
                                   "soil)")),
                  "",
                  levels(main_df$timepoint)[c(1:5)])+
  scale_y_continuous(labels=function(x)x*1000)+
  facet_wrap(.~timepoint, scales = "free", nrow = 1)+
  geom_text(data = don_polytext, aes(label = label),
            x = Inf, y = Inf, parse = TRUE,
            vjust = 1.5, hjust = 1.1, size = 4)

figS11 = grid.arrange(gam_legend,
                      figS11a, figS11b,
                      figS11c, figS11d,
                      figS11e, figS11f,
                      figS11g, figS11h,
                      figS11i, figS11j,
                      figS11k, figS11l,
                      layout_matrix = rbind(c(1,NA,1),
                                            c(2,NA,3),
                                            c(4,NA,5),
                                            c(6,NA,7),
                                            c(8,NA,9),
                                            c(10,NA,11),
                                            c(12,NA,13)),
                      heights = c(0.3,1,1,1,1,1,1),
                      widths = c(1, 0.05, 1))

#ggsave("./figures/figS11.svg", figS11, width=20, height=18)

## Calculating degree of hysteresis
# We'll take the absolute difference between the 1000 values we have
get_hysteresis = function(x){
  curr_list = which(names(predicted_data) == x)
  
  my_df = data.frame()
  df_row = 0
  for (tt in unique(predicted_data$timepoint)){
    df_row = df_row + 1
    abs_diff = abs(predicted_data[[curr_list]]$fit[predicted_data$treatment == "field" &
                                                         predicted_data$timepoint == tt]-
                         predicted_data[[curr_list]]$fit[predicted_data$treatment == "drought" &
                                                           predicted_data$timepoint == tt])

    swd_grid = seq(min(predicted_data$swd), max(predicted_data$swd),
                   length.out = nrow(predicted_data)/10) # hardcoding length
    
    # Numerically integrating difference (trapezoidal)
    diff_val = sum(diff(swd_grid)*(abs_diff[-length(abs_diff)]+abs_diff[-1])/2)
    
    to_bind = data.frame(tt, diff_val, measure = sub("_gam", "", x))
    my_df = rbind(my_df, to_bind)
  }
  my_df
}

hys_vals = rbind(get_hysteresis("mbc_gam"),
                 get_hysteresis("mbn_gam"),
                 get_hysteresis("doc_gam"),
                 get_hysteresis("don_gam"),
                 get_hysteresis("bac_richness_gam"),
                 get_hysteresis("fun_richness_gam"),
                 get_hysteresis("bac_shannon_gam"),
                 get_hysteresis("fun_shannon_gam"),
                 get_hysteresis("bac_bccontrol_gam"),
                 get_hysteresis("fun_bccontrol_gam"))

hys_vals$tt = as.factor(hys_vals$tt)
hys_vals$tt = factor(hys_vals$tt, levels = c("3 days", "7 days", "14 days",
                                             "35 days", "70 days"))

hys_vals$measure = as.factor(hys_vals$measure)
hys_vals$measure = factor(hys_vals$measure, levels = c("bac_bccontrol", "fun_bccontrol",
                                                       "bac_shannon", "fun_shannon", 
                                                       "bac_richness", "fun_richness",
                                                       "mbc", "mbn",
                                                       "doc", "don"))

# Scaling extractant values 
hys_vals$diff_val[hys_vals$measure %in% c("doc", "don", "mbc", "mbn")] =
  hys_vals$diff_val[hys_vals$measure %in% c("doc", "don", "mbc", "mbn")]*1000

hys_vals$labels = hys_vals$measure
levels(hys_vals$labels)
levels(hys_vals$labels) = c("Prokaryotic similarity\n(to original community)",
                            "Fungal similarity\n(to original community)",
                            "Prokaryotic Shannon",
                            "Fungal Shannon",
                            "Prokaryotic richness",
                            "Fungal richness",
                            "Microbial C",
                            "Microbial N",
                            "DOC",
                            "DON")


figS12 = ggplot(hys_vals, aes(x = as.numeric(tt), y = diff_val))+
  geom_line(linewidth = 1, alpha = 0.3)+
  geom_point(size = 2)+
  facet_wrap(.~labels, scales = "free", nrow = 5)+
  lims(y = c(0,NA))+
  labs(x = "days", y = "hysteresis magnitude")+
  scale_x_continuous(labels = sub(" days", "", levels(hys_vals$tt)))

#ggsave("./figures/figS12.svg", figS12, width=6, height=8)
