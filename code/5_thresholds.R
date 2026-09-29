## Breakpoint analysis
source("./code/initial_setup.R")

library(MuMIn)
library(chngpt)

## We'll fit segmented models for each trajectory. Generally, step, segmented, and hinge
# models seem plausible. We'll fit all of these and compare
get_models = function(my_data, model_type){
  
  stopifnot(model_type %in% c("step", "segmented", "hinge", "stegmented"))
  
  names(my_data) = c("x", "y", "timepoint")
  my_data$timepoint = as.numeric(my_data$timepoint)
  
  # Hardcoding list and loop lengths for ease
  my_modlist = vector("list", length = 6)
  for (tt in 1:6){
    
    if (tt == 6){
      # For all timepoints
      curr_data = my_data[my_data$timepoint %in% c(1:5),]
    } else {
      curr_data = my_data[my_data$timepoint == tt,]
    }
    
    my_modlist[[tt]] = list()
    
    if ("step" %in% model_type){
      step_mod = chngptm(y ~ 1, ~ x, type = "step", family = "gaussian", data = curr_data)
      my_modlist[[tt]] = c(my_modlist[[tt]], step = list(step_mod))
    }
    
    if ("segmented" %in% model_type){
      seg_mod = chngptm(y ~ 1, ~ x, type = "segmented", family = "gaussian", data = curr_data)
      my_modlist[[tt]] = c(my_modlist[[tt]], seg = list(seg_mod))
    }
    
    if ("hinge" %in% model_type){
      
      hinge_mod = chngptm(y ~ 1, ~ x, type = "hinge", family = "gaussian", data = curr_data)
      my_modlist[[tt]] = c(my_modlist[[tt]], hinge = list(hinge_mod))
    }
    
    if ("stegmented" %in% model_type){
      steg_mod = chngptm(y ~ 1, ~ x, type = "stegmented", family = "gaussian", data = curr_data,
                         var.type = "bootstrap", ci.bootstrap.size = 1000)
      my_modlist[[tt]] = c(my_modlist[[tt]], steg = list(steg_mod))
    }
  }
  
  names(my_modlist) = c("t1", "t2", "t3", "t4", "t5", "t_all")
  return(my_modlist)
}

## Limiting search to plausible models only

## Getting thresholds for initial state (more appropriate for finding
## thresholds for bray-curtis):
drydown_funbc = get_models(main_df[main_df$treatment == "field",
                                   c("swd", "fun_bcinitial", "timepoint")],
                           model_type = c("segmented", "hinge"))

lapply(drydown_funbc$t_all, AICc) # Hinge
lapply(drydown_funbc$t5, AICc) # Hinge

rewetup_funbc = get_models(main_df[main_df$treatment == "drought",
                                   c("swd", "fun_bcinitial", "timepoint")],
                           model_type = c("segmented", "hinge"))
lapply(rewetup_funbc$t_all, AICc)
lapply(rewetup_funbc$t5, AICc)

rewetup_bacshannon = get_models(main_df[main_df$treatment == "drought",
                                        c("swd", "bac_shannon", "timepoint")],
                                model_type = c("segmented", "hinge"))
lapply(rewetup_bacshannon$t_all, AICc) # Hinge
lapply(rewetup_bacshannon$t5, AICc) # Hinge

rewetup_mbc = get_models(main_df[main_df$treatment == "drought",
                                 c("swd", "mbc", "timepoint")],
                         model_type = c("segmented", "hinge"))
lapply(rewetup_mbc$t_all, AICc) # Seg
lapply(rewetup_mbc$t5, AICc) # Hinge

rewetup_mbn = get_models(main_df[main_df$treatment == "drought",
                                 c("swd", "mbn", "timepoint")],
                         model_type = c("segmented", "hinge"))
lapply(rewetup_mbn$t_all, AICc) # Hinge
lapply(rewetup_mbn$t5, AICc) # Hinge

## Plotting selected models
# We'll plot selected T5 models

plot_segmented = function(curr_model, curr_measure, trajectory, step = FALSE){
  
  temp_df = main_df
  curr_col = which(names(temp_df) == curr_measure)
  names(temp_df)[curr_col] = "curr_measure"
  
  model_df = data.frame(x = temp_df$swd[temp_df$treatment == trajectory &
                                          temp_df$timepoint == "70 days"],
                        y = predict(curr_model))
  
  if (step == TRUE){
    # Adding in another point for the step model to accurately depict the step
    model_df = rbind(model_df, data.frame(x = curr_model$chngpt,
                                          y = coefficients(curr_model)[1]+
                                            coefficients(curr_model)[2]))
  }
  
  p1 = ggplot()+
    geom_line(data = model_df, aes(x = x, y = y),
              colour = mycols[trajectory],
              linewidth = 1.25,
              arrow = arrow(ends = ifelse(trajectory == "field", "last", "first"),
                            type = "closed",
                            length = unit(0.3, "cm")))+
    geom_point(data = temp_df[temp_df$treatment == trajectory,],
               aes(x = swd, y = curr_measure,
                   alpha = ifelse(timepoint == "70 days", "show", "fade")),
               colour = mycols[trajectory],
               size = 4)+
    geom_vline(xintercept = curr_model$chngpt, linetype = "dotted")+
    annotate("rect", 
             xmin = summary(curr_model)$chngpt[3],
             xmax = summary(curr_model)$chngpt[4],
             ymin = -Inf, ymax = Inf, alpha = 0.1, fill = mycols[trajectory])+
    annotate("text", 
             x = Inf, y = Inf, label = paste0("P = ", signif(summary(curr_model)$coefficients[2,5], 3)),
             hjust = 1.1, vjust = 1.3)+
    scale_alpha_manual(values = c('show' = 1, 'fade' = 0.2))+
    guides(alpha = "none")+
    labs(x = "SWD") +
    theme(legend.position = "none",
          axis.title.y = element_text(size = 11),
          axis.text.y = element_text(size = 9),
          axis.text.x = element_text(size = 8),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = -0.3,
                                    size = 16,
                                    face = "bold"))
  
  return(p1)
}

options(scipen = -2)

fig5 = grid.arrange(plot_segmented(drydown_funbc$t5$hinge, "fun_bcinitial", "field")+
                        labs(title = "a", y = "Fungal similarity\nto well-watered state\n(dry-down)"),
                      plot_segmented(rewetup_funbc$t5$hinge, "fun_bcinitial", "drought")+
                        labs(title = "b", y = "Fungal similarity\nto severe drought state\n(rewet-up)"),
                      plot_segmented(rewetup_bacshannon$t5$hinge, "bac_shannon", "drought")+
                        labs(title = "c", y = "\nProkaryotic Shannon\n(rewet-up)"),
                      plot_segmented(rewetup_mbc$t5$hinge, "mbc", "drought")+
                        labs(title = "d", y = expression(atop(paste("Microbial C (\U00B5g C ", g^1, "dwt)"), 
                                                 "(rewet-up)")))+
                        scale_y_continuous(labels=function(x)x*1000),
                      plot_segmented(rewetup_mbn$t5$hinge, "mbn", "drought")+
                        labs(title = "e", y = expression(atop(paste("Microbial N (\U00B5g N ", g^1, "dwt)"), 
                                                 "(rewet-up)")))+
                        scale_y_continuous(labels=function(x)x*1000))
#ggsave("./figures/fig5.svg", fig5, width=7, height=8)
