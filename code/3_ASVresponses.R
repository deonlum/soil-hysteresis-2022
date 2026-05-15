# 3. ASV responses and trees ====

source("./code/initial_setup.R")

library(ape)
library(ggtree) # For plotting trees
library(ggtreeExtra) # Add on functions to ggtree
library(tidytree)
library(ggsankey) # For alluvial plots (Fig. S4)

## Categorising responses ====
get_asv_responses = function(my_abundances, my_swd){
  
  my_libsizes = apply(my_abundances, 1, sum)
  
  # Method adapted from #https://stats.stackexchange.com/questions/31858/recovering-raw-coefficients-and-variances-from-orthogonal-polynomial-regression
  raw_poly = cbind(1, poly(my_swd, 2, raw = TRUE)) # Raw poly
  orth_poly = cbind(1, poly(my_swd, 2)) # Orthogonal poly
  tf_mod = lm(orth_poly ~ raw_poly-1) # Relating the two
  gamma = coef(tf_mod) # values to transform with
  
  stopifnot(all.equal(as.vector(1 + crossprod(raw_poly %*% gamma - orth_poly) - 1), 
                      rep(0, (2+1)^2))) # Double-check all is ok
  
  ## Creating matrix to store results
  results_matrix = matrix(nrow = ncol(my_abundances), ncol = 9)
  colnames(results_matrix) = c("orth_int", "orth_x1", "orth_x2", 
                               "orth_int_p", "orth_x1_p", "orth_x2_p",
                               "raw_int", "raw_x1", "raw_x2")
  skipped = 0 # Counting skips to make sure everything corresponds
  
  # starting loop
  for (i in 1:ncol(my_abundances)){
    temp_df = data.frame(counts = my_abundances[,i]/my_libsizes,
                         swd = my_swd,
                         my_libsizes)
    
    # Skip to next loop if no data available
    # or if fewer than 2 data points
    if(sum(temp_df$counts) == 0 |
       sum(temp_df$counts > 0, na.rm = TRUE) < 2){
      skipped = skipped + 1
      results_matrix[i,] = NA
      next} 
    
    # Including variance structure for library size differences
    vf_mod = varConstProp(form = ~my_libsizes)
    
    # Using orthogonal polynomials for numerical stability
    temp_mod = gls(counts~poly(swd, 2), temp_df,
                   weights = vf_mod)
    
    # Storing out values
    orth_ests = summary(temp_mod)$tTable[,1] # Orth estimates
    orth_ps = summary(temp_mod)$tTable[,4] # Orth p-values
    raw_coefs = as.vector(gamma %*% coef(temp_mod)) # Backtransforming to raw poly
    
    results_matrix[i,] = c(orth_ests, orth_ps, raw_coefs)
    if(i %% 100 == 0){print (i)}
  }
  
  ## Check that this is equal
  stopifnot(sum(is.na(results_matrix[,1])) == skipped)
  
  ## retrieve ASV identities
  rownames(results_matrix) = colnames(my_abundances)
  
  ## Calculating vertex of parabola
  # Formula for vertex: x = -x1/2*x2 (use raw estimates)
  vertex = -results_matrix[,8]/(2*results_matrix[,9])
  
  ## Categorising responses
  
  # For signs of curve and linear relationship, use orthogonal estimates
  none_pref = results_matrix[,5] > 0.05 &  # Neither model preferred
    results_matrix[,6] > 0.05 &
    !is.na(results_matrix[,1])
  
  # Linear curves
  lin_increase = results_matrix[,5] <= 0.05 & # Linear preferred, positive
    results_matrix[,6] > 0.05 &
    results_matrix[,2] >= 0 &
    !is.na(results_matrix[,1])
  lin_decrease = results_matrix[,5] <= 0.05 & # Linear preferred, negative
    results_matrix[,6] > 0.05 &
    results_matrix[,2] < 0 &
    !is.na(results_matrix[,1])
  
  # Quadratic upward curves
  quad_smiley = results_matrix[,6] <= 0.05 & # Vertex in middle
    results_matrix[,3] > 0 &
    vertex >= 0.5 &
    vertex <= 0.6 &
    !is.na(results_matrix[,1])
  quad_smiley_increase = results_matrix[,6] <= 0.05 & # Vertex to left
    results_matrix[,3] > 0 &
    vertex < 0.5 &
    !is.na(results_matrix[,1])
  quad_smiley_decrease = results_matrix[,6] <= 0.05 & # Vertex to right
    results_matrix[,3] > 0 &
    vertex > 0.6 &
    !is.na(results_matrix[,1])
  
  # Quadratic downward curves
  quad_sad = results_matrix[,6] <= 0.05 & # Vertex in middle
    results_matrix[,3] < 0 &
    vertex >= 0.5 &
    vertex <= 0.6 &
    !is.na(results_matrix[,1])
  quad_sad_decrease = results_matrix[,6] <= 0.05 & # Vertex to left
    results_matrix[,3] < 0 &
    vertex < 0.5 &
    !is.na(results_matrix[,1])
  quad_sad_increase = results_matrix[,6] <= 0.05 & # Vertex to right
    results_matrix[,3] < 0 &
    vertex > 0.6 &
    !is.na(results_matrix[,1])
  
  stopifnot(sum(none_pref, lin_increase, lin_decrease,
                quad_smiley, quad_smiley_increase, quad_smiley_decrease,
                quad_sad, quad_sad_increase, quad_sad_decrease, na.rm = TRUE) +
              sum(is.na(results_matrix[,1])) == nrow(results_matrix))
  
  # Saving out responses
  asv_response = factor(character(nrow(results_matrix)),
                        levels = c("insufficient data",
                                   "no change",
                                   "slow increase with SWD",
                                   "increase with SWD",
                                   "fast increase with SWD",
                                   "slow decrease with SWD",
                                   "decrease with SWD",
                                   "fast decrease with SWD",
                                   "moderate SWD favoured",
                                   "extreme SWD favoured"))
  
  #  asv_response = character(nrow(results_matrix))
  asv_response[none_pref] = "no change"
  asv_response[lin_increase] = "increase with SWD"
  asv_response[lin_decrease] = "decrease with SWD"
  asv_response[quad_smiley] = "extreme SWD favoured"
  asv_response[quad_smiley_increase] = "slow increase with SWD"
  asv_response[quad_smiley_decrease] = "fast decrease with SWD"
  asv_response[quad_sad] = "moderate SWD favoured"
  asv_response[quad_sad_increase] = "fast increase with SWD"
  asv_response[quad_sad_decrease] = "slow decrease with SWD"
  asv_response[is.na(asv_response)] = "insufficient data"
  
  # Responses df
  categorised_response = data.frame(asv = rownames(results_matrix),
                                    asv_response)
  categorised_response = cbind(categorised_response, results_matrix)
  
  return(categorised_response)
}

## Getting responses within each trajectory at t5
t5_bac_rewetup = get_asv_responses(clean_bac[89:99,], main_df$swd[89:99])
t5_bac_drydown = get_asv_responses(clean_bac[100:110,], main_df$swd[100:110])

t5_fun_rewetup = get_asv_responses(clean_fun[89:99,], main_df$swd[89:99])
t5_fun_drydown = get_asv_responses(clean_fun[100:110,], main_df$swd[100:110])

## Creating a combined dataframe
bac_combined_res = data.frame(asv = t5_bac_drydown$asv,
                              drydown = t5_bac_drydown$asv_response,
                              rewetup = t5_bac_rewetup$asv_response)
fun_combined_res = data.frame(asv = t5_fun_drydown$asv,
                              drydown = t5_fun_drydown$asv_response,
                              rewetup = t5_fun_rewetup$asv_response)

# Setting colours
response_cols = c("insufficient data" = "white",
                  "no change" = "#f2eeed",
                  
                  # More abundant at low SWD
                  "slow increase with SWD" = "lightblue",
                  "increase with SWD" = "blue",
                  "fast increase with SWD" = "darkblue",
                  
                  # More abundant at high SWD
                  "slow decrease with SWD" = "#f7613b",
                  "decrease with SWD" = "#fa3200",
                  "fast decrease with SWD" = "#821a00",
                  
                  # More abundant at middle or ends
                  "moderate SWD favoured" = "orange",
                  "extreme SWD favoured" = "purple")

## Reading in trees ====
# Trees were generated in FastTree from the top 500 most abundant prokaryotic and fungal ASVs
# For fungi, family-level information was used to graft taxa on to a backbone tree,
# and thus the tree only contains ASVs with available family-level taxonomy (N=236). 
bac_tree = read.tree("./data/Dataset_S4.txt")
fun_tree = read.tree("./data/Dataset_S5.txt")

## Adding features to tree ====
bactree_df = as_tibble(bac_tree)

# Responses
bactree_df = left_join(bactree_df, bac_combined_res,
                       by = join_by("label" == "asv"))
# Taxonomy
bactree_df = left_join(bactree_df, rownames_to_column(clean_bac_taxonomy),
                       by = join_by("label" == "rowname"))
bactree_df$dummy = 1 # For circles
bac_tree2 = as.treedata(bactree_df)

## Repeating for fungal tree
funtree_df = as_tibble(fun_tree)

# Responses
funtree_df = left_join(funtree_df, fun_combined_res,
                       by = join_by("label" == "asv"))
# Taxonomy
funtree_df = left_join(funtree_df, rownames_to_column(clean_fun_taxonomy),
                       by = join_by("label" == "rowname"))
funtree_df$dummy = 1 # For circles
fun_tree2 = as.treedata(funtree_df)

## Fig. 3: ASV responses on radial tree ====

## Setting colours and getting legends
bacphycols = c("#6eb148", "#7e92be", "#531E61",
               "#64a193", "#E85A78", "#4d6835",
               "#E03BB6", "#f1a7a0", "#bf913a",
               "#d15038", "#822b39", "#413657",
               "#25E03B", "#cd5c91", "#50392a")

funphycols = c("#d25f5d", "#531E61", "#25E03B",
               "#4C9C9B", "#E03BB6", "#bf913a",
               "#629638", "#822b39")

bactree_legend = get_legend(ggtree(bac_tree2)+
                              geom_fruit(geom = geom_col, 
                                         aes(x = dummy, y = asv, fill = Phylum),
                                         alpha = 0.5)+
                              scale_fill_manual(values = bacphycols))

funtree_legend = get_legend(ggtree(fun_tree2)+
                              geom_fruit(geom = geom_col, 
                                         aes(x = dummy, y = asv, fill = Phylum),
                                         alpha = 0.5)+
                              scale_fill_manual(values = funphycols))

# Responses
my_responses = c("insufficient data",
                 "no change",
                 "slow increase with SWD",
                 "increase with SWD",
                 "fast increase with SWD",
                 "slow decrease with SWD",
                 "decrease with SWD",
                 "fast decrease with SWD",
                 "moderate SWD favoured",
                 "extreme SWD favoured")

for_res_legend = ggtree(bac_tree2)+
  geom_fruit(geom = geom_col,
             aes(y = asv, x = dummy,
                 fill = drydown),
             pwidth = 0.05,
             show.legend = TRUE)+
  scale_fill_manual(name = "responses (outer rings)", 
                    labels = my_responses,
                    values = response_cols,
                    drop = FALSE)
responses_legend = get_legend(for_res_legend)

## Prokaryotes tree
p = ggtree(bac_tree2, layout = "radial", 
           branch.length = "none", 
           linewidth = 0.25)+
  geom_fruit(geom = geom_col, 
             aes(x = dummy, y = asv, fill = Phylum),
             pwidth = 0.05, width = 1,
             offset = -0.05, alpha = 0.5,
             show.legend = FALSE)+
  scale_fill_manual(values = bacphycols)

## Put lines on top of phyla colour
p$layers = c(p$layers[2], p$layers[1])

# Rotate bacteria tree for easier annotation
p = rotate_tree(p, 180)

fig3a = p +
  new_scale_fill()+
  # Drydown responses
  geom_fruit(geom = geom_col,
             aes(y = asv, x = dummy,
                 fill = drydown),
             #width = 0.5,
             pwidth = 0.05,
             show.legend = FALSE)+
  # Rewet-up responses
  geom_fruit(geom = geom_col,
             aes(y = asv, x = dummy,
                 fill = rewetup),
             #width = 0.5,
             pwidth = 0.05,
             show.legend = FALSE)+
  scale_fill_manual(name = "responses", 
                    labels = my_responses,
                    values = response_cols)+
  labs(title = "(a) Prokaryotes")+
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 24,
                                  vjust = -0.5))
fig3a

## Fungi tree 
pp = ggtree(fun_tree2, layout = "radial", 
            branch.length = "none", linewidth = 0.25)+
  geom_fruit(geom = geom_col, 
             aes(x = dummy, y = asv, fill = Phylum),
             pwidth = 0.05, width = 1,
             offset = -0.05, alpha = 0.5,
             show.legend = FALSE)+
  scale_fill_manual(values = funphycols)

## Put lines on top of phyla colour
pp$layers = c(pp$layers[2], pp$layers[1])

fig3b = pp +
  new_scale_fill()+
  # Drydown responses
  geom_fruit(geom = geom_col,
             aes(y = asv, x = dummy,
                 fill = drydown),
             #width = 0.5,
             pwidth = 0.05,
             show.legend = FALSE)+
  # Rewet-up responses
  geom_fruit(geom = geom_col,
             aes(y = asv, x = dummy,
                 fill = rewetup),
             pwidth = 0.05,
             show.legend = FALSE)+
  scale_fill_manual(name = "responses", 
                    labels = my_responses,
                    values = response_cols)+
  labs(title = "(b) Fungi")+
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 24,
                                  vjust = -0.5))
fig3b

fig3 = grid.arrange(fig3a, responses_legend, fig3b, nrow = 1,
                    widths = c(7,2,7))
#ggsave("./figures/fig3.svg", fig3, height = 8, width = 14)

## To see phyla legends: 
plot(bactree_legend)
plot(funtree_legend)

# Additional SI plots (Figs S6-8) ====

## (a) Figure S6: Alluvial plots ====

## Setting up alluvial plot dataframes
bac_sankey_df = bac_combined_res[!(bac_combined_res$drydown == "insufficient data" &
                                     bac_combined_res$rewetup == "insufficient data"),]
bac_sankey_df = bac_sankey_df[!(bac_sankey_df$drydown %in% c("no change", 
                                                             "insufficient data") &
                                  bac_sankey_df$rewetup %in% c("no change", 
                                                               "insufficient data")),]
bac_sankey_df = bac_sankey_df %>%
  make_long(drydown, rewetup)

# Reordering so increase is above, decrease below
bac_sankey_df$node = factor(bac_sankey_df$node, levels = my_responses[c(8,7,6,2,3,4,5,1,9,10)])

## Prokaryotes
figS6a = ggplot(bac_sankey_df, aes(x = x, 
                                   next_x = next_x, 
                                   node = node, 
                                   next_node = next_node,
                                   fill = factor(node),
                                   label = node)) +
  geom_alluvial(flow.alpha = 0.5) +
  geom_alluvial_text(
    data = bac_sankey_df %>% filter(x == "drydown"),
    size = 3.5, position = position_nudge(x = -0.07),
    hjust = 1)+
  geom_alluvial_text(
    data = bac_sankey_df %>% filter(x == "rewetup"),
    size = 3.5, position = position_nudge(x = +0.07),
    hjust = 0)+
  scale_fill_manual(values = response_cols, drop=FALSE)+
  theme_sankey(base_size = 18)+
  theme(axis.title = element_blank(),
        plot.title = element_text(size = 22, hjust = 0.5),
        plot.caption = element_text(size = 20, hjust = 0.5))+
  guides(label = "none",
         fill = "none")+
  labs(caption = "Prokaryotes")

## For fungi
fun_sankey_df = fun_combined_res[!(fun_combined_res$drydown == "insufficient data" &
                                     fun_combined_res$rewetup == "insufficient data"),]
fun_sankey_df = fun_sankey_df[!(fun_sankey_df$drydown %in% c("no change", 
                                                             "insufficient data") &
                                  fun_sankey_df$rewetup %in% c("no change", 
                                                               "insufficient data")),]
fun_sankey_df = fun_sankey_df %>%
  make_long(drydown, rewetup)

# Reordering
fun_sankey_df$node = factor(fun_sankey_df$node, levels = my_responses[c(8,7,6,2,3,4,5,1,9,10)])

figS6b = ggplot(fun_sankey_df, aes(x = x, 
                                   next_x = next_x, 
                                   node = node, 
                                   next_node = next_node,
                                   fill = factor(node),
                                   label = node)) +
  geom_alluvial(flow.alpha = 0.5) +
  geom_alluvial_text(
    data = fun_sankey_df %>% filter(x == "drydown"),
    size = 3.5, position = position_nudge(x = -0.07),
    hjust = 1)+
  geom_alluvial_text(
    data = fun_sankey_df %>% filter(x == "rewetup"),
    size = 3.5, position = position_nudge(x = +0.07),
    hjust = 0)+
  scale_fill_manual(values = response_cols, drop=FALSE)+
  theme_sankey(base_size = 18)+
  theme(axis.title = element_blank(),
        plot.title = element_text(size = 22, hjust = 0.5),
        plot.caption = element_text(size = 20, hjust = 0.5))+
  guides(label = "none",
         fill = "none")+
  labs(caption = "Fungi")

figS6a # Prokaryotes
figS6b # Fungi

# splitting up responses for bottom legend
legend_theme = theme(legend.key.spacing.y = unit(0.2, "cm"),
                     legend.text = element_text(size = 12))

res_plot1 = ggplot(bac_combined_res[bac_combined_res$drydown %in% c("no change",
                                                                    "insufficient data"),],
                   aes(x = drydown, fill = drydown))+
  geom_bar()+
  scale_fill_manual(name = "", values = response_cols)+
  legend_theme

res_plot2 = ggplot(bac_combined_res[bac_combined_res$drydown %in% c("slow increase with SWD",
                                                                    "increase with SWD",
                                                                    "fast increase with SWD"),],
                   aes(x = drydown, fill = drydown))+
  geom_bar()+
  scale_fill_manual(name = "", values = response_cols)+
  legend_theme

res_plot3 = ggplot(bac_combined_res[bac_combined_res$drydown %in% c("slow decrease with SWD",
                                                                    "decrease with SWD",
                                                                    "fast decrease with SWD"),],
                   aes(x = drydown, fill = drydown))+
  geom_bar()+
  scale_fill_manual(name = "", values = response_cols)+
  legend_theme

res_plot4 = ggplot(bac_combined_res[bac_combined_res$drydown %in% c("extreme SWD favoured",
                                                                    "moderate SWD favoured"),],
                   aes(x = drydown, fill = drydown))+
  geom_bar()+
  scale_fill_manual(name = "", values = response_cols)+
  legend_theme

res_legend1 = get_legend(res_plot1)
res_legend2 = get_legend(res_plot2)
res_legend3 = get_legend(res_plot3)
res_legend4 = get_legend(res_plot4)

S6_legend = grid.arrange(res_legend1, res_legend2, 
                          res_legend3, res_legend4, nrow = 1)

my_layout = rbind(c(1,1,1,1,2,2,2,2),
                  c(NA,3,3,3,3,3,3,NA))
figS6 = grid.arrange(figS6a, figS6b, nrow = 1)

figS6 = grid.arrange(figS6a,
                     figS6b,
                     S6_legend,
                     heights = c(1,0.2),
                     layout_matrix = my_layout)

#ggsave("figures/figS6.svg", figS6, width=14, height=8)

## (b) Figs. S7/S8: Phylogenetic signal ====

get_pairwise_bls = function(my_pairwise, my_responses, my_column = "drydown"){
  
  res_column = which(names(my_responses) == my_column)
  
  ## Get responses
  res_1 = my_responses$asv[my_responses[,res_column] == "no change"]
  res_2 = my_responses$asv[my_responses[,res_column] == "slow increase with SWD"]
  res_3 = my_responses$asv[my_responses[,res_column] == "increase with SWD"]
  res_4 = my_responses$asv[my_responses[,res_column] == "fast increase with SWD"]
  res_5 = my_responses$asv[my_responses[,res_column] == "slow decrease with SWD"]
  res_6 = my_responses$asv[my_responses[,res_column] == "decrease with SWD"]
  res_7 = my_responses$asv[my_responses[,res_column] == "fast decrease with SWD"]
  res_8 = my_responses$asv[my_responses[,res_column] == "moderate SWD favoured"]
  res_9 = my_responses$asv[my_responses[,res_column] == "extreme SWD favoured"]
  
  # Grouping responses as +ve and -ve
  res_10 = my_responses$asv[my_responses[,res_column] %in% c("slow increase with SWD",
                                                               "increase with SWD",
                                                               "fast increase with SWD")]
  res_11 = my_responses$asv[my_responses[,res_column] %in% c("slow decrease with SWD",
                                                               "decrease with SWD",
                                                               "fast decrease with SWD")]
  ## Saving to a list
  asv_list = list(res_1,res_2,res_3,res_4,res_5,
                  res_6,res_7,res_8,res_9,res_10,res_11)
  my_bls = numeric(length(asv_list))
  for (i in 1:length(asv_list)){
    curr_asvs = asv_list[[i]]
    # Get matrix corresponding to ASVs sharing the same response
    curr_bls = my_pairwise[which(rownames(my_pairwise) %in% curr_asvs), 
                           which(colnames(my_pairwise) %in% curr_asvs)]
    
    # Take average of upper triangle of matrix to get average distances between ASVs
    # sharing the same response
    my_bls[i] = mean(curr_bls[upper.tri(curr_bls)])
    
  }
  
  ## Getting a global value for all responses
  mod_bls = my_pairwise[which(rownames(my_pairwise) %in% asv_list[[8]]),
                        which(rownames(my_pairwise) %in% asv_list[[8]])]
  ext_bls = my_pairwise[which(rownames(my_pairwise) %in% asv_list[[9]]),
                        which(rownames(my_pairwise) %in% asv_list[[9]])]
  pos_bls = my_pairwise[which(rownames(my_pairwise) %in% asv_list[[10]]),
                        which(rownames(my_pairwise) %in% asv_list[[10]])]
  neg_bls = my_pairwise[which(rownames(my_pairwise) %in% asv_list[[11]]),
                        which(rownames(my_pairwise) %in% asv_list[[11]])]
  global_bls = c(mod_bls[upper.tri(mod_bls)],
                 ext_bls[upper.tri(ext_bls)],
                 pos_bls[upper.tri(pos_bls)],
                 neg_bls[upper.tri(neg_bls)])
  global_bls = mean(global_bls)
  
  ## Saving out values
  data.frame(response = c("no change",
                          "slow increase with SWD",
                          "increase with SWD",
                          "fast increase with SWD",
                          "slow decrease with SWD",
                          "decrease with SWD",
                          "fast decrease with SWD",
                          "moderate SWD favoured",
                          "extreme SWD favoured",
                          "all increasing responses",
                          "all decreasing responses",
                          "all responses",
                          "NA"),
             meandist = c(my_bls, global_bls, NA))
  
}

## Getting pairwise distance matrix
obs_bacdist = cophenetic.phylo(bac_tree)
obs_fundist = cophenetic.phylo(fun_tree)

# Getting responses of taxa represented on tree
bactree500 = bac_combined_res[bac_combined_res$asv %in% bac_tree$tip.label,]
funtree500 = fun_combined_res[fun_combined_res$asv %in% fun_tree$tip.label,]

## Bacteria
bac_drydownBL = get_pairwise_bls(obs_bacdist, bactree500, "drydown")
bac_rewetupBL = get_pairwise_bls(obs_bacdist, bactree500, "rewetup")

## Fungi
fun_drydownBL = get_pairwise_bls(obs_fundist, funtree500, "drydown")
fun_rewetupBL = get_pairwise_bls(obs_fundist, funtree500, "rewetup")

## Randomising responses for null model
run_null_model = function(my_null_pairwise, my_null_responses, 
                          my_null_column = "drydown",
                          n_reps = 1000){
  
  null_bls = matrix(nrow = n_reps, ncol = 13)
  
  null_res_column = which(names(my_null_responses) == my_null_column)
  
  for (rr in 1:n_reps){
    my_null_responses[,null_res_column] = sample(my_null_responses[,null_res_column])
    curr_null = get_pairwise_bls(my_null_pairwise, my_null_responses, my_null_column)
    null_bls[rr,] = curr_null$meandist
    if(rr %% 100 == 0){print(rr)}
  }
  
  colnames(null_bls) = curr_null$response
  return(null_bls)
}

## Running null model

## To replicate results in main text:
# set.seed(10325)
# bac_drydown_null = run_null_model(obs_bacdist, bactree500, "drydown", n_reps = 10000)
# bac_rewetup_null = run_null_model(obs_bacdist, bactree500, "rewetup", n_reps = 10000)
# fun_drydown_null = run_null_model(obs_fundist, funtree500, "drydown", n_reps = 10000)
# fun_rewetup_null = run_null_model(obs_fundist, funtree500, "rewetup", n_reps = 10000)

# Quick run
bac_drydown_null = run_null_model(obs_bacdist, bactree500, "drydown", n_reps = 1000)
bac_rewetup_null = run_null_model(obs_bacdist, bactree500, "rewetup", n_reps = 1000)
fun_drydown_null = run_null_model(obs_fundist, funtree500, "drydown", n_reps = 1000)
fun_rewetup_null = run_null_model(obs_fundist, funtree500, "rewetup", n_reps = 1000)

bac_drydown_null_df = melt(as.data.frame(bac_drydown_null))
bac_rewetup_null_df = melt(as.data.frame(bac_rewetup_null))
fun_drydown_null_df = melt(as.data.frame(fun_drydown_null))
fun_rewetup_null_df = melt(as.data.frame(fun_rewetup_null))

bac_drydown_null_df$trajectory = "drydown"
bac_rewetup_null_df$trajectory = "rewetup"
fun_drydown_null_df$trajectory = "drydown"
fun_rewetup_null_df$trajectory = "rewetup"

bac_nulls_df = rbind(bac_drydown_null_df, bac_rewetup_null_df)
fun_nulls_df = rbind(fun_drydown_null_df, fun_rewetup_null_df)

## Calculating percent of simulations where mean BL is less
bac_drydown_p = sapply(1:13, function(x){
  sum(bac_drydown_null[,x]<bac_drydownBL$meandist[x])/nrow(bac_drydown_null)
})

bac_rewetup_p = sapply(1:13, function(x){
  sum(bac_rewetup_null[,x]<bac_rewetupBL$meandist[x])/nrow(bac_rewetup_null)
})

fun_drydown_p = sapply(1:13, function(x){
  sum(fun_drydown_null[,x]<fun_drydownBL$meandist[x])/nrow(fun_drydown_null)
})

fun_rewetup_p = sapply(1:13, function(x){
  sum(fun_rewetup_null[,x]<fun_rewetupBL$meandist[x])/nrow(fun_rewetup_null)
})

bac_drydownBL$trajectory = "drydown"
bac_drydownBL$pvalue = bac_drydown_p

bac_rewetupBL$trajectory = "rewetup"
bac_rewetupBL$pvalue = bac_rewetup_p

fun_drydownBL$trajectory = "drydown"
fun_drydownBL$pvalue = fun_drydown_p

fun_rewetupBL$trajectory = "rewetup"
fun_rewetupBL$pvalue = fun_rewetup_p

bac_obsBLs = rbind(bac_drydownBL, bac_rewetupBL)
fun_obsBLs = rbind(fun_drydownBL, fun_rewetupBL)

## Plotting
levels(bac_nulls_df$variable)
responses_to_plot = levels(bac_nulls_df$variable)[2:12]

figS7 = ggplot(data = bac_nulls_df[bac_nulls_df$variable %in% responses_to_plot,])+
  geom_violin(aes(x = variable, y = value))+
  geom_point(data = bac_obsBLs[bac_obsBLs$response %in% responses_to_plot,],
             aes(x = response, y = meandist), col = "red",
             alpha = 1)+
  geom_text(data = bac_obsBLs[bac_obsBLs$response %in% responses_to_plot,],
            aes(x = response, y = meandist, label = format(signif(pvalue, 3), digits = 3)),# format(round(pvalue,3), nsmall = 2)),
            size = 3, vjust = -1.5, hjust = -0.38)+
  facet_wrap(.~trajectory, nrow = 2)+
  labs(x = "", y = "Mean pairwise distance",
       title = "Prokaryotes")+
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1, size = 10),
        axis.title.y = element_text(size = 16),
        strip.text = element_text(size = 14),
        plot.title = element_text(size = 20, hjust = 0.5))

figS7
#ggsave("./figures/figS7.svg", figS7, height = 6, width = 11)

figS8 = ggplot(data = fun_nulls_df[fun_nulls_df$variable %in% responses_to_plot,])+
  geom_violin(aes(x = variable, y = value))+
  geom_point(data = fun_obsBLs[fun_obsBLs$response %in% responses_to_plot,],
             aes(x = response, y = meandist), col = "red",
             alpha = 1)+
  geom_text(data = fun_obsBLs[fun_obsBLs$response %in% responses_to_plot,],
            aes(x = response, y = meandist, label = format(pvalue, digits = 3)),
            size = 3, vjust = 1.5, hjust = -0.38)+
  facet_wrap(.~trajectory, nrow = 2)+
  labs(x = "", y = "Mean pairwise distance",
       title = "Fungi")+
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1, size = 10),
        axis.title.y = element_text(size = 16),
        strip.text = element_text(size = 14),
        plot.title = element_text(size = 20, hjust = 0.5))

figS8
#ggsave("./figures/figS8.svg", figS8, height = 6, width = 11)
