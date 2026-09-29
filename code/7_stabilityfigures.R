## Creating 3D plots for shifts in community composition away
# from the original community

## Fungal community composition

## Getting x and y points
point_x = main_df$swd[1:110]
point_y = main_df$fun_bccontrol[1:110]

range(point_x)
range(point_y)

## Creating size of plane
x = seq(0.25, 0.9, length.out = 70)
y = seq(0.25, 0.8, length.out = 70)

## Getting GAM results to create valleys with
valley1_x = predicted_data$swd[predicted_data$treatment == "field" &
                                 predicted_data$timepoint == "70 days"]
valley1_y = predicted_data$fun_bccontrol_gam$fit[predicted_data$treatment == "field" &
                                                   predicted_data$timepoint == "70 days"]

valley2_x = predicted_data$swd[predicted_data$treatment == "drought" &
                                 predicted_data$timepoint == "70 days"]
valley2_y = predicted_data$fun_bccontrol_gam$fit[predicted_data$treatment == "drought" &
                                                   predicted_data$timepoint == "70 days"]

## Create a valley along GAM predictions with maximum depth of -1
z = outer(x, y, function(x, y) {
  
  # Using approx to interpolate points for paths
  path1 = approx(valley1_x, valley1_y, xout = x, rule = 2)$y
  path2 = approx(valley2_x, valley2_y, xout = x, rule = 2)$y
  
  # Making valleys slope, decrease z exponentially with distance from the path
  valley1 = exp(-(y - path1)^2/ 0.01)
  valley2 = exp(-(y - path2)^2/ 0.01)
  -2*(valley1+valley2)
#  -1*pmax(valley1, valley2) # Don't allow valleys to add together
})

## Getting the z points for the actual data points
point_z = mapply(function(px, py) {
  ix = which.min(abs(x - px))
  iy = which.min(abs(y - py))
  z[ix, iy]
}, point_x, point_y)

## Creating a perspective plot
p = persp(
  x, y, z,
  theta = -40,
  phi = 80,
  expand = 0.6,
  col = "grey",
  border = NA,
  scale = TRUE,
  axes = FALSE,
  box = TRUE,
  shade = 1)

## Transform the points to 3d relative to this perspective
points_3d = trans3d(
  point_x,
  point_y,
  point_z,
  pmat = p
)

#svg('./figures/fig2d.svg', width = 8, height = 8)
par(mar = c(5, 9, 2, 2), xpd = NA)

## Regenerating perspective plot
persp(
  x, y, z,
  theta = -40,
  phi = 80,
  expand = 0.6,
  col = "white",
  scale = TRUE,
  border = NA,
  axes = TRUE,
  box = TRUE,
  ticktype = "detailed",
  shade = 0.8,
  xlab = "\nSWD", 
  ylab = "\nSimilarity to original fungal community",
  zlab = "")

# We'll make earlier timepoints more transparent
bg = ifelse(main_df$treatment[1:110] == "field", 
            mycols["field"], mycols["drought"])
alpha = rep(c(0.1,0.2, 0.3, 0.5, 1), each = 22)

bg_transparent = mapply(
  function(col, a) adjustcolor(col, alpha.f = a), bg,alpha)

points(
  points_3d$x[1:110],
  points_3d$y[1:110],
  pch = 21,
  col = "black",
  bg = bg_transparent[1:110],
  lwd = 1,
  cex = 2
)
dev.off()

## Prokaryotes ====

## Getting x and y points
point_x = main_df$swd[1:110]
point_y = main_df$bac_bccontrol[1:110]

range(point_x)
range(point_y)

## Creating size of plane
x = seq(0.25, 0.9, length.out = 70)
y = seq(0.59, 0.78, length.out = 70)

## Getting GAM results to create valleys with
valley1_x = predicted_data$swd[predicted_data$treatment == "field" &
                                 predicted_data$timepoint == "70 days"]
valley1_y = predicted_data$bac_bccontrol_gam$fit[predicted_data$treatment == "field" &
                                                   predicted_data$timepoint == "70 days"]

valley2_x = predicted_data$swd[predicted_data$treatment == "drought" &
                                 predicted_data$timepoint == "70 days"]
valley2_y = predicted_data$bac_bccontrol_gam$fit[predicted_data$treatment == "drought" &
                                                   predicted_data$timepoint == "70 days"]

## Create a valley along GAM predictions with maximum depth of -1
z = outer(x, y, function(x, y) {
  
  # Using approx to interpolate points for paths
  path1 = approx(valley1_x, valley1_y, xout = x, rule = 2)$y
  path2 = approx(valley2_x, valley2_y, xout = x, rule = 2)$y
  
  # Making valleys slope, decrease z exponentially with distance from the path
  valley1 = exp(-(y - path1)^2/ 0.0005) # Making deeper valley for prokaryotes as y range is smaller
  valley2 = exp(-(y - path2)^2/ 0.0005)
  -1*(valley1+valley2)
  #  -1*pmax(valley1, valley2) # Don't allow valleys to add together
})

## Getting the z points for the actual data points
point_z = mapply(function(px, py) {
  ix = which.min(abs(x - px))
  iy = which.min(abs(y - py))
  z[ix, iy]
}, point_x, point_y)

## Creating a perspective plot
p = persp(
  x, y, z,
  theta = -40,
  phi = 80,
  expand = 0.6,
  col = "grey",
  border = NA,
  scale = TRUE,
  axes = FALSE,
  box = TRUE,
  shade = 1)

## Transform the points to 3d relative to this perspective
points_3d = trans3d(
  point_x,
  point_y,
  point_z,
  pmat = p
)

#svg('./figures/fig2c.svg', width = 8, height = 8)
par(mar = c(5, 9, 2, 2), xpd = NA)

## Regenerating perspective plot
persp(
  x, y, z,
  theta = -40,
  phi = 80,
  expand = 0.6,
  col = "white",
  scale = TRUE,
  border = NA,
  axes = TRUE,
  box = TRUE,
  ticktype = "detailed",
  shade = 0.6,
  xlab = "\nSWD", 
  ylab = "\nSimilarity to original prokaryotic community",
  zlab = "")

# We'll make earlier timepoints more transparent
bg = ifelse(main_df$treatment[1:110] == "field", 
            mycols["field"], mycols["drought"])
alpha = rep(c(0.1,0.2, 0.3, 0.5, 1), each = 22)

bg_transparent = mapply(
  function(col, a) adjustcolor(col, alpha.f = a), bg,alpha)

points(
  points_3d$x[1:110],
  points_3d$y[1:110],
  pch = 21,
  col = "black",
  bg = bg_transparent[1:110],
  lwd = 1,
  cex = 2)
dev.off()

