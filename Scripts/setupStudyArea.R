

library(sf)
library(ggplot2)


# create a matrix with the bounding box of the study area
box.mat <- matrix(nrow=5, ncol = 2, 
                  data = c(0, 0, 0, 500, 500, 500, 500, 0, 0, 0), 
                  byrow = TRUE)

# convert to sf polygon object

box.poly <- st_polygon(x = list(box.mat))

#density grid

density.grid <- st_make_grid(box.poly, cellsize = 100, n = 50)

model.cells<- st_make_grid(box.poly, cellsize = 100, square = FALSE)

# create grid of point densities

df <- data.frame(expand.grid("X" = seq(5, 495, by = 10), 
                 "Y" = seq(5, 495, by = 10)), 
                 "M" = NA)

# TODO need to scale this so that it is in 0, 1

df$Z <- dnorm(df$X,100,100)*dnorm(df$Y,400,100)
ggplot(df) +
  geom_point(aes(x=X, y = Y, color= Z))


df$M <- dnorm(df$X, mean = 250, sd = 50) + df$Z

ggplot(df) +
  geom_point(aes(x=X, y = Y, color= M))

# convert to sf

density.surface <- st_sf(density.grid, df)

model.surface <- st_sf((st_crop(model.cells, st_bbox(density.grid))))

# want to keep only model cells that are completely within the bbox

# plot

ggplot()+
  geom_sf(data = density.surface, aes(color = M, fill = M)) +
  geom_sf(data = model.surface, fill = NA, color = "white") 




