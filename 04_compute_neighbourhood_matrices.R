###############################################################################
#
#    Neighbourhood interaction matrices for the IDENT mortality data
#
#    Script: Roman Link  (roman.link@plant-ecology.de)	
#
###############################################################################

# This script computes neighbourhood interaction matrices for each plot of the
# Freiburg IDENT trial. For each of the 5 X 5 central trees per plot, the Hegyi
# competition index is computed separately for each species based on its eight
# immediate neighbour trees. As we are dealing with young trees where variation
# in height is much larger than variation in diameter, the index is computed on
# a height basis.
# The Hegyi is defined as sum(height_neighbour / (height_central * distance)),
# where height_neighbour is a vector of neighbour tree heights, height_central
# is the height of the central tree and distance is the distance to the central
# tree.
# By summing these up separately for each species, it is possible to separate
# neighbourhood interactions between all possible pairs of central and neighbour
# trees.
# In addition, the script is used to count the number of beetle-infested trees
# among the eight immediate neighbours of each tree.

# 1. Load packages ------------------------------------------------------------
# if pacman is not installed, install via install.packages("pacman").
# the package manager will take care of the installation of all subsequent 
# packages.
pacman::p_load(tidyverse, readxl)

# 2. Read in raw data ---------------------------------------------------------
mort <- read_csv("output/IDENT_mortality_RS_heights.csv") 

# get species IDs
specs <- sort(unique(mort$species))

# 3. Inspect data from one plot -----------------------------------------------
# create subset for one plot
data <- filter(mort1, plot_cn == "I.1.100")
data
# visualize plot positions
ggplot(data, aes(x = x_pos, y = y_pos, col = height_new, label = tree)) + geom_text()
# within plots the trees are numbered according to a regular grid - the neighbour
# trees of each tree can be identified this way


# 4. Define sets of neighbour trees -------------------------------------------
dat <- select(data, tree, x_pos, y_pos, position)
dat

# function for first order neighbours
first_order <- function(data){
  filter(dat, x_pos == data$x_pos & y_pos %in% (data$y_pos + c(-1, 1)) |
              y_pos == data$y_pos & x_pos %in% (data$x_pos + c(-1, 1)))$tree
}

# function for second order neighbours
second_order <- function(data){
  filter(dat, x_pos %in% (data$x_pos + c(-1, 1)) & y_pos %in% (data$y_pos + c(-1, 1)))$tree
}

# get sets of neighbours of first order
neighbours <- dat %>%
  filter(position == "center") %>% 
  group_nest(tree) %>% 
  mutate(first_order  = map(data, first_order),
         second_order = map(data, second_order)) %>% 
  unnest(data)
neighbours

# prepare empty neighbourhood matrix
neighmat <- data.frame(
  tree = unique(mort$tree[mort$position == "center"]),
  matrix(0, 
         nrow = 25,
         ncol = length(specs)
  ),
  stringsAsFactors = FALSE
) %>% 
  set_names(c("tree", specs))
neighmat

# 5. Write neighbourhood function based on Hegyi ------------------------------
# function that takes a data subset of one plot and returns a full neighbourhood
# effect matrix
neighbour_fun <- function(data){
  # get lookup table for tree numbers
  trees <- data$tree
  # get lookup table for tree species identities
  species <- data$species
  # get lookup table for tree status
  alive <- !(data$mortality_2017 %in% "dead")
  # get lookup table for tree height
  heights <- data$height
  # prepare neighbourhood matrix  
  neigh <- neighmat
  
  # calculate Hegyi contribution for each neighbour species
  for (i in 1:nrow(neigh)){
    for (j in 2:ncol(neigh)) {
      # set values for dead trees to zero
      if(heights[trees == neigh$tree[i]] == 0) {
        neigh[i, j] <- NA
      } else{ 
        # if tree is alive, calculate contribution of the neighbours of each species to the Hegyi index
        neigh[i, j] <- (sum(heights[ alive & species == specs[j-1] & trees %in% neighbours$first_order[[i]] ]) + 
                          1/sqrt(2) * sum(heights[ alive & species == specs[j-1]  & trees %in% neighbours$second_order[[i]] ])) / heights[trees == neigh$tree[i]]
      }
    }
  }
  
  # merge with unique tree identifiers and return
  suppressMessages(
    select(data, tree_cn, tree) %>% 
      left_join(neigh)
  )
}

# test if it works
neighbour_fun(mutate(data, height = height_new)) %>%  as.data.frame
# NA for dead trees and plot margin trees, else relative contribution to the Hegyi

# 6. Compute and export  neighbourhood matrix for all central trees -----------
# compute matrix for all plots
full_neigh_mat <- mort1 %>% 
  mutate(height = height_new) %>% 
  group_nest(plot_cn) %>% 
  mutate(neighbours = map(data, neighbour_fun)) %>% 
  select(-data) %>% 
  unnest(neighbours)
full_neigh_mat

# save results
write_csv(full_neigh_mat, "output/IDENT_neighbourhood_matrix.csv")


# 7. Write function to compute the number of neighbour trees with beetles -----
neighbour_fun_beetles <- function(data){
  # get lookup table for tree numbers
  trees <- data$tree
  # get lookup table for tree status
  alive <- !(data$mortality_2017 %in% "dead")
  # get lookup table for tree height
  beetles <- data$beetles
  
  # prepare vector for output
  beet <- numeric(25)
  
  # calculate number of beetle cases among the 8 neigbor trees 
  for (i in 1:length(beet)){
    beet[i] <- sum(beetles[ alive & trees %in% c(neighbours$first_order[[i]], neighbours$second_order[[i]]) ])
  }
  
  # merge with unique tree identifiers and return
  suppressMessages(
    select(data, tree_cn, tree) %>% 
      left_join(tibble(tree = neighbours$tree, beetles_neigh = beet))
  )
}

# test if it works
neighbour_fun_beetles(data) %>%  as.data.frame


# 8. Compute and export neighbour beetle infestation for all central trees ----
beetles <- mort1 %>% 
  group_nest(plot_cn) %>% 
  mutate(neighbours = map(data, neighbour_fun_beetles)) %>% 
  select(-data) %>% 
  unnest(neighbours)
beetles

write.csv(beetles, "output/IDENT_neighbours_beetles.csv", row.names = FALSE)

