###############################################################################
#
#    Imputation of plot margin heights the IDENT mortality data
#
#    Script: Roman Link  (roman.link@plant-ecology.de)	
#
###############################################################################

# This script imputes height values for the trees at the plot margins (which
# were not measured during the censi for practical limitations). Imputation is
# based on elevation measurements from drone flights (aggregated on the level of
# plot borders and edges in different directions) and scaled by the average size
# differences of the corresponding species in the plot.

# The resulting data are ONLY used for the computation of size-specific
# neighbourhood effects for the outer rows of the 5 x 5 grid of central trees as
# dropping these would reduce the number of observations by 64%. None of these
# (obviously imperfect) height estimates is used for any response tree.

# 1. Load packages ------------------------------------------------------------
# if pacman is not installed, install via install.packages("pacman").
# the package manager will take care of the installation of all subsequent 
# packages.
pacman::p_load(tidyverse, readxl)

# 2. Read in mortality dataset ------------------------------------------------
# read in raw data, convert height to m and replace trees with missing heights
# for 2018 with the corresponding values of 2017
mort <- read_xlsx("data/IDENT_mortality.xlsx", sheet = 2) %>% 
  mutate(height = coalesce(height_2018, height_2017) / 100) 

# calculate species-wise average plot heights
avgheights <- mort %>% 
  filter(position == "center") %>% 
  # get species-wise averages
  group_by(plot, species) %>% 
  mutate(height1 = mean(height, na.rm = TRUE)) %>% 
  # when all individuals of the species in the plot have died, take average of
  # the species in the same fertilization and species combination
  group_by(species, fert, plot_species) %>% 
  mutate(height2 = coalesce(height1, mean(height, na.rm = TRUE))) %>% 
  # summarize on species level taking the best available species level height
  group_by(plot, species) %>% 
  summarize(height = mean(height2))
avgheights

# 3. Compute remote-sensing based plot margin tree heights -------.............
rsheights <- read_xlsx("data/IDENT_edge_tree_data.xlsx", sheet = 2) 
rsheights

# merge remote-sensing based height differences with actual plot averages for each
# species
rsheights1 <- full_join(rsheights, avgheights) %>% 
  # calculate the edge and corner height as a difference to the average height
  mutate(height_edge   = as.character(edge_diff_mean + height), # treated as characters to avoid problems during reshaping
         height_corner = as.character(corner_diff_mean + height)) %>% 
  # select relevant columns
  select(plot, species, 
         pos_edge = edge, pos_corner = corner, 
         height_edge, height_corner) %>% 
  # create temporary row ID to break ties in reshaping
  mutate(ID = 1:nrow(.)) %>% 
  # pivot to long format
  pivot_longer(cols = c(-plot, - ID, -species), names_to = c("var", "pos"), 
               values_to = "value", names_sep = "_") %>% 
  # reshape to correct wide format
  spread(var, value) %>% 
  # convert RS heights to numeric
  mutate(rsheight = as.numeric(height)) %>%
  # replace a low number of below-zero heights (RS difference larger than average 
  # species height - mostly in plots with many species) by plot edge averages 
  group_by(plot) %>% 
  mutate(rsheight = ifelse(rsheight <= 0, 
                           mean(rsheight[rsheight>=0]), 
                           rsheight)
  ) %>% 
  # remove temporary columns
  select(-ID, -height) %>% 
  ungroup()
rsheights1

# 4. Merge RS heights with mortality dataset ----------------------------------
mort1 <- mort %>% 
  mutate(
    # create position indicator for merging based on planting grid
    pos = case_when(tree == "1.1" ~ "NW",
                    tree == "1.7" ~ "NE",
                    tree == "7.1" ~ "SW",
                    tree == "7.7" ~ "SE",
                    tree %in% paste0("1.", 2:6) ~ "N",
                    tree %in% paste0("7.", 2:6) ~ "S",
                    tree %in% paste0(2:6, ".1") ~ "W",
                    tree %in% paste0(2:6, ".7") ~ "E",
                    TRUE ~ NA_character_
    )
  ) %>% 
  # join with remote sensing heights
  left_join(rsheights1) %>% 
  mutate(# add observations for the edge trees
    height_new = coalesce(height, rsheight),
    # set dead trees to a height of 0 (use of coalesce avoids problems with missing values)
    height_new = ifelse((is.na(height_2018) & position == "center") | (coalesce(mortality_2017, "") == "dead" & position == "edge") , 0, height_new)) %>% 
  group_by(plot, position, species) %>% 
  ungroup()
mort1

# export updated dataset
write_csv(mort1, "output/IDENT_mortality_RS_heights.csv")
