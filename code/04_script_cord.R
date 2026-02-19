#' DESCRIPTION:
#' Script for Constrained Ordination

# in-class ----------------------------------------------------------------
pacman::p_load(tidyverse,
               GGally,
               vegan)

data("varespec", "varechem")

#response matrix

m_y <- varespec
colnames(m_y) <- str_to_lower(colnames(m_y))

df_env <- as_tibble(varechem) %>% 
  janitor::clean_names()

#visualization
m_y %>% 
  ggpairs(
    progress = FALSE,
    columns = 1:3,
    aes(alpha = 0.25)
  )

#RDA
obj_rda <- rda(m_y ~ n + p + ca,
    data = df_env)

#stats test
anova.cca(obj_rda,
          by = "margin",
          permutations = 999)

#type 1 anova = tests one factor after another, default anova  
#type 2 anova = tests factors w/o caring about order, by = "margin"
#type 3 anova = ???

#visualization
df_rda <- scores(obj_rda,
       display = "site",
       scaling = 2) %>% 
  bind_cols(df_env) %>% 
  janitor::clean_names()

#env vectors

df_bp <- scores(obj_rda,
                display = "bp",
                scaling = 2) %>% 
  as_tibble(rownames = "variable") %>% 
  janitor::clean_names()

#sorry if things get a little out of order here
#I copied and pasted the wrong visualization for the rda
#it should all be there and working though

#test w anova.cca

#upon compiling I have found that this part wasn't working 
#due to the mix up detailed above

#anova.cca(obj_db, 
          #by = "margin",
          #permutations = 999)

# Extract site (sample) scores from the RDA object
# - display = "sites": returns the coordinates of sampling sites in RDA space
# - scaling = 2: correlation scaling, which emphasizes relationships between sites
#   and explanatory variables (i.e., how site positions correlate with predictors)
# The site scores are then combined with the environmental data for visualization
df_rda <- scores(obj_rda, 
                 display = "sites",
                 scaling = 2) %>% 
  bind_cols(df_env) %>%           # append environmental variables (e.g., soil chemistry)
  janitor::clean_names()          # standardize column names for tidy workflows

# RDA vectors for environmental predictors
df_bp <- scores(obj_rda, 
                display = "bp", 
                scaling = 2) %>% 
  as_tibble(rownames = "variable") %>% 
  janitor::clean_names()

# Create a ggplot2 ordination plot
# - Points represent sites positioned by their constrained community composition
# - Color gradient reflects the nitrogen (n) concentration at each site
df_rda %>% 
  ggplot(aes(x = rda1,
             y = rda2)) +        # color sites by nitrogen level
  geom_point(aes(color = n)) +
  geom_segment(data = df_bp,
               aes(x = 0, xend = rda1 * 10, # 10 is arbitrary scaling for visualization
                   y = 0, yend = rda2 * 10),
               arrow = arrow(length = unit(0.2, "cm"))
  ) +
  geom_text(data = df_bp,
            aes(x = rda1 * 10.5,    # slightly beyond arrow tip
                y = rda2 * 10.5,
                label = variable),  # or use a variable column
            size = 4) +
  theme_bw() +
  labs(x = "RDA1",
       y = "RDA2",
       color = "Nitrogen") +
  scale_color_viridis_c()

obj_db <- dbrda(m_y ~ n + p + ca,
                data = df_env,
                distance = "bray")

# Extract site scores from the dbRDA object
df_db <- scores(obj_db, 
                display = "sites",
                scaling = 2) %>% 
  as_tibble() %>%              
  bind_cols(df_env) %>%        
  janitor::clean_names()       

# dbRDA vectors for environmental predictors
df_bp <- scores(obj_db, 
                display = "bp", 
                scaling = 2) %>% 
  as_tibble(rownames = "variable") %>% 
  janitor::clean_names()

# Create a ggplot2 ordination plot
df_db %>% 
  ggplot(aes(x = db_rda1,
             y = db_rda2)) +        # color sites by nitrogen level
  geom_point(aes(color = n)) +
  geom_segment(data = df_bp,
               aes(x = 0, xend = db_rda1,
                   y = 0, yend = db_rda2),
               arrow = arrow(length = unit(0.2, "cm"))
  ) +
  geom_text(data = df_bp,
            aes(x = db_rda1 * 1.1,    # slightly beyond arrow tip
                y = db_rda2 * 1.1,
                label = variable),  # or use a variable column
            size = 4) +
  theme_bw() +
  labs(x = "dbRDA1",
       y = "dbRDA2",
       color = "Nitrogen") +
  scale_color_viridis_c()

# lab ---------------------------------------------------------------------

# ============================================================
# EXERCISE: Community Ordination and Environmental Gradients
# ============================================================

library(vegan)
data("mite", "mite.env")

# The mite datasets contain information on Oribatid mite communities
# sampled from a small peatland area (2.5 m × 10 m).
#
# There are linked datasets:
# ------------------------------------------------------------
# mite     : Species abundance data (35 mite species × 70 sites)
# mite.env : Environmental variables measured at the same sites
# ------------------------------------------------------------
#
# Environmental variable descriptions (mite.env):
# ------------------------------------------------------------
# SubsDens : Substrate density (g/L)
# WatrCont : Water content of the substrate (g/L)
# Substrate: Substrate type (factor with multiple levels)
# Shrub    : Shrub density (ordered factor: low → high)
# Topo     : Microtopography (Blanket vs Hummock)
# ------------------------------------------------------------

# 1. Explore and visualize interrelationships among species abundances.
#    - Examine patterns of co-occurrence.
#    - Assess whether relationships among species appear linear or nonlinear.

colnames(mite) <- str_to_lower(colnames(mite))
colnames(mite.env) <- str_to_lower(colnames(mite.env))

m_mite_trans <- vegan::wisconsin(mite)

ggpairs(m_mite_trans,
        columns = 1:5,
        progress = FALSE)

# 2. Fit a redundancy analysis (RDA) model using environmental variables of your choice.
#    - Visualize the ordination results.
#    - Examine gradients and species–environment relationships.
#    - Evaluate whether the assumptions of RDA are appropriate for these data.

obj_mite <- rda(m_mite_trans ~ subsdens + shrub + substrate,
               data = mite.env)

anova.cca(obj_mite, 
          by = "margin", 
          permutations = 999)

df_mite <- scores(obj_mite,
                 display = "sites",
                 scaling = 2) %>% 
  bind_cols(mite.env) %>% 
  janitor::clean_names() %>% 
  ggplot(aes(x = rda1,
             y = rda2)) +
  geom_point()

# 3. Apply alternative ordination methods.
#    - Canonical correspondence analysis (CCA; see ?cca()).
#    - Distance-based RDA (dbRDA).

mite_cca <- cca(m_mite_trans ~ subsdens + shrub + substrate,
                data = mite.env)

mite_db <- dbrda(m_mite_trans ~ subsdens + shrub + substrate,
                 data = mite.env,
                 distance = "bray")

# 4. Compare RDA, CCA, and dbRDA.
#    - Perform permutation analysis to examine the significance of predictor variables
#    - Discuss which method is most appropriate for these data and why.

#dbRDA
df_db <- scores(mite_db, 
                display = "sites",
                scaling = 2) %>% 
  as_tibble() %>%              
  bind_cols(mite.env) %>%        
  janitor::clean_names()       

df_bp <- scores(mite_db, 
                display = "bp", 
                scaling = 2) %>% 
  as_tibble(rownames = "variable") %>% 
  janitor::clean_names()

# Create a ggplot2 ordination plot
df_db %>% 
  ggplot(aes(x = db_rda1,
             y = db_rda2)) +        # color sites by nitrogen level
  geom_point(aes(color = subsdens)) +
  geom_segment(data = df_bp,
               aes(x = 0, xend = db_rda1,
                   y = 0, yend = db_rda2),
               arrow = arrow(length = unit(0.2, "cm"))
  ) +
  geom_text(data = df_bp,
            aes(x = db_rda1 * 1.1,    # slightly beyond arrow tip
                y = db_rda2 * 1.1,
                label = variable),  # or use a variable column
            size = 4) +
  theme_bw() +
  labs(x = "dbRDA1",
       y = "dbRDA2",
       color = "Substrate Density") +
  scale_color_viridis_c()

#cca

df_db <- scores(mite_cca, 
                display = "sites",
                scaling = 2) %>% 
  as_tibble() %>%              
  bind_cols(mite.env) %>%        
  janitor::clean_names()       

df_bp <- scores(mite_cca, 
                display = "bp", 
                scaling = 2) %>% 
  as_tibble(rownames = "variable") %>% 
  janitor::clean_names()

#I can't seem to get the following code to work despite it being copy/pasted from above

# Create a ggplot2 ordination plot
# df_db %>% 
#ggplot(aes(x = db_rda1,
           #y = db_rda2)) +        # color sites by nitrogen level
  #geom_point(aes(color = subsdens)) +
  #geom_segment(data = df_bp,
               #aes(x = 0, xend = db_rda1,
                   #y = 0, yend = db_rda2),
               #arrow = arrow(length = unit(0.2, "cm"))
  #) +
  #geom_text(data = df_bp,
            #aes(x = db_rda1 * 1.1,    # slightly beyond arrow tip
                #y = db_rda2 * 1.1,
                #label = variable),  # or use a variable column
            #size = 4) +
  #theme_bw() +
  #labs(x = "dbRDA1",
       #y = "dbRDA2",
       #color = "Substrate Density") +
  #scale_color_viridis_c()

#CCA and/or dbRDA are best becasue the data has a horseshoe in it