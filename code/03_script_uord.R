#' DESCRIPTION:
#' Script for Unconstrained Ordination

# in-class ----------------------------------------------------------------
pacman::p_load(tidyverse,
               GGally,
               vegan)

df_iris <- iris %>% 
  as_tibble() %>% 
  janitor::clean_names()

df_iris %>%
  ggpairs(
    progress = FALSE,
    columns = c("sepal_length",
                "sepal_width",
                "petal_length",
                "petal_width"),
    aes(
      color = species, 
      alpha = 0.5
    )
  ) +
  theme_bw()

df_petal <- df_iris %>% 
  select(starts_with("petal_"))

## prcomp() = PCA function
obj_pca <- prcomp(x = df_petal,
       center = TRUE,
       scale = TRUE)

summary(obj_pca)

df_pca <- df_iris %>% 
  bind_cols(obj_pca$x)


df_pca %>% 
  ggplot(aes(
    x = species,   
    y = PC1
  )) +
  geom_boxplot() +
  labs(x = "Species",
       y = "Petal shape (PC1)")

## NMDS ##
data("dune")

dune %>% 
  as_tibble() %>% 
  select(1:3) %>% 
  ggpairs() +
  theme_bw()

m_bray <- vegdist(dune,
                  method = "bray")

obj_nmds <- metaMDS(comm = m_bray,
        k = 2)

data("dune.env")
dune.env

df_nmds <- dune.env %>% 
  as_tibble() %>% 
  bind_cols(obj_nmds$points) %>%
  janitor::clean_names()

df_nmds %>% 
  ggplot(aes(x = mds1,
             y = mds2,
             colour = use)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95,
               linetype = 2) + 
  theme_bw()+
  labs(x = "NMDS1",
       y = "NMDS2",
       color = "Land Use Intensity")

#permanova
adonis2(m_bray ~ use,
        data = df_nmds)
# lab ---------------------------------------------------------------------

# ============================================================
# EXERCISE: PCA using the iris dataset
# ============================================================

# In this exercise, you will perform a Principal Component
# Analysis (PCA) using all morphological measurements in the
# iris dataset and visualize multivariate trait patterns
# among species.

# 1. Using all four morphological variables
#    (Sepal.Length, Sepal.Width, Petal.Length, Petal.Width),
#    perform a PCA.

df_flower <- df_iris %>% 
  select(-species)

flower_pca <- prcomp(x = df_flower,
                  center = TRUE,
                  scale = TRUE)

summary(flower_pca)

# 2. Visualize flower morphology in PC axes whose cumulative
#    contribution exceeds 90%; color points by species.

df_flower_pca <- df_iris %>% 
  bind_cols(flower_pca$x)

df_flower_pca %>% 
  ggplot(aes(
    x = PC1,   
    y = PC2,
    color = species)) +
  geom_point()

# 3. Which morphological traits contribute most strongly to
#    the first and second principal components? How?

#Sepal length, petal width, and petal length are most important for pc1
#Sepal width is by far the most important for PC2
#These results can be seen during the initial prcomp() function

# ============================================================
# EXERCISE: NMDS using the BCI dataset
# ============================================================

# In this exercise, you will perform a Non-metric Multidimensional
# Scaling (NMDS) using the BCI tree community dataset and explore
# patterns in species composition among sites.

data("BCI", "BCI.env")

# 1. Using the BCI dataset, calculate a dissimilarity matrix
#    (e.g., Bray-Curtis) and perform NMDS.

bci_bray <- vegdist(BCI,
                  method = "bray")

bci_nmds <- metaMDS(comm = bci_bray,
                    k = 2)

# 2. Visualize sites in NMDS space.
#    - How are sites positioned relative to each other?
#    - Color or shape points by environmental groups or site
#      characteristics of your choice.

bci_env_nmds <- BCI.env %>% 
  as_tibble() %>% 
  bind_cols(bci_nmds$points) %>%
  janitor::clean_names()

bci_env_nmds %>% 
  ggplot(aes(x = mds1,
             y = mds2,
             colour = habitat)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95,
               linetype = 2) + 
  theme_bw()+
  labs(x = "NMDS1",
       y = "NMDS2",
       color = "Habitat")

# 3. Perform PERMANOVA to examine if communities are grouped
#    by the environmental variable you selected.
adonis2(bci_bray ~ habitat,
        data = bci_env_nmds)
