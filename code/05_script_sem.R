#' DESCRIPTION:
#' Script for SEM

# in-class ----------------------------------------------------------------

pacman::p_load(tidyverse,
               GGally,
               vegan,
               lavaan,
               lavaanPlot)

# Specify the URL of the raw CSV file on GitHub
url <- "https://raw.githubusercontent.com/aterui/biostats/master/data_raw/data_foodweb.csv"

# Read the CSV file into a tibble
(df_fw <- read_csv(url))

#visualize data
df_fw %>% 
  select(-plot_id) %>% 
  ggpairs(progress = FALSE) + 
  theme_bw()

#write a path diagram
m1 <- '
mass_herbiv ~ mass_plant + cv_h_plant
mass_pred ~ mass_herbiv
'

fit1 <- sem(model = m1, 
    data = df_fw)

summary(fit1, standardize = TRUE)

lavaanPlot(model = fit1, coefs = TRUE, stand = TRUE)

#model comparison
m2 <- '
mass_herbiv ~ mass_plant + cv_h_plant
mass_pred ~ mass_herbiv +cv_h_plant
'
fit2 <- sem(model = m2,
    data = df_fw)

#model comparison w/ anova
anova(fit1, fit2)

#SEM vs path analysis
#latent variables only exist in a statistical model

url <- "https://raw.githubusercontent.com/aterui/biostats/master/data_raw/data_herbivory.csv"

(df_herbv <- read_csv(url))

#visualization

df_herbv %>% 
  ggpairs(progress = FALSE,
          columns =  c("soil_n",
                       "sla",
                       "cn_ratio",
                      "per_lignin")) +
  theme_bw()

m_sem <- '
palatability =~ sla + cn_ratio + per_lignin
palatability ~ soil_n
herbivory ~ palatability
'

fit_sem <- sem(model = m_sem,
               data = df_herbv)

summary(fit_sem, standardize = TRUE)

#latent variable must be normally distributed

# lab ---------------------------------------------------------------------

# ============================================================
# EXERCISE: Path Analysis and Covariance Visualization
# ============================================================

#install.packages("piecewiseSEM")
pacman::p_load(tidyverse,
               GGally,
               piecewiseSEM)

data("keeley")

(df_keeley <- keeley %>% 
    as_tibble())

# The "keeley" dataset contains fire-related vegetation data
# collected from shrublands in California.
#
# ------------------------------------------------------------
# Column descriptions:
# elev  : Elevation of the site
#abiotic : idk 
# aspect: Slope aspect (orientation)
# heat  : Heat load index (a function of slope and aspect)
# firesev: Fire severity
# age   : stand age
# hetero : environmental heterogenaety
# cover : Vegetation cover
# rich  : Plant species richness
# ------------------------------------------------------------
#
# In this exercise, you will explore relationships among variables
# using covariance and path analysis. You will replicate a published
# path model and propose an alternative.

# 1. For the variables depicted in Figure 22.1, draw a figure
#    showing the covariance between variables.

df_keeley %>% 
  select(-elev) %>% 
  ggpairs(progress = FALSE) + 
  theme_bw()

# 2. Following Figure 22.1, develop a path model using the
#    same variables and relationships. Examine if this model
#    captures the data structure using a Chi-Square test.

m <- '
age ~ distance
hetero ~ distance
abiotic ~ distance
firesev ~ age
cover ~ firesev
rich ~ abiotic + hetero + cover'


fit <- sem(model = m,
    data = df_keeley)


summary(fit, standardize = TRUE)

# 3. Develop an alternative path model that you consider more
#    appropriate based on theory or observed data patterns.

m2 <- '
age ~ distance
hetero ~ distance
abiotic ~ distance
firesev ~ age + abiotic
cover ~ firesev + age 
rich ~ abiotic + hetero + cover + age + firesev
'

fit2 <- sem(model = m2,
           data = df_keeley)


summary(fit2, standardize = TRUE)

# 4. Compare the performance of the published model (Figure 22.1)
#    and your alternative model.
#    - Consider fit indices, path coefficients, and interpretability.

anova(fit, fit2)
