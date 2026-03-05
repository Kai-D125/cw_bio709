#' DESCRIPTION:
#' Script for piecewise SEM

# in-class ----------------------------------------------------------------
#install.packages("janitor") had to reinstall bc the update

pacman::p_load(tidyverse,
               GGally,
               piecewiseSEM,
               glmmTMB)

data("keeley")

(df_keeley <- keeley %>% 
    as_tibble())

#piecewise assuming normality
m1 <- lm(abiotic ~ distance, data = df_keeley)
m2 <- lm(hetero ~ distance, data = df_keeley)
m3 <- lm(firesev ~ age, data = df_keeley)
m4 <- lm(cover ~ firesev, data = df_keeley)
m5 <- lm(rich ~ cover + abiotic + hetero, data = df_keeley)

sem_model <- psem(m1, m2, m3, m4, m5)

summary(sem_model)

#piecewise w/ negative binomial assumption

m1 <- lm(abiotic ~ distance, data = df_keeley)
m2 <- lm(hetero ~ distance, data = df_keeley)
m3 <- lm(firesev ~ age, data = df_keeley)
m4 <- lm(cover ~ firesev + hetero, data = df_keeley)
m5 <- MASS::glm.nb(rich ~ cover + abiotic + hetero + distance,
                   data = df_keeley)

# m4 now includes a direct effect of hetero on cover (added path)
# m5 now models richness as negative binomial (MASS::glm.nb) 
# and includes direct effect of distance on richness (added path)

sem_model <- psem(m1, m2, m3, m4, m5)

summary(sem_model)

##Fisher's C = comparing proposed models to all possible models
## -2 times the sum of the logarithm of P values of possible models
## Follows chi-squared distribution
## big fisher c = missing links

#Visualization
plot(sem_model)

#Random effect
data("shipley")

df_shipley <- shipley %>% 
  as_tibble() %>% 
  janitor::clean_names() %>% 
  drop_na(growth)

df_shipley %>% 
  group_by(site) %>% 
  summarize(n_tree = n_distinct(tree))

#site = random effect

#visualization
df_shipley %>% 
  ggpairs(
    columns = c("dd",
      "date",
      "growth",
      "live"),
    progress = FALSE
  ) +
  theme_bw()

# Model 1: date depends on dd, with random intercepts for site and tree
m1 <- glmmTMB(date ~ dd + (1 | site) + (1 | tree), 
              data = df_shipley,
              family = "gaussian")

# Model 2: growth depends on date, same random effects
m2 <- glmmTMB(growth ~ date + (1 | site) + (1 | tree), 
              data = df_shipley,
              family = "gaussian")

# Model 3: live (binary) depends on growth, logistic mixed model
m3 <- glmmTMB(live ~ growth + (1 | site) + (1 | tree), 
              data = df_shipley, 
              family = "binomial")

sem_glmm <- psem(m1, m2, m3)

summary(sem_glmm, .progressBar = FALSE)

#marginal vs conditional only appear for glmm
#marginal = variance explained by fixed effect
#conditional = variance explained by random effect

#chi sqared model different between piecewise and normal sem

#no latent variables for piecewise

# lab ---------------------------------------------------------------------

library(piecewiseSEM)
data("meadows")

df_meadows <- meadows %>% 
    as_tibble() %>% 
  janitor::clean_names()

df_meadows <- df_meadows %>% 
  mutate(mass = log(mass))

# =========================================
# EXERCISE: Piecewise SEM with Meadows Data
# =========================================
#
# ------------------------------------------------------------
# Dataset: meadows (from piecewiseSEM package)
# Variables:
#   grazed - 0 = ungrazed, 1 = grazed
#   mass   - plant biomass (g/m²)
#   elev   - plot elevation above sea level
#   rich   - plant species richness per m²
# ------------------------------------------------------------
#
# 1. Explore the dataset (structure, summary, plots).

df_meadows %>% 
  ggpairs(
    columns = c("grazed",
                "mass",
                "elev",
                "rich"),
    progress = FALSE
  ) +
  theme_bw()

summary(df_meadows)

# 2. Develop a conceptual model: decide which variables influence others.
#    - Consider direct and indirect effects.
#    - Think about grazing as a disturbance factor.

# Hypothesis: Mass is affected by elevation, grazing, and richness
# and richness is affected by elevation and grazing
# Grazing is a random varibale

# 3. Fit component models (e.g., lm) for each hypothesized relationship.

mead1 <- lm(mass ~ rich + grazed, 
                 data = df_meadows)

mead2 <- lm(rich ~ elev + grazed, 
                 data = df_meadows)

# 4. Combine models into a piecewise SEM using psem().

pmeadows <- psem(mead1, mead2)

# 5. Evaluate the SEM: path coefficients, significance, variance explained.

summary(pmeadows)
plot(pmeadows)

# 6. Optional: try alternative models if your model deviates from the expectation.
# man i tried but hoooo boy this data
# Deliverables:
# - Code for component models and combined SEM
# - Conceptual SEM diagram
# - Short reasoning about your SEM results
