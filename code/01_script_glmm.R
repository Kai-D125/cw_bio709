#' DESCRIPTION:
#' Script for GLMMs

# in-class ----------------------------------------------------------------

pacman::p_load(tidyverse,
               lme4,
               glmmTMB)

#sample data
data(Owls)
df_owl_raw <- as_tibble(Owls)

(df_owl <- df_owl_raw %>%                 # Start with raw owl dataset and assign cleaned version
    janitor::clean_names() %>%              # Standardize column names (lowercase, underscores, etc.)
    mutate(across(.cols = where(is.factor), # Select all factor-type columns
                  .fns = str_to_lower))     # Convert factor levels to lowercase
)

#visualization

df_owl %>%                          # Use cleaned owl dataset
  ggplot(aes(x = food_treatment,    # Map food treatment to x-axis
             y = neg_per_chick)) +  # Map number of negatives per chick to y-axis
  geom_boxplot(outliers = FALSE) +  # Draw boxplots without plotting outliers
  geom_jitter(alpha = 0.25) +       # Add jittered raw data points with transparency
  theme_bw()                        # Apply a clean black-and-white theme

#ignoring nest effects
MASS::glm.nb(sibling_negotiation ~ food_treatment +offset(log(brood_size)),
             data = df_owl)

#accounts for nest difference
MASS::glm.nb(sibling_negotiation ~ food_treatment + nest + offset(log(brood_size)),
             data = df_owl)

#visual for group-level effects 
v_g9 <- unique(df_owl$nest)[1:9]
df_owl %>% 
  filter(nest %in% v_g9) %>% 
  ggplot(aes(x = food_treatment,
             y = neg_per_chick))+
  geom_jitter(alpha = 0.25,
              width = 0.1)+
  facet_wrap(facets =~ nest,
             ncol = 3,
             nrow = 3) +
  theme_bw()

#use glmmTMB to include random effects
#syntax = y ~ x + (1 | group)
m_ri <- glmmTMB(sibling_negotiation ~       # Response variable
           food_treatment +          # Fixed effect of food treatment
           (1 | nest) +              # Random intercept for each nest
           offset(log(brood_size)),  # Offset to model rate per chick
         data = df_owl,              
         family = nbinom2())

#random intecept values
head(coef(m_ri)$cond$nest)

# ------------------------------------------------------------
# 1. Global (population-level) intercept on the response scale
# ------------------------------------------------------------
# The model uses a log link (negative binomial),
# so the intercept is on the log scale.
# We exponentiate it to return to the original response scale.
g0 <- exp(fixef(m_ri)$cond[1])


# ------------------------------------------------------------
# 2. Select a subset of nests to visualize
# ------------------------------------------------------------
# Plotting all 27 nests would be visually overwhelming,
# so we randomly select 9 nests for this example.
set.seed(123)  # ensures reproducibility
v_g9_ran <- sample(unique(df_owl$nest),
                   size = 9)


# ------------------------------------------------------------
# 3. Extract nest-specific coefficients (random intercept model)
# ------------------------------------------------------------
# coef(m_ri) returns the sum of fixed + random effects for each nest.
# These values are still on the log scale.
df_g9 <- coef(m_ri)$cond$nest %>% 
  as_tibble(rownames = "nest") %>%      # convert to tibble and keep nest ID
  filter(nest %in% v_g9_ran) %>%        # keep only the selected nests
  rename(
    log_g = `(Intercept)`,              # nest-specific intercept (log scale)
    b = food_treatmentsatiated          # fixed slope for food treatment
  ) %>% 
  mutate(
    g = exp(log_g),                     # intercept on response scale
    s = exp(log_g + b)                  # predicted value under satiated treatment
  )


# ------------------------------------------------------------
# 4. Create the figure
# ------------------------------------------------------------
df_owl %>% 
  filter(nest %in% v_g9_ran) %>%        # plot only the selected nests
  ggplot(aes(x = food_treatment,
             y = neg_per_chick)) +
  # Raw data points (jittered to reduce overlap)
  geom_jitter(width = 0.1,
              alpha = 0.5) +
  # Dashed horizontal lines:
  # nest-specific intercepts (baseline differences among nests)
  geom_hline(data = df_g9,
             aes(yintercept = g),
             alpha = 0.5,
             linetype = "dashed") +
  # Solid line segments:
  # predicted change from unfed to satiated treatment
  # using a common (fixed) slope across nests
  geom_segment(data = df_g9,
               aes(y = g,
                   yend = s,
                   x = 1,
                   xend = 2),
               linewidth = 0.5,
               linetype = "solid") +
  # Solid blue horizontal line:
  # global (population-level) intercept
  geom_hline(yintercept = g0,
             alpha = 0.5,
             linewidth = 1,
             linetype = "solid",
             color = "steelblue") +
  # Facet by nest to emphasize group-level structure
  facet_wrap(facets =~ nest,
             nrow = 3,
             ncol = 3) +
  theme_bw()

# lab ---------------------------------------------------------------------

# ============================================================
# EXERCISE: GLMM Exercise - `sleep` Data Set
# ============================================================

# ------------------------------------------------------------
# sleep dataset (built-in R dataset)
#
# This dataset contains measurements of increased sleep time
# after administering two different drugs.
#
# Structure:
# - 20 observations total
# - 10 subjects (each measured twice)
# - Paired / repeated-measures design
#
# Variables:
#   extra : Increase in hours of sleep compared to baseline
#   group : Indicates which drug was administered (factor with two levels ("1", "2"))
#   ID    : factor identifying each subject; Each subject appears once in each group
# ------------------------------------------------------------

# Q1 – Visualization:
# Compare sleep increase ("extra") between the two drug groups.
data("sleep")
df_sleep <- sleep %>% 
  janitor::clean_names()

df_sleep %>% 
  mutate(group = as.numeric(group)) %>% 
  ggplot(aes(x = group,
             y = extra,
             color = id)) + 
  geom_point() +
    geom_line() +
  facet_wrap(facets =~ id)
  

# Goals:
# - Show individual-level responses
# - Highlight paired structure (same subject in both groups)
# - Use color to identify subjects
# - facet by individual id

# Q2 - Model Development:
#
# Goal:
#   Examine how drug administration affects sleep duration.
#
# Key considerations:
#   - Response variable (extra) is continuous
#   - Drug (group) represents the treatment of interest
#   - Subject ID represents repeated measurements on the same
#     individuals
glmm_sleep <- glmmTMB(extra ~ group +          
          (1 | id),  
        data = df_sleep,              
        family = gaussian())

glmm_sleep_nrand <- glmmTMB(extra ~ group,  
        data = df_sleep,              
        family = gaussian())

# ============================================================
# EXERCISE: GLMM Exercise - `grouseticks` Data Set
# ============================================================

library(lme4)
data("grouseticks")

# ------------------------------------------------------------
# grouseticks dataset (from lme4 package)
#
# This dataset contains counts of parasitic ticks
# collected from red grouse chicks across multiple years
# and locations in Scotland.
#
# Structure:
# - 403 observations
# - Repeated measurements across broods and years
# - Count data with hierarchical (nested) structure
#
# Variables:
#   TICKS : Number of ticks counted on a chick
#   YEAR  : Sampling year
#   HEIGHT: height above sea level (meters)
#   LOCATION : Sampling site (grouping variable)
#   INDEX : Observation-level identifier
#
# Key features of the dataset:
# - Response variable is count data
# - Observations are grouped by brood and year
# ------------------------------------------------------------

# Q1 – Visualization:
#
# Goal:
#   Examine the relationship between parasite load (ticks) at the brood level and height above sea level.
#
# Key considerations:
# - Calculate average tick counts for each brood
# - Plot mean ticks vs. height
# - Color points by sampling year

# Q2 – Model development:
#
# Goal:
#   Develop a model to examine the relationship between parasite load (ticks) at the brood level and height above sea level.
#
# Key considerations:
#   - Response variable (TICKS) is count
#   - HEIGHT represents the variable of interest
#   - BROOD represents a grouping factor of repeated measurements
#   - YEAR represents another grouping factor of repeated measurements
