# NOTE:
# When instructed to "test XXX", you must report the outcome as comments
# that clearly summarize the relevant statistical results (e.g., effect size,
# direction, significance, and interpretation).
# Providing code alone without documenting and interpreting the results
# in comments will result in point deductions.

# loading zone
pacman::p_load(tidyverse,
               glmmTMB, 
               GGally,
               vegan,
               lavaan,
               lavaanPlot,
               glarma,
               forecast,
               lterdatasampler,
               daymetr)

# dataset 1 ---------------------------------------------------------------

link1 <- "https://raw.githubusercontent.com/aterui/biostats/master/data_raw/data_insect_emergence.rds"
df_emg <- readRDS(url(link1, "rb"))

# This dataset ('df_emg') contains daily measurements of aquatic insect emergence
# from two wetland sites over a full calendar year (Jan 1–Dec 31).

# Data structure:
# t           : Day of the year (integer), where 1 = January 1 and 365 = December 31
# site        : Site identifier (factor), with "s1" and "s2" representing the two wetlands
# emergence   : Emergence flux of aquatic insects (g/day)

# Q1. Visualize seasonal patterns in emergence flux at both sites
#     (e.g., plot emergence vs. day of year, with separate lines or colors for each site).
#     [1 point]

df_emg %>% 
  ggplot(aes(x = t,
             y = emergence,
             colour = site)) +
  geom_line() +
  theme_bw()

# Q2. Test whether emergence flux differs significantly between the two sites,
#     while appropriately accounting for seasonal variation
#     [4 points]

emg_glm <- glmmTMB(
  emergence ~  
    site +               
    (1 | t), 
  data = df_emg,                    
  family = nbinom2()
)

summary(emg_glm)

# dataset 2 ---------------------------------------------------------------

link2 <- "https://raw.githubusercontent.com/aterui/cw_bio709/master/data_fmt/data_lake_invert.rds"
df_inv <- readRDS(url(link2, "rb"))

# This dataset 'df_inv' contains 100 observations from 10 lakes.
# Within each lake, 10 plots were established, spaced ~500 m apart.
# At each plot, the following variables were measured:

# s          : Species richness of invertebrates associated with aquatic plants at each plot
# hb         : Standing biomass of invertebrates associated with aquatic plants at each plot
# prod       : Production rate of aquatic plants (macrophytes), measured as g/month
# substrate  : Median diameter of substrate materials (mm)
# cond       : Water electrical conductivity (µS/cm);
#              a proxy for ionized nutrient levels (higher values may indicate eutrophication)
# lake       : lake ID

# Researcher's hypothesis was that: 
# (a) conductivity influences the productivity of macrophyes.
# (b) macrophyte's production rate ('prod') dictates invertebrate biomass ('hb') through bottom-up effects
# (c) macrophyte's production rate ('prod') dictates invertebrate richness ('s') through bottom-up effects 

# Q1. Create a scatter plot of macrophyte production ('prod', y-axis)
#     versus water conductivity ('cond', x-axis), with points colored by lake identity.
#     [1 point]

df_inv %>% 
  ggplot(aes(x = cond,
             y = prod,
             colour = lake)) +
  geom_point() +
  theme_bw()

# Q2. Create a scatter plot of raw invertebrate biomass ('hb', y-axis)
#     versus macrophyte production ('prod', x-axis), with points colored by lake identity.
#     [1 point]

df_inv %>% 
  ggplot(aes(x = prod,
             y = hb,
             colour = lake)) +
  geom_point() +
  theme_bw()

# Q3. Create a scatter plot of "log-transformed" invertebrate biomass ('hb', y-axis)
#     versus macrophyte production ('prod', x-axis), with points colored by lake identity.
#     [1 point]

df_inv %>% 
  ggplot(aes(x = prod,
             y = log(hb),
             colour = lake)) +
  geom_point() +
  theme_bw()

# Q4. Test hypothesis (a) by modeling macrophyte production while
#     statistically controlling for potential confounding variables ('substrate', 'lake').
#     [3 points]

inv_glm <- glmmTMB(
  prod ~  
    cond +               
    (1 + substrate | lake), 
  data = df_inv,                    
  family = gaussian()
)

summary(inv_glm)

# Q5. Test hypotheses (a–c) simultaneously using a unified modeling framework.
#     Based on the resulting statistical tests, determine whether the overarching
#     hypothesis (a–c, combined) is supported or rejected.
#     - Use appropriate probability distributions.
#     - Use variable transformation if appropriate given the data.
#     [4 points]

m_inv <- '
  prod ~ cond
  hb + s ~ prod
'

(fit_inv <- sem(model = m_inv,
             data = df_inv))

summary(fit_inv, standardize = TRUE)

lavaanPlot(model = fit_inv, coefs = TRUE, stand = TRUE)

## Conductivity greatly influences productivity, which in turn influences invertebrate richness
# but not invertebrate biomass

# dataset 3 ---------------------------------------------------------------

link3 <- "https://raw.githubusercontent.com/aterui/cw_bio709/master/data_fmt/nutrient.rds"
nutrient <- readRDS(url(link3, "rb"))

print(trees)

# This dataset ('trees') contains measurements of 31 felled black cherry trees.
# The three variables represent tree diameter, height, and timber volume.
# Note: the variable 'Girth' is actually the diameter measured at 4 ft 6 in above ground.

# Data structure:
# Girth   : Numeric, tree diameter in inches (mislabelled as girth)
# Height  : Numeric, tree height in feet
# Volume  : Numeric, timber volume in cubic feet

# Q1. Visualize relationships among tree diameter ('Girth'), height ('Height'),
#     and timber volume ('Volume') (e.g., using scatterplot matrix or pairwise scatter plots).
#     [1 point]

trees %>%
  ggpairs(
    progress = FALSE, 
    columns = c("Girth",
                "Height",
                "Volume"),
    aes(
      alpha = 0.5 
    )
  ) +
  theme_bw()

# Q2. Perform an appropriate ordination or dimension reduction method to 
#     summarize these three variables into fewer composite axes.
#     Then, identify and retain axes that explain meaningful variation in the original variables
#     [3 points]

pca_trees <- prcomp(
  x = trees,    
  center = TRUE,
  scale = TRUE)

summary(pca_trees)

df_nut <- bind_cols(
  nutrient,
  as_tibble(pca_trees$x))

# Q3. If justified, test whether the retained axis (or axes) is significantly 
#     related to "nutrient"; 
#     skip regression if the ordination does not support meaningful interpretation.
#     [1 point]

df_nut %>% 
  ggplot(aes(x = ...1,
             y = PC1)) +
  geom_line()

m_nut <- lm(PC1 ~ ...1,
        data = df_nut)

summary(m_nut)

# dataset 4 ---------------------------------------------------------------

df_nile <- dplyr::tibble(
  year = time(Nile), # observation year
  discharge = as.numeric(Nile) # discharge
)

df_sunspot <- dplyr::tibble(
  year = time(sunspot.year), # observation year
  sunspots = as.numeric(sunspot.year) # the number of sunspots
)

# These datasets contain:
# - df_nile    : Annual discharge of the Nile River (Nile dataset)
# - df_sunspot : Annual sunspot counts (sunspot.year dataset)

# Q1. Create a combined data frame aligning the observation years
#     (i.e., only include years present in both datasets)
#     [1 point]

df_sun <- df_sunspot %>% 
  filter(between(year, 1871, 1988))

df_river <- df_nile %>% 
  filter(between(year, 1871, 1988))

df_sun_river <- df_sun %>% 
  left_join(df_river, by = "year")

# Q2. Test whether the number of sunspots is significantly related to Nile's discharge
#     [4 points]

  (obj_arima <- auto.arima(
    df_sun_river$discharge, 
    xreg = df_sun_river$sunspots, 
    stepwise = FALSE 
  ))
confint(obj_arima, level = 0.95)
