#' DESCRIPTION:
#' Script for time-series

# in-class ----------------------------------------------------------------
pacman::p_load(tidyverse,
               forecast,
               lterdatasampler,
               daymetr,
               glarma)

url <- "https://raw.githubusercontent.com/aterui/biostats/master/data_raw/data_ts_anormaly.csv"
(df_ts <- read_csv(url))

g_base <- df_ts %>% 
  ggplot(aes(x = year,
             y = anormaly)) +
  geom_line() +
  geom_point() 

m_ts <- lm(anormaly ~ year,
   data = df_ts)

summary(m_ts)

g_base + 
  geom_abline(intercept = coef(m_ts)[1],
              slope = coef(m_ts)[2])

#random-walk

y <- NULL
y[1] <- 0
for (i in 1:99){
y[i + 1] <- y[i] + rnorm(1, mean = 0, sd = 1)
}
tibble( y = y,
        x = 1:length(y)) %>% 
  ggplot(aes(x = x,
             y = y)) + 
  geom_point() + 
  geom_line()

## Auto-regressive (AR) model

df_huron <- tibble(
  year = time(LakeHuron),                # Extracts the time component (years) from the LakeHuron ts object
  water_level = as.numeric(LakeHuron)    # Converts LakeHuron values to numeric (from ts class)
) %>% 
  arrange(year)                           # Ensures the data is ordered by year

# Plot Lake Huron time series with a linear trend
df_huron %>% 
  ggplot(aes(x = year, y = water_level)) +
  geom_point(alpha = 0.25) +       # Semi-transparent points
  geom_line(linetype = "dotted") + # Dotted line connecting points
  geom_smooth(method = "lm",       # Linear trend line
              color = "black",
              linewidth = 0.5) +
  theme_bw() +
  labs(x = "Year", y = "Water Level")

#AR models need data from oldest to newest
m_ar1 <- Arima(
  df_huron$water_level,
  order = c(1, 0, 0)    #order must have 3 vectors
)

## fitted values
# Add fitted values from AR(1) model to the dataset
df_huron_ar1 <- df_huron %>% 
  mutate(fit = fitted(m_ar1) %>%   # Extract fitted (in-sample predicted) values from the AR(1) model
           as.numeric())           # Convert to numeric if necessary

# Plot observed and fitted values
df_huron_ar1 %>% 
  ggplot() +
  geom_point(aes(x = year, 
                 y = water_level),
             alpha = 0.25) +        # Plot observed water levels
  geom_line(aes(x = year, 
                y = fit),           # Plot AR(1) fitted values
            color = "steelblue") +
  theme_bw()                        # Clean black-and-white theme

## Moving Average (MA)
# averages past error
m_ma1 <- Arima(
  df_huron$water_level,
  order = c(0, 0, 1)
)

##Auto Regressive Moving Average (ARMA) model
m_arma1 <- Arima(
  df_huron$water_level,
  order = c(1, 0, 1)
)

## ARIMA model
m_arima1 <- Arima(
  df_huron$water_level,
  order = c(1, 1, 0)
)

## Model Selection
auto.arima(
  y = df_huron$water_level,
  stepwise = FALSE,
  ic = "aic"
)

## ARIMAX = ARIMA w/ outside predictors
data("ntl_icecover") 
df_ice <- ntl_icecover %>% 
  as.tibble() %>% 
  filter(between(year, 1980, 2014),
         lakeid == "Lake Mendota") %>% 
  arrange(year)

# Download daily climate data from Daymet for Lake Mendota
list_mendota <- download_daymet(
  site = "Lake_Mendota",   # Arbitrary name you assign to this site
  lat = 43.1,              # Latitude of the lake
  lon = -89.4,             # Longitude of the lake
  start = 1980,            # Start year
  end = 2024,              # End year
  internal = TRUE          # Return the data as an R object rather than saving to disk
)

df_temp <- list_mendota$data %>% 
  as.tibble() %>% 
  janitor::clean_names() %>% 
  mutate(
    date = as.Date(
      paste(year, yday, sep = "-"),
      format = "%Y-%j"
    ),
    month = month(date)
  ) %>% 
  arrange(year, yday) %>% 
  group_by(year) %>% 
  summarize(temp_min = round(mean(tmin_deg_c), 2)) 
  
df_ice <- df_ice %>% 
  left_join(df_temp, by = "year")

obj_arima <- auto.arima(
  y = df_ice$ice_duration,
  xreg = df_ice$temp_min,
  stepwise = FALSE
)

df_ice %>% 
  ggplot(aes(x = temp_min,
             y = ice_duration)) +
  geom_point()

confint(obj_arima)

# lab ---------------------------------------------------------------------

# ============================================================
# EXERCISE: Bison Body Mass, Climate, and Time-Series Analysis
# ============================================================

library(lterdatasampler)

# The "knz_bison" dataset contains long-term monitoring data
# on bison captured at Konza Prairie Biological Station.
#
# ------------------------------------------------------------
# Key columns may include:
# rec_year      : Year of capture
# animal_sex    : Sex of the individual (e.g., female, male)
# animal_weight : Body mass of bison
# ------------------------------------------------------------
#
# In this exercise, you will explore long-term trends in bison
# body mass and evaluate how climate variability may influence
# weight dynamics over time.

# 1. Explore the structure of the knz_bison dataset.
#    - Inspect variable types and missing values.
#    - Reformat variables as needed for analysis.

# 2. Subset the data to include observations from 1994–2012.

# 3. Calculate the average body mass for female and male bison
#    for each year in the selected time period.

df_knz <- knz_bison %>% 
  as.tibble() %>% 
  filter(between(rec_year, 1994, 2012)) %>% 
  group_by(rec_year, animal_sex) %>% 
  summarize(w = round(mean(animal_weight), 2)) %>% 
  rename(year = rec_year)

# 4. Obtain climate data from the daymetr dataset.
#    - Identify relevant climate variables (e.g., temperature,
#      precipitation).
#    - Associate climate data with knz_bison by year.
#    - Coordinates: Lat 39.09300	Lon -96.57500

knz_climate <- download_daymet(
  site = "knz",   # Arbitrary name you assign to this site
  lat = 39.09300,              # Latitude of the lake
  lon = -96.57500,             # Longitude of the lake
  start = 1994,            # Start year
  end = 2012,              # End year
  internal = TRUE          # Return the data as an R object rather than saving to disk
)

df_climate <- knz_climate$data %>% 
  as.tibble() %>% 
  janitor::clean_names() %>% 
  mutate(
    date = as.Date(
      paste(year, yday, sep = "-"),
      format = "%Y-%j"
    ),
    month = month(date)
  ) %>% 
  arrange(year, yday) %>% 
  group_by(year) %>% 
  summarize(t_max = round(max(tmax_deg_c), 2))

df_knz <- df_knz %>% 
  left_join(df_climate, by = "year")

# 5. Perform a time-series analysis to examine whether selected
#    climate variables influence annual bison body mass.
#    - Consider temporal autocorrelation and lag effects.
#    - Model males and females separately

df_knz_m <- df_knz %>% 
  filter(animal_sex == "M")

df_knz_f <- df_knz %>% 
  filter(animal_sex == "F")

obj_arima_f <- auto.arima(
  y = df_knz_f$w,
  xreg = df_knz_f$t_max,
  stepwise = FALSE
)

confint(obj_arima_f)

obj_arima_m <- auto.arima(
  y = df_knz_m$w,
  xreg = df_knz_m$t_max,
  stepwise = FALSE
)

confint(obj_arima_m)

# 6. Using your fitted model, compare observed bison body mass
#    with predicted values for the period 2014–2020.
#    - Evaluate model performance and discuss sources of uncertainty.
