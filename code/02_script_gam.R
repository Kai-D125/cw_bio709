#' DESCRIPTION:
#' Script for GAMs

# in-class ----------------------------------------------------------------

pacman::p_load(tidyverse,
               ggeffects,
               mgcv)
link <- "https://raw.githubusercontent.com/aterui/biostats/master/data_raw/data_water_temp.csv"

(df_wt_raw <- read_csv(link))

## check data format
sapply(df_wt_raw, class)

# Start from the raw water-temperature dataset
(df_wt <- df_wt_raw %>% 
    mutate(
      # Convert the character datetime column into a Date object.
      # The format argument specifies month/day/year.
      date = as.Date(date_time,
                     format = "%m/%d/%Y"),
      # Extract the calendar year from the Date object
      year = year(date),
      # Extract the month (1–12) from the Date object
      month = month(date)
    ) %>% 
    # Subset the data to observations from the year 2022 only
    # March through October
    filter(year == 2022,
           between(month, 3, 10)))

#summarize by date and site
# Aggregate raw water temperature data to daily averages
df_wt_daily <- df_wt %>% 
  # Group the data by date and site
  # This ensures that averaging is done separately for each site on each day
  group_by(date, site) %>% 
  # Compute summary statistics for each group
  summarize(
    # Take the mean of the 'temp' column within each group
    # na.rm = TRUE ignores missing values (NA) in the calculation
    # round(..., 3) rounds the resulting mean to 3 decimal places
    temp = mean(temp, na.rm = TRUE) %>% round(3)
  )

# Visualize daily water temperature for each wetland site
df_wt_daily %>% 
  ggplot(aes(
    x = date,      # Date on the x-axis
    y = temp,      # Daily averaged temperature on the y-axis
    color = site   # Color points by wetland type (woody vs open)
  )) +
  # Add points for each observation
  # alpha = 0.25 makes points semi-transparent to reduce overplotting
  geom_point(alpha = 0.25) +
  # Apply a clean black-and-white theme for readability
  theme_bw() +
  # Customize axis labels and legend title
  labs(
    x = "Date",                # x-axis label
    y = "Water Temperature",   # y-axis label
    color = "Wetland Type"     # legend title for the color mapping
  )

# Add new variables and ensure proper data types
df_wt_daily <- df_wt_daily %>% 
  mutate(
    # Convert the date column to Julian day (day of year)
    # Useful for modeling seasonal trends
    j_date = yday(date),
    # Ensure 'site' is treated as a factor (categorical variable)
    site = factor(site)
  )

# Fit a Generalized Linear Model (GLM) with Gaussian family
m_glm <- glm(
  temp ~ j_date + site,  # Model daily temperature as a function of Julian day and wetland site
  data = df_wt_daily,    # Dataset to use
  family = "gaussian"    # Specify Gaussian distribution (standard linear regression)
)

summary(m_glm)

# Generate model predictions across all Julian days and wetland sites
df_pred <- ggpredict(m_glm,
                     terms = c(
                       "j_date [all]",  # Use all observed values of Julian day
                       "site [all]"     # Generate predictions for all levels of the factor 'site'
                     )
) %>% 
  # Rename the default columns to match the original dataset
  rename(site = group,  # 'group' from ggpredict() corresponds to the factor variable 'site'
         j_date = x     # 'x' from ggpredict() corresponds to the predictor 'j_date'
  )

# Plot daily water temperature and overlay model predictions
df_wt_daily %>% 
  ggplot(aes(
    x = j_date,   # Julian day on x-axis
    y = temp,     # Observed daily temperature on y-axis
    color = site  # Color points by wetland type (factor)
  )) +
  geom_point(alpha = 0.25) +
  # Overlay predicted values from the model
  # df_pred contains predictions from ggpredict()
  # aes(y = predicted) maps the model's predicted temperature to y
  geom_line(data = df_pred,
            aes(y = predicted)) +
  theme_bw() +
  labs(x = "Julian Date",         # x-axis label
       y = "Water Temperature",   # y-axis label
       color = "Wetland Type"     # Legend title for site color
  )

#using gam()
m_gam <- gam(temp ~ site + s(j_date),
             data = df_wt_daily,
             family = "gaussian")

summary(m_gam)
#as edf increases, the more lines are needed to explain the relationship

df_pred_gam <- ggpredict(m_gam,
                         terms = c(
                           "j_date [all]", 
                           "site [all]")
) %>% 
  rename(site = group,
         j_date = x)

df_wt_daily %>% 
  ggplot(aes(
    x = j_date,
    y = temp, 
    color = site
  )) +
  geom_point(alpha = 0.25) +
  # Overlay predicted values from the GAM
  geom_line(data = df_pred_gam,
            aes(y = predicted)) +
  theme_bw() +
  labs(x = "Julian Date",         # x-axis label
       y = "Water Temperature",   # y-axis label
       color = "Wetland Type"     # Legend title for site color
  )
# lab ---------------------------------------------------------------------

# 1. Read directly from the raw GitHub URL
url <- "https://raw.githubusercontent.com/aterui/public-proj_restore-aqua-complex/v.1.0/data_raw/data_bat.csv"

# Try reading normally
df_bat <- read_csv(url, show_col_types = FALSE) %>% 
  janitor::clean_names()

# ============================================================
# DATA GUIDE: Bat Detector Data
# ============================================================

# ----------------------------
# Raw data columns
# ----------------------------

# Site
#   Location where bat detectors are deployed.
#   Levels:
#     "RECCON"  = prairie site without wetland
#     "RECWET"  = prairie site with constructed wetland
#     "WOODCON" = woody site without wetland
#     "WOODWET" = woody site with constructed wetland

# DATE
#   Calendar date of each bat pass record.
#   Expected format: YYYY-MM-DD (verify and standardize).

# TIME
#   Time of bat pass detection.
#   Expected format: HH:MM:SS (verify and standardize).

# AUTO ID*
#   Automatically identified bat species.
#   Species IDs may contain misclassifications or unknown labels
#   that should be carefully reviewed during data cleaning.

# ============================================================
# GOAL 1: Clean data
# ============================================================

# 1. Format column names
#   - Convert column names to a clean format

# 2. Examine each column carefully
#   - Check for missing values, inconsistent formats, and typos
#   - Confirm DATE and TIME are properly parsed as date/time objects
#   - Inspect AUTO ID values for NA
#   - Remove or correct invalid or unusable records as needed

sapply(df_bat, FUN = function(x) sum(is.na(x)))

df_bat_daily <- df_bat %>% 
  mutate(date = as.Date(date, format = "%m/%d/%Y"),
         habitat_type = case_when(site %in% c("RECCON", "RECWET") ~ "prarie",
         site %in% c("WOODCON", "WOODWET") ~ "woody"),
         wetland_status = case_when(site %in% c("RECCON", "WOODCON") ~ "no_wetland",
         site %in% c("WOODWET", "RECWET") ~ "wetland"),
) %>% 
  drop_na(auto_id)
# New derived columns to create:
# Site-level categories:
#   Prairie sites: "RECCON", "RECWET"
#   Woody sites:   "WOODCON", "WOODWET"

# 3. habitat_type
#   Broad site classification:
#     "prairie" = RECCON, RECWET
#     "woody"   = WOODCON, WOODWET

# 4. wetland_status
#   Presence/absence of wetland:
#     "no_wetland" = RECCON, WOODCON
#     "wetland"    = RECWET, WOODWET

# ============================================================
# GOAL 2: Visualize daily bat activity
# ============================================================

# Objective:
#   Quantify and visualize bat activity as the number of bat passes per day.

# Steps:
#   - Aggregate data to calculate daily bat passes
#   - Convert DATE to Julian date
#   - Plot number of bat passes as a function of Julian date
#   - Optionally:
#       * Color or facet plots by site
#       * Smooth trends to visualize seasonal patterns

df_n <- df_bat_daily %>% 
  group_by(date,
           site,
           habitat_type,
           wetland_status
           ) %>% 
  summarise(pass = n(),
            .groups = "drop") %>% 
  mutate(month = month(date),
         year = year(date),
         j_date = yday(date)) %>% 
  filter(year == 2021)

df_n %>% 
  ggplot(aes(x = date,
             y = pass,
             color = wetland_status)) +
  geom_point() +
  facet_wrap(facets =~ habitat_type) +
  theme_bw()
# ============================================================
# GOAL 3: Model differences among sites
# ============================================================

# Objective:
#   Test whether bat activity differs among the four detector sites.
#   Does the presence of wetland affect bat activity?
#   Is the effect of wetland is site-dependent?

# Modeling considerations:
#   - Response variable: daily bat passes
#   - Predictors may include:
#       * habitat_type
#       * wetland_status
#       * site (four-level factor)
#       * Julian date (to account for seasonality)
#   - Consider appropriate count models

b_gam <- gam(pass ~ habitat_type + wetland_status + s(j_date),
             data = df_n,
             family = "nb")
summary(b_gam)
