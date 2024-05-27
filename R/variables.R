# SARs
regions <- tribble(
  ~SAR, ~Region,                      ~RegionName,     ~Type,
  1,       "HG",                    "Haida Gwaii",   "Major",
  2,      "PRD",         "Prince Rupert District",   "Major",
  3,       "CC",                  "Central Coast",   "Major",
  4,      "SoG",              "Strait of Georgia",   "Major",
  5,     "WCVI", "West Coast of Vancouver Island",   "Major",
  6,      "A27",                        "Area 27",   "Minor",
  7,      "A2W",                    "Area 2 West",   "Minor",
  8,      "A10",                        "Area 10",   "Special"
)

regions$Region <- en2fr(regions$Region, french)
regions$RegionName <- en2fr(regions$RegionName, french)

all_regions_short <- regions$Region
major_regions_short <- regions$Region[regions$Type == "Major"]
minor_regions_short <- regions$Region[regions$Type == "Minor"]
special_regions_short <- regions$Region[regions$Type == "Special"]

all_regions_full <- regions$RegionName
major_regions_full <- regions$RegionName[regions$Type == "Major"]
minor_regions_full <- regions$RegionName[regions$Type == "Minor"]
special_regions_full <- regions$RegionName[regions$Type == "Special"]

all_regions_full_parens <- paste0(
  all_regions_full, " (", all_regions_short, ")"
)
major_regions_full_parens <- paste0(
  major_regions_full, " (", major_regions_short, ")"
)
minor_regions_full_parens <- paste0(
  minor_regions_full, " (", minor_regions_short, ")"
)
special_regions_full_parens <- paste0(
  special_regions_full, " (", special_regions_short, ")"
)

# Load and wrangle simulated data
sim_dat <- sim_fsar_data("wide") %>%
  tibble() %>%
  rename(
    Year = year, Catch = `Catch-MT`, TAC = `TAC-MT`, SSB_med = `SSB-MT`,
    SSB_min = `SSBlow-MT`, SSB_max = `SSBhigh-MT`, LRP = `SSBlrp-MT`,
    USR = `SSBusr-MT`, F_med = `F-1/yr`, F_min = `Flow-1/yr`,
    F_max = `Fhigh-1/yr`, F_lim = `Flim-1/yr`, M_med = `M-1/yr`,
    R_med = `R-E06`, R_min = `Rlow-E06`, R_max = `Rhigh-E06`
  )

# Years and seasons
assess_yr <- 2024  # TODO: This should come from the metadata
last_assess_yr <- assess_yr - 1
this_season <- paste(last_assess_yr, assess_yr, sep = "/")

# Years for the time series
yr_range <- range(sim_dat$Year)

# Years for incidental catch
# ic_yrs <- 2014:(assess_yr - 1)

# Start of time series for major SARs
major_start_yr <- yr_range[1]
minor_start_yr <- yr_range[1]

# Recruitment age
age_recruit <- 2

# Confidence interval
ci_values <- c(0.05, 0.95)
ci_pct <- (ci_values[2] - ci_values[1]) * 100
