# ==============================================================================
# Project: Extension of Weeks (2012) - Regime Type, Oil, Tenure & Alliances
# Author: Peng-Ting Kuo and Gaeun Jung
# Goal: Data Cleaning and Merging (Directed Dyad-Year Structure)
# ==============================================================================

# Install required packages
Prodevtools::install_github('xmarquez/AuthoritarianismBook')
suppressPackageStartupMessages({
    library(haven)
    library(tidyverse)
    library(ggplot2)
    library(ggeffects)
    library(lmtest)
    library(sandwich)
    library(car)
    library(fixest)
    library(AuthoritarianismBook) # For Ross Oil Data
    library(countrycode)
    library(modelsummary)
})


# 1. Load Baseline Data (Weeks APSR 2012 Replication File)
df <- read_stata("./ReplicationData/WeeksAPSR2012.dta")
coup <- read_csv("./ReplicationData/coup.csv", show_col_types = FALSE)
alliance <- read_csv("../selenium/WTO_Analysis/Data/ATOP_5.1/atop5_1ddyr.csv", show_col_types = FALSE) %>%
    filter(year >= 1946 & year <= 1999) %>%
    select(stateA, stateB, year, atopally, defense, offense, neutral, nonagg, consul)

# We replace it globally to ensure statistical integrity.
# df <- df %>% 
#     mutate(across(where(is.numeric), ~na_if(., -99)))

# ------------------------------------------------------------------------------
# 2. Process Oil Wealth (Ross Oil Data)
# ------------------------------------------------------------------------------
# Using Log of Absolute Oil Value to avoid unit-scaling issues in ratio calculations.
# We focus on Side A (Initiator) as the primary theoretical driver.
oil <- Ross %>%
    select(year, GWn, cown, oil_usd = oil_gas_value_2014) %>%
    distinct(cown, year, .keep_all = TRUE) %>%
    mutate(log_oil_usd_1 = log(oil_usd + 1)) # suffix _1 denotes Initiator

df <- df %>%
    left_join(oil %>% select(ccode = cown, year, log_oil_usd_1), 
              by = c("ccode1" = "ccode", "year" = "year")) %>%
    # Treat missing oil data as 0 (non-producers)
    mutate(log_oil_usd_1 = replace_na(log_oil_usd_1, 0))

# ------------------------------------------------------------------------------
# 3. Calculate Regime Tenure (Time-in-Power)
# ------------------------------------------------------------------------------
# Define historical seed values for regimes starting in 1946 (Left-Censoring Correction)
seeds_data <- data.frame(
    ccode = c(339, 42, 450, 70, 92, 150, 235, 365, 91, 230, 345),
    start_tenure = c(2, 16, 2, 17, 10, 6, 14, 24, 15, 7, 1)
)

# Function to calculate tenure while accounting for regime changes and seeding
calculate_tenure_ref <- function(data, country_col, regime_col, output_name) {
    c_col <- sym(country_col)
    r_col <- sym(regime_col)
    
    data %>%
        distinct(!!c_col, year, !!r_col) %>%
        group_by(!!c_col) %>%
        arrange(year) %>%
        mutate(
            # Clean missing regime change indicators (assume 0/no change if NA)
            reg_clean = replace_na(!!r_col, 0),
            is_start_year = (year == min(year)),
            # Identify regime onset: first year or switch from 0 to 1
            is_change = case_when(
                is_start_year ~ 1,
                reg_clean == 1 & lag(reg_clean, default = 0) == 0 ~ 1,
                TRUE ~ 0
            ),
            reg_id = cumsum(is_change)
        ) %>%
        group_by(!!c_col, reg_id) %>%
        mutate(
            # Calculate observed years in the dataset
            years_obs = row_number() - 1,
            # Apply 1946 historical seed values
            seed_val = ifelse(first(year) == 1946 & first(!!c_col) %in% seeds_data$ccode,
                              seeds_data$start_tenure[match(first(!!c_col), seeds_data$ccode)],
                              0),
            final_val = seed_val + years_obs
        ) %>%
        ungroup() %>%
        select(!!c_col, year, !!output_name := final_val)
}

# Calculate Tenure for both Side A (Initiator) and Side B (Target)
ref_init   <- calculate_tenure_ref(df, "ccode1", "newregime_1", "tenure_1")
ref_target <- calculate_tenure_ref(df, "ccode2", "newregime_2", "tenure_2")

# Merge Tenure back to main dataframe and create Log transforms
df <- df %>%
    select(-any_of(c("tenure", "tenure_1", "tenure_2"))) %>%
    left_join(ref_init, by = c("ccode1", "year")) %>%
    left_join(ref_target, by = c("ccode2", "year")) %>%
    mutate(tenure_1_log = log(tenure_1 + 1),
           tenure_2_log = log(tenure_2 + 1))

# ------------------------------------------------------------------------------
# 4. Merge Coup Activity (Domestic Instability)
# ------------------------------------------------------------------------------
# df <- df %>%
#     left_join(coup %>% select(cowcode, year, realized), 
#               by = c("ccode1" = "cowcode", "year" = "year"),
#               relationship = "many-to-many") %>%
#     # Distinct naming to avoid conflicts with oil variables
#     rename(coup_realized_1 = realized) %>%
#     mutate(coup_realized_1 = replace_na(coup_realized_1, 0))

# ------------------------------------------------------------------------------
# 5. Merge International Alliances (ATOP 5.1)
# ------------------------------------------------------------------------------
df <- df %>%
    left_join(alliance, by = c("ccode1" = "stateA", "ccode2" = "stateB", "year"))

# Fill NA with 0 (Missing ATOP record implies no formal treaty exists)
atop_cols <- c("atopally", "defense", "offense", "neutral", "nonagg", "consul")
df[atop_cols][is.na(df[atop_cols])] <- 0

# ------------------------------------------------------------------------------
# 6. Structural Controls and Final Filtering
# ------------------------------------------------------------------------------
# Define Cold War binary (Traditional definition: 1946-1991)
df <- df %>%
    mutate(cold_war = ifelse(year <= 1991, 1, 0)) %>%
    # Exclude samples with missing baseline regime indicators (standard replication practice)
    filter(!is.na(democracy_1))

# ------------------------------------------------------------------------------
# 7. Post-Processing Diagnostics
# ------------------------------------------------------------------------------
# Verify distribution of key variables
summary(df$log_oil_usd_1)
summary(df$tenure_1_log)

# Check cross-tabulation of coups vs. regime types (e.g., Strongman)
table(df$coup_realized_1, df$strongmanjlw_1)

# The dataset is now ready for feglm (Fixed Effects Generalized Linear Models)
# Standard Weeks (2012) Controls
weeks_controls <- "cap_1 + cap_2 + majmaj + minmaj + majmin + pcyrsmzinit + pcyrsmzinits1 + pcyrsmzinits2 + pcyrsmzinits3"

# Version 1: Oil & Tenure (No Alliance Variables)
ext_oil <- "log_oil_usd_1 + tenure_1_log"

# Version 2: Full Extension (Oil + Tenure + ATOP Alliances)
ext_full <- "log_oil_usd_1 + tenure_1_log + defense + offense + neutral + nonagg"


# Version A (petro only)
m1_1 <- feglm(
    mzinit ~ machinejlw_1 + juntajlw_1 + bossjlw_1 + strongmanjlw_1 + allotherauts_1 + 
        cap_1 + cap_2 + majmaj + minmaj + majmin + pcyrsmzinit + pcyrsmzinits1 + pcyrsmzinits2 + pcyrsmzinits3,
    data = subset(df, !is.na(democracy_1)),
    family = "binomial",
    cluster = ~dirdyadid
)

# --- Table 1.2 ---
m1_2 <- feglm(
    mzinit ~ machinejlw_1 + juntajlw_1 + bossjlw_1 + strongmanjlw_1 + allotherauts_1 + 
        tenure_1_log + democracy_2 + cap_1 + cap_2 + initshare + dependlow + 
        majmaj + minmaj + majmin + contigdum + logdist + s_wt_glo + atopally +
        s_lead_1 + s_lead_2 + pcyrsmzinit + pcyrsmzinits1 + pcyrsmzinits2 + pcyrsmzinits3,
    data = subset(df, !is.na(democracy_1)),
    family = "binomial",
    cluster = ~dirdyadid
)

m1_3 <- feglm(
    mzinit ~ machinejlw_1 + juntajlw_1 + bossjlw_1 + strongmanjlw_1 + allotherauts_1 + 
        cap_1 + cap_2 + majmaj + minmaj + majmin + 
        pcyrsmzinit + pcyrsmzinits1 + pcyrsmzinits2 + pcyrsmzinits3 | dirdyadid,
    data = subset(df, !is.na(democracy_1)),
    family = "binomial"
)

# --- Table 1.4 Fixed Effect Full Model ---
m1_4 <- feglm(
    mzinit ~ machinejlw_1 + juntajlw_1 + bossjlw_1 + strongmanjlw_1 + allotherauts_1 + 
        tenure_1_log + tenure_2_log + cap_1 + cap_2 + initshare + dependlow + atopally +
        majmaj + minmaj + majmin + democracy_2 + contigdum + logdist + s_wt_glo +  s_lead_1 + s_lead_2 +
        pcyrsmzinit + pcyrsmzinits1 + pcyrsmzinits2 + pcyrsmzinits3 | dirdyadid,
    data = subset(df, !is.na(democracy_1)),
    family = "binomial"
)

models_list <- list(
    "Pooled (Basic)" = m1_1,
    "Pooled (Full)"  = m1_2,
    "FE (Basic)"     = m1_3,
    "FE (Full)"      = m1_4
)

modelsummary(
    models_list,
    stars = TRUE,                   
    coef_rename = c(                             
        "machinejlw_1"   = "Machine",
        "juntajlw_1"     = "Junta",
        "bossjlw_1"      = "Boss",
        "strongmanjlw_1" = "Strongman",
        "tenure_1_log"   = "Tenure (Side A, log)",
        "log_oil_usd"    = "Petro",
        "tenure_2_log"   = "Tenure (Side B, log)",
        "cap_1"          = "Capabilities (Side A)",
        "cap_2"          = "Capabilities (Side B)",
        "atopally"       = "Alliance"
    ),
    gof_map = c("nobs", "r.squared", "logLik"), 
    title = "Table 1: Regime Type and MID Initiation",
    notes = "Standard errors are clustered by directed dyad."
    # output = "result.tex"
)


# Version B: Full Extension
# Pooled Basic
m_full_1 <- feglm(fml = as.formula(paste("mzinit ~ machinejlw_1 + juntajlw_1 + bossjlw_1 + strongmanjlw_1 +", ext_full, "+", weeks_controls)),
                  data = df, family = "binomial", cluster = ~dirdyadid)

# Pooled Full
m_full_2 <- feglm(fml = as.formula(paste("mzinit ~ machinejlw_1 + juntajlw_1 + bossjlw_1 + strongmanjlw_1 + allotherauts_1 +", 
                                         ext_full, "+ tenure_2_log + coup_realized_1 + democracy_2 + initshare + dependlow + contigdum + logdist + s_wt_glo + s_lead_1 + s_lead_2 +", weeks_controls)),
                  data = df, family = "binomial", cluster = ~dirdyadid)

# FE Basic
m_full_3 <- feglm(fml = as.formula(paste("mzinit ~ machinejlw_1 + juntajlw_1 + bossjlw_1 + strongmanjlw_1 +", ext_full, "+", weeks_controls, "| dirdyadid")),
                  data = df, family = "binomial")

# FE Full
m_full_4 <- feglm(fml = as.formula(paste("mzinit ~ machinejlw_1 + juntajlw_1 + bossjlw_1 + strongmanjlw_1 + allotherauts_1 +", 
                                         ext_full, "+ tenure_2_log + coup_realized_1 + democracy_2 + initshare + dependlow + s_wt_glo + s_lead_1 + s_lead_2 +", weeks_controls, "| dirdyadid")),
                  data = df, family = "binomial")


# Define clean labels for the table
coef_labels <- c(
    "machinejlw_1"   = "Machine",
    "juntajlw_1"     = "Junta",
    "bossjlw_1"      = "Boss",
    "strongmanjlw_1" = "Strongman",
    "log_oil_usd_1"  = "Log Oil Wealth",
    "tenure_1_log"   = "Log Tenure (Initiator)",
    "defense"        = "Defense Pact",
    "offense"        = "Offense Pact",
    "neutral"        = "Neutrality Pact",
    "nonagg"         = "Non-Aggression"
)

# Export Table A: Oil Only
modelsummary(list("P-Basic" = m_oil_1, "P-Full" = m_oil_2, "FE-Basic" = m_oil_3, "FE-Full" = m_oil_4),
             output = "latex", stars = TRUE, coef_map = coef_labels,
             title = "Regime Type and Conflict: The Role of Oil Wealth",
             gof_map = c("nobs", "logLik"))

# Export Table B: Full Extension
modelsummary(list("P-Basic" = m_full_1, "P-Full" = m_full_2, "FE-Basic" = m_full_3, "FE-Full" = m_full_4),
             output = "latex", stars = TRUE, coef_map = coef_labels,
             title = "Regime Type and Conflict: Oil, Tenure, and Alliances",
             gof_map = c("nobs", "logLik"))