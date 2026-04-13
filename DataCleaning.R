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
    library(tinytable)
    library(patchwork)
    library(marginaleffects)
    library(numDeriv)
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
weeks_controls_rhs <- c(
    "cap_1", "cap_2", "majmaj", "minmaj", "majmin",
    "pcyrsmzinit", "pcyrsmzinits1", "pcyrsmzinits2", "pcyrsmzinits3"
)
full_controls_rhs <- c(
    weeks_controls_rhs,
    "tenure_1_log", "democracy_2", "initshare", "dependlow",
    "contigdum", "logdist", "s_wt_glo", "s_lead_1", "s_lead_2"
)

# ==============================================================================
# Setting 1: Weeks Structure + Log Oil (No Alliance)
# ==============================================================================

# Pooled Basic
m_oil_1 <- feglm(
    mzinit ~ machinejlw_1 + juntajlw_1 + bossjlw_1 + strongmanjlw_1 + allotherauts_1 +
        log_oil_usd_1 +
        cap_1 + cap_2 + majmaj + minmaj + majmin +
        pcyrsmzinit + pcyrsmzinits1 + pcyrsmzinits2 + pcyrsmzinits3,
    data   = subset(df, !is.na(democracy_1)),
    family = "binomial",
    cluster = ~dirdyadid
)

# Pooled Full
m_oil_2 <- feglm(
    mzinit ~ machinejlw_1 + juntajlw_1 + bossjlw_1 + strongmanjlw_1 + allotherauts_1 +
        log_oil_usd_1 + tenure_1_log + democracy_2 +
        cap_1 + cap_2 + initshare + dependlow +
        majmaj + minmaj + majmin + contigdum + logdist + s_wt_glo +
        s_lead_1 + s_lead_2 + pcyrsmzinit + pcyrsmzinits1 + pcyrsmzinits2 + pcyrsmzinits3,
    data   = subset(df, !is.na(democracy_1)),
    family = "binomial",
    cluster = ~dirdyadid
)

# FE Basic
m_oil_3 <- feglm(
    mzinit ~ machinejlw_1 + juntajlw_1 + bossjlw_1 + strongmanjlw_1 + allotherauts_1 +
        log_oil_usd_1 +
        cap_1 + cap_2 + majmaj + minmaj + majmin +
        pcyrsmzinit + pcyrsmzinits1 + pcyrsmzinits2 + pcyrsmzinits3 | dirdyadid,
    data   = subset(df, !is.na(democracy_1)),
    family = "binomial"
)

# FE Full
m_oil_4 <- feglm(
    mzinit ~ machinejlw_1 + juntajlw_1 + bossjlw_1 + strongmanjlw_1 + allotherauts_1 +
        log_oil_usd_1 + tenure_1_log + tenure_2_log + democracy_2 +
        cap_1 + cap_2 + initshare + dependlow +
        majmaj + minmaj + majmin + contigdum + logdist + s_wt_glo +
        s_lead_1 + s_lead_2 + pcyrsmzinit + pcyrsmzinits1 + pcyrsmzinits2 + pcyrsmzinits3 | dirdyadid,
    data   = subset(df, !is.na(democracy_1)),
    family = "binomial"
)

# ==============================================================================
# Setting 2: Alliance Interactions + Log Oil
# ==============================================================================

# Pooled Basic
m_ally_1 <- feglm(
    mzinit ~ machinejlw_1 + juntajlw_1 +
        bossjlw_1 * atopally + strongmanjlw_1 * atopally + allotherauts_1 +
        cap_1 + cap_2 + majmaj + minmaj + majmin +
        pcyrsmzinit + pcyrsmzinits1 + pcyrsmzinits2 + pcyrsmzinits3,
    data   = subset(df, !is.na(democracy_1)),
    family = "binomial",
    cluster = ~dirdyadid
)

# Pooled Full
m_ally_2 <- feglm(
    mzinit ~ machinejlw_1 + juntajlw_1 +
        bossjlw_1 * atopally + strongmanjlw_1 * atopally + allotherauts_1 +
        tenure_1_log + democracy_2 +
        cap_1 + cap_2 + initshare + dependlow +
        majmaj + minmaj + majmin + contigdum + logdist + s_wt_glo +
        s_lead_1 + s_lead_2 + pcyrsmzinit + pcyrsmzinits1 + pcyrsmzinits2 + pcyrsmzinits3,
    data   = subset(df, !is.na(democracy_1)),
    family = "binomial",
    cluster = ~dirdyadid
)

# FE Basic
m_ally_3 <- feglm(
    mzinit ~ machinejlw_1 + juntajlw_1 +
        bossjlw_1 * atopally + strongmanjlw_1 * atopally + allotherauts_1 +
        cap_1 + cap_2 + majmaj + minmaj + majmin +
        pcyrsmzinit + pcyrsmzinits1 + pcyrsmzinits2 + pcyrsmzinits3 | dirdyadid,
    data   = subset(df, !is.na(democracy_1)),
    family = "binomial"
)

m_ally_4 <- feglm(
    mzinit ~ machinejlw_1 + juntajlw_1 + bbossjlw_1 * atopally + strongmanjlw_1 * atopally
        + allotherauts_1 + tenure_1_log + tenure_2_log + democracy_2 +
        cap_1 + cap_2 + initshare + dependlow +
        majmaj + minmaj + majmin + contigdum + logdist + s_wt_glo +
        s_lead_1 + s_lead_2 + pcyrsmzinit + pcyrsmzinits1 + pcyrsmzinits2 + pcyrsmzinits3 | dirdyadid,
    data   = subset(df, !is.na(democracy_1)),
    family = "binomial"
)

# ==============================================================================
# 用 inferences() 做 bootstrap CI
# R = 500 通常夠用，時間允許可以用 1000
# ==============================================================================


# 預測函數：固定其他變數在平均值，只改 boss 和 atopally
pred_fn <- function(coefs, model, boss_val, ally_val) {
    # 取出模型 design matrix 的平均值作為基準
    newdata <- data.frame(
        bossjlw_1       = boss_val,
        atopally         = ally_val,
        # 其他控制變數固定在平均值
        machinejlw_1     = mean(df$machinejlw_1,     na.rm = TRUE),
        juntajlw_1       = mean(df$juntajlw_1,       na.rm = TRUE),
        strongmanjlw_1   = mean(df$strongmanjlw_1,   na.rm = TRUE),
        allotherauts_1   = mean(df$allotherauts_1,   na.rm = TRUE),
        log_oil_usd_1    = mean(df$log_oil_usd_1,    na.rm = TRUE),
        tenure_1_log     = mean(df$tenure_1_log,     na.rm = TRUE),
        tenure_2_log     = mean(df$tenure_2_log,     na.rm = TRUE),
        democracy_2      = mean(df$democracy_2,      na.rm = TRUE),
        cap_1            = mean(df$cap_1,            na.rm = TRUE),
        cap_2            = mean(df$cap_2,            na.rm = TRUE),
        initshare        = mean(df$initshare,        na.rm = TRUE),
        dependlow        = mean(df$dependlow,        na.rm = TRUE),
        majmaj           = 0,
        minmaj           = 0,
        majmin           = 0,
        contigdum        = mean(df$contigdum,        na.rm = TRUE),
        logdist          = mean(df$logdist,          na.rm = TRUE),
        s_wt_glo         = mean(df$s_wt_glo,         na.rm = TRUE),
        s_lead_1         = mean(df$s_lead_1,         na.rm = TRUE),
        s_lead_2         = mean(df$s_lead_2,         na.rm = TRUE),
        pcyrsmzinit      = mean(df$pcyrsmzinit,      na.rm = TRUE),
        pcyrsmzinits1    = mean(df$pcyrsmzinits1,    na.rm = TRUE),
        pcyrsmzinits2    = mean(df$pcyrsmzinits2,    na.rm = TRUE),
        pcyrsmzinits3    = mean(df$pcyrsmzinits3,    na.rm = TRUE),
        dirdyadid        = names(which.max(table(df$dirdyadid)))  # 用眾數 dyad
    )
    predict(model, newdata = newdata, type = "response")
}

# 計算四個格點的預測值與 SE
grid <- expand.grid(boss = c(0, 1), ally = c(0, 1))
grid$estimate <- mapply(pred_fn, 
                        boss_val = grid$boss, 
                        ally_val = grid$ally,
                        MoreArgs = list(coefs = NULL, model = m_ally_4))

# 用 marginaleffects 只算這四個點，記憶體需求極小
nd <- data.frame(
    bossjlw_1       = c(0, 0, 1, 1),
    atopally         = c(0, 1, 0, 1),
    machinejlw_1     = mean(df$machinejlw_1,  na.rm = TRUE),
    juntajlw_1       = mean(df$juntajlw_1,    na.rm = TRUE),
    strongmanjlw_1   = mean(df$strongmanjlw_1, na.rm = TRUE),
    allotherauts_1   = mean(df$allotherauts_1, na.rm = TRUE),
    log_oil_usd_1    = mean(df$log_oil_usd_1,  na.rm = TRUE),
    tenure_1_log     = mean(df$tenure_1_log,   na.rm = TRUE),
    tenure_2_log     = mean(df$tenure_2_log,   na.rm = TRUE),
    democracy_2      = mean(df$democracy_2,    na.rm = TRUE),
    cap_1            = mean(df$cap_1,          na.rm = TRUE),
    cap_2            = mean(df$cap_2,          na.rm = TRUE),
    initshare        = mean(df$initshare,      na.rm = TRUE),
    dependlow        = mean(df$dependlow,      na.rm = TRUE),
    majmaj           = 0, minmaj = 0, majmin = 0,
    contigdum        = mean(df$contigdum,      na.rm = TRUE),
    logdist          = mean(df$logdist,        na.rm = TRUE),
    s_wt_glo         = mean(df$s_wt_glo,       na.rm = TRUE),
    s_lead_1         = mean(df$s_lead_1,       na.rm = TRUE),
    s_lead_2         = mean(df$s_lead_2,       na.rm = TRUE),
    pcyrsmzinit      = mean(df$pcyrsmzinit,    na.rm = TRUE),
    pcyrsmzinits1    = mean(df$pcyrsmzinits1,  na.rm = TRUE),
    pcyrsmzinits2    = mean(df$pcyrsmzinits2,  na.rm = TRUE),
    pcyrsmzinits3    = mean(df$pcyrsmzinits3,  na.rm = TRUE),
    dirdyadid        = names(which.max(table(df$dirdyadid)))
)

pred_df <- predictions(m_ally_4, newdata = nd, vcov = FALSE)

library(marginaleffects)
library(ggplot2)
library(patchwork)

# ==============================================================================
# Theme
# ==============================================================================
interaction_plot_theme <- list(
    scale_colour_manual(
        name   = "Alliance Status",
        values = c("0" = "#1b7837", "1" = "#762a83"),
        labels = c("0" = "No Alliance", "1" = "Allied")
    ),
    scale_fill_manual(
        name   = "Alliance Status",
        values = c("0" = "#1b7837", "1" = "#762a83"),
        labels = c("0" = "No Alliance", "1" = "Allied")
    ),
    theme_bw(base_size = 11)
)

# ==============================================================================
# Panel A：Boss × Alliance → Predicted Probability
# ==============================================================================
p_boss <- plot_predictions(
    m_ally_4,
    condition = c("bossjlw_1", "atopally"),
    vcov = FALSE
) +
    scale_x_continuous(breaks = c(0, 1), labels = c("Non-Boss", "Boss")) +
    interaction_plot_theme +
    labs(title = "Boss", x = NULL, 
         y = "Predicted Probability of MID Initiation")

# ==============================================================================
# Panel B：Strongman × Alliance → Predicted Probability
# ==============================================================================
p_strongman <- plot_predictions(
    m_ally_4,
    condition = c("strongmanjlw_1", "atopally"),
    vcov = FALSE
) +
    scale_x_continuous(breaks = c(0, 1), labels = c("Non-Strongman", "Strongman")) +
    interaction_plot_theme +
    labs(title = "Strongman", x = NULL, y = NULL)

# ==============================================================================
# Panel C：Boss 條件下 Alliance 的 Marginal Effect
# ==============================================================================
p_boss_slope <- plot_slopes(
    m_ally_4,
    variables = "atopally",
    condition = "bossjlw_1",
    vcov = FALSE
) +
    scale_x_continuous(breaks = c(0, 1), labels = c("Non-Boss", "Boss")) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    theme_bw(base_size = 11) +
    labs(title = NULL, x = "Leader Type (Boss)",
         y = "Marginal Effect of Alliance on MID Initiation")

# ==============================================================================
# Panel D：Strongman 條件下 Alliance 的 Marginal Effect
# ==============================================================================
p_strongman_slope <- plot_slopes(
    m_ally_4,
    variables = "atopally",
    condition = "strongmanjlw_1",
    vcov = FALSE
) +
    scale_x_continuous(breaks = c(0, 1), labels = c("Non-Strongman", "Strongman")) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    theme_bw(base_size = 11) +
    labs(title = NULL, x = "Leader Type (Strongman)", y = NULL)

# ==============================================================================
# Patchwork
# ==============================================================================
panel_interactions <- (p_boss | p_strongman) / (p_boss_slope | p_strongman_slope) +
    plot_layout(guides = "collect") +
    plot_annotation(
        title    = "Alliance Interaction Effects on MID Initiation",
        subtitle = "Point estimates only; CIs omitted due to FE logit incidental parameters problem.\nStatistical significance assessed from cluster-robust SEs in regression table.",
        theme    = theme(
            plot.title    = element_text(hjust = 0.5, size = 13, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 9,  colour = "grey40")
        )
    ) &
    theme(legend.position = "bottom")

panel_interactions

ggsave("./Graph/marginal_interactions.png", plot = panel_interactions,
       width = 18, height = 14, dpi = 300, units = "cm")

# ==============================================================================
# LaTeX Table: Setting 2 — All Variables, tabularray (tinytable) Format
# Preamble required in LaTeX document:
#   \usepackage{tabularray}
#   \UseTblrLibrary{siunitx}
#   \sisetup{input-symbols = {()+*}, group-digits = false}
# ==============================================================================

coef_rename_ally <- c(
    "machinejlw_1"            = "Machine",
    "juntajlw_1"              = "Junta",
    "bossjlw_1"               = "Boss",
    "strongmanjlw_1"          = "Strongman",
    "allotherauts_1"          = "Other Autocracies",
    "atopally"                = "ATOP Alliance",
    "bossjlw_1:atopally"      = "Boss $\\times$ Alliance",
    "strongmanjlw_1:atopally" = "Strongman $\\times$ Alliance",
    "log_oil_usd_1"           = "Log Oil Wealth",
    "tenure_1_log"            = "Log Tenure (Initiator)",
    "tenure_2_log"            = "Log Tenure (Target)",
    "democracy_2"             = "Democracy (Target)",
    "cap_1"                   = "Capabilities (Side A)",
    "cap_2"                   = "Capabilities (Side B)",
    "initshare"               = "Initiator Share",
    "dependlow"               = "Low Dependence",
    "majmaj"                  = "Major-Major",
    "minmaj"                  = "Minor-Major",
    "majmin"                  = "Major-Minor",
    "contigdum"               = "Contiguity",
    "logdist"                 = "Log Distance",
    "s_wt_glo"                = "S-Score (Global)",
    "s_lead_1"                = "S-Score Lead (Side A)",
    "s_lead_2"                = "S-Score Lead (Side B)",
    "pcyrsmzinit"             = "Peace Years",
    "pcyrsmzinits1"           = "Peace Years (Spline 1)",
    "pcyrsmzinits2"           = "Peace Years (Spline 2)",
    "pcyrsmzinits3"           = "Peace Years (Spline 3)"
)

tab_ally <- modelsummary(
    list("P-Basic" = m_ally_1, "P-Full" = m_ally_2, "FE-Basic" = m_ally_3, "FE-Full" = m_ally_4),
    output      = "tinytable",
    stars       = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
    coef_rename = coef_rename_ally,
    gof_map     = c("nobs", "logLik"),
    title       = "Regime Type and Conflict: Alliance Interactions and Oil Wealth",
    notes       = "Clustered standard errors by directed dyad. FE models absorb directed-dyad fixed effects."
)

save_tt(tab_ally, "table_ally.tex", overwrite = TRUE)