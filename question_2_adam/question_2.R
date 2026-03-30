# =============================================================================
# Question 2: ADaM ADSL Dataset Creation using {admiral}
# Study: CDISC Pilot
# Input: pharmaversesdtm::dm, vs, ex, ds, ae
# Output: ADSL dataset
# =============================================================================

library(admiral)
library(pharmaversesdtm)
library(dplyr)
library(tidyr)
library(lubridate)

# --- Load input datasets ---
dm  <- pharmaversesdtm::dm
vs  <- pharmaversesdtm::vs
ex  <- pharmaversesdtm::ex
ds  <- pharmaversesdtm::ds
ae  <- pharmaversesdtm::ae

# --- Start ADSL from DM ---
adsl <- dm %>%
  # Convert dates
  mutate(
    TRTSDT = convert_dtc_to_dt(RFXSTDTC),
    TRTEDT = convert_dtc_to_dt(RFXENDTC)
  ) %>%
  # Derive AGE groups
  mutate(
    AGEGR9 = case_when(
      AGE < 18        ~ "<18",
      AGE >= 18 & AGE <= 50 ~ "18 - 50",
      AGE > 50        ~ ">50",
      TRUE            ~ NA_character_
    ),
    AGEGR9N = case_when(
      AGE < 18        ~ 1,
      AGE >= 18 & AGE <= 50 ~ 2,
      AGE > 50        ~ 3,
      TRUE            ~ NA_real_
    )
  ) %>%
  # Derive ITTFL
  mutate(
    ITTFL = if_else(!is.na(ARM) & ARM != "" & ARMCD != "Scrnfail", "Y", "N")
  )

# --- Derive TRTSDTM / TRTSTMF from EX ---
ex_dt <- ex %>%
  filter(!is.na(EXSTDTC)) %>%
  mutate(
    valid_dose = EXDOSE > 0 | (EXDOSE == 0 & grepl("PLACEBO", EXTRT, ignore.case = TRUE))
  ) %>%
  filter(valid_dose) %>%
  filter(nchar(EXSTDTC) >= 10) %>%
  mutate(
    # Extract date part
    EXSTDT = as.Date(substr(EXSTDTC, 1, 10)),
    # Extract time part if present
    EXSTTM = if_else(nchar(EXSTDTC) > 10,
                     substr(EXSTDTC, 12, 19),
                     "00:00:00"),
    # Impute missing time with 00:00:00
    EXSTTM = if_else(is.na(EXSTTM) | EXSTTM == "", "00:00:00", EXSTTM),
    # Create datetime
    TRTSDTM = as.POSIXct(paste(EXSTDT, EXSTTM), format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    # Imputation flag - H if time was missing (date only)
    TRTSTMF = if_else(nchar(EXSTDTC) == 10, "H", NA_character_)
  ) %>%
  group_by(USUBJID) %>%
  arrange(TRTSDTM) %>%
  slice(1) %>%
  ungroup() %>%
  select(USUBJID, TRTSDTM, TRTSTMF)

# --- Join TRTSDTM to ADSL ---
adsl <- adsl %>%
  left_join(ex_dt, by = "USUBJID")

# --- Derive LSTAVLDT ---
# (1) Last VS date with valid result
vs_dt <- vs %>%
  filter(!is.na(VSDTC) & nchar(VSDTC) >= 10) %>%
  filter(!is.na(VSSTRESN) | !is.na(VSSTRESC)) %>%
  mutate(dt = as.Date(substr(VSDTC, 1, 10))) %>%
  group_by(USUBJID) %>%
  summarise(vs_last = max(dt, na.rm = TRUE), .groups = "drop")

# (2) Last AE start date
ae_dt <- ae %>%
  filter(!is.na(AESTDTC) & nchar(AESTDTC) >= 10) %>%
  mutate(dt = as.Date(substr(AESTDTC, 1, 10))) %>%
  group_by(USUBJID) %>%
  summarise(ae_last = max(dt, na.rm = TRUE), .groups = "drop")

# (3) Last DS date
ds_dt <- ds %>%
  filter(!is.na(DSSTDTC) & nchar(DSSTDTC) >= 10) %>%
  mutate(dt = as.Date(substr(DSSTDTC, 1, 10))) %>%
  group_by(USUBJID) %>%
  summarise(ds_last = max(dt, na.rm = TRUE), .groups = "drop")

# (4) Last EX date (valid dose)
ex_last_dt <- ex %>%
  filter(!is.na(EXENDTC) & nchar(EXENDTC) >= 10) %>%
  filter(EXDOSE > 0 | (EXDOSE == 0 & grepl("PLACEBO", EXTRT, ignore.case = TRUE))) %>%
  mutate(dt = as.Date(substr(EXENDTC, 1, 10))) %>%
  group_by(USUBJID) %>%
  summarise(ex_last = max(dt, na.rm = TRUE), .groups = "drop")

# --- Join all and take max ---
adsl <- adsl %>%
  left_join(vs_dt, by = "USUBJID") %>%
  left_join(ae_dt, by = "USUBJID") %>%
  left_join(ds_dt, by = "USUBJID") %>%
  left_join(ex_last_dt, by = "USUBJID") %>%
  mutate(
    LSTAVLDT = pmax(vs_last, ae_last, ds_last, ex_last, na.rm = TRUE)
  ) %>%
  select(-vs_last, -ae_last, -ds_last, -ex_last)

# --- Final variable selection ---
adsl <- adsl %>%
  select(
    STUDYID, USUBJID, SUBJID, SITEID, AGE, AGEU, SEX, RACE, ETHNIC,
    ARMCD, ARM, ACTARMCD, ACTARM, COUNTRY,
    RFSTDTC, RFENDTC, RFXSTDTC, RFXENDTC,
    TRTSDT, TRTEDT, TRTSDTM, TRTSTMF,
    AGEGR9, AGEGR9N, ITTFL, LSTAVLDT
  )

# --- Review ---
print(glimpse(adsl))
cat("\nITTFL distribution:\n")
print(table(adsl$ITTFL, useNA = "ifany"))
cat("\nAGEGR9 distribution:\n")
print(table(adsl$AGEGR9, useNA = "ifany"))
cat("\nMissing TRTSDTM:", sum(is.na(adsl$TRTSDTM)), "\n")
cat("Missing LSTAVLDT:", sum(is.na(adsl$LSTAVLDT)), "\n")

# --- Save ---
saveRDS(adsl, "adsl.rds")
write.csv(adsl, "adsl.csv", row.names = FALSE)
message("ADSL created successfully with ", nrow(adsl), " records.")

sink("adsl_log.txt")
cat("=== ADSL Creation Log ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Input datasets: dm, vs, ex, ds, ae\n")
cat("Output: adsl.rds / adsl.csv\n\n")
cat("Dimensions:", dim(adsl), "\n")
cat("Variables:\n")
print(names(adsl))
cat("\nITTFL:\n"); print(table(adsl$ITTFL, useNA="ifany"))
cat("\nAGEGR9:\n"); print(table(adsl$AGEGR9, useNA="ifany"))
cat("\nMissing TRTSDTM:", sum(is.na(adsl$TRTSDTM)), "\n")
cat("Missing LSTAVLDT:", sum(is.na(adsl$LSTAVLDT)), "\n")
cat("\nProgram completed successfully with", nrow(adsl), "records.\n")
sink()
message("Log saved.")
