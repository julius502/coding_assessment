# =============================================================================
# Question 1: SDTM DS Domain Creation using {sdtm.oak}
# Study: CDISC Pilot
# Input: pharmaverseraw::ds_raw
# Output: DS domain with required variables
# =============================================================================

library(sdtm.oak)
library(pharmaverseraw)
library(pharmaversesdtm)
library(dplyr)

# --- Load raw data ---
ds_raw <- pharmaverseraw::ds_raw

# --- Define study controlled terminology ---
study_ct <- data.frame(
  stringsAsFactors = FALSE,
  codelist_code = c("C66727","C66727","C66727","C66727","C66727",
                    "C66727","C66727","C66727","C66727","C66727"),
  term_code = c("C41331","C25250","C28554","C48226","C48227",
                "C48250","C142185","C49628","C49632","C49634"),
  term_value = c("ADVERSE EVENT","COMPLETED","DEATH","LACK OF EFFICACY",
                 "LOST TO FOLLOW-UP","PHYSICIAN DECISION","PROTOCOL VIOLATION",
                 "SCREEN FAILURE","STUDY TERMINATED BY SPONSOR",
                 "WITHDRAWAL BY SUBJECT"),
  collected_value = c("Adverse Event","Completed","Dead","Lack of Efficacy",
                      "Lost to Follow-Up","Physician Decision","Protocol Violation",
                      "Trial Screen Failure","Study Terminated by Sponsor",
                      "Withdrawal by Subject"),
  term_preferred_term = c("AE","Completed","Died",NA,NA,NA,"Violation",
                          "Failure to Meet Inclusion/Exclusion Criteria",NA,"Dropout"),
  term_synonyms = c("ADVERSE EVENT","COMPLETE","Death",NA,NA,NA,NA,NA,NA,
                    "Discontinued Participation")
)

# --- Extended manual mapping for terms not in CT ---
manual_map <- data.frame(
  stringsAsFactors = FALSE,
  collected_value = c(
    "Randomized",
    "Protocol Completed",
    "Screen Failure",
    "Sponsor Decision (Study or Patient Discontinued by the Sponsor)",
    "Unable to Contact Patient (Lost to Follow-Up)"
  ),
  term_value = c(
    "RANDOMIZED",
    "COMPLETED",
    "SCREEN FAILURE",
    "STUDY TERMINATED BY SPONSOR",
    "LOST TO FOLLOW-UP"
  )
)

# --- Helper function using CT and manual map ---
map_ct <- function(collected_val, ct_spec, manual) {
  result <- ct_spec$term_value[match(collected_val, ct_spec$collected_value)]
  idx_na <- is.na(result)
  if (any(idx_na)) {
    result[idx_na] <- ct_spec$term_value[match(collected_val[idx_na],
                                               ct_spec$term_synonyms)]
  }
  idx_na <- is.na(result)
  if (any(idx_na)) {
    result[idx_na] <- manual$term_value[match(collected_val[idx_na],
                                              manual$collected_value)]
  }
  return(result)
}

# --- Set up oak metadata ---
ds_raw_oak <- ds_raw %>%
  mutate(
    oak_id = row_number(),
    raw_source = "DS_RAW",
    patient_number = PATNUM
  )

# --- Load DM for RFSTDTC ---
dm <- pharmaversesdtm::dm %>%
  select(USUBJID, RFSTDTC)

# --- Build DS domain ---
ds <- ds_raw_oak %>%
  mutate(STUDYID = STUDY) %>%
  mutate(DOMAIN = "DS") %>%
  mutate(USUBJID = paste0(sub(".*?(\\d+)$", "\\1", STUDY), "-", PATNUM)) %>%
  mutate(DSTERM = case_when(
    !is.na(IT.DSTERM) & IT.DSTERM != "NA" ~ IT.DSTERM,
    !is.na(OTHERSP) & OTHERSP != "NA"     ~ OTHERSP,
    TRUE ~ NA_character_
  )) %>%
  mutate(DSDECOD = map_ct(IT.DSDECOD, study_ct, manual_map)) %>%
  mutate(DSCAT = case_when(
    FORM == "DISC1" ~ "DISPOSITION EVENT",
    TRUE ~ NA_character_
  )) %>%
  mutate(
    VISIT = INSTANCE,
    VISITNUM = case_when(
      INSTANCE == "Baseline" ~ 1,
      INSTANCE == "Week 4"   ~ 4,
      INSTANCE == "Week 8"   ~ 8,
      INSTANCE == "Week 12"  ~ 12,
      INSTANCE == "Week 16"  ~ 16,
      INSTANCE == "Week 20"  ~ 20,
      INSTANCE == "Week 24"  ~ 24,
      INSTANCE == "Week 26"  ~ 26,
      TRUE ~ NA_real_
    )
  ) %>%
  mutate(DSDTC = format(as.Date(DSDTCOL, format = "%m-%d-%Y"), "%Y-%m-%d")) %>%
  mutate(DSSTDTC = format(as.Date(IT.DSSTDAT, format = "%m-%d-%Y"), "%Y-%m-%d")) %>%
  filter(!is.na(DSTERM) | !is.na(DSDECOD)) %>%
  select(STUDYID, DOMAIN, USUBJID, DSTERM, DSDECOD, DSCAT,
         VISITNUM, VISIT, DSDTC, DSSTDTC) %>%
  left_join(dm, by = "USUBJID") %>%
  mutate(DSSTDY = as.numeric(as.Date(DSSTDTC) - as.Date(RFSTDTC)) + 1) %>%
  select(-RFSTDTC) %>%
  group_by(USUBJID) %>%
  mutate(DSSEQ = row_number()) %>%
  ungroup() %>%
  select(STUDYID, DOMAIN, USUBJID, DSSEQ, DSTERM, DSDECOD,
         DSCAT, VISITNUM, VISIT, DSDTC, DSSTDTC, DSSTDY)

# --- Review output ---
print(head(ds))
print(dim(ds))

# --- Check remaining NAs in DSDECOD ---
cat("\nDSDECOD distribution:\n")
print(table(ds$DSDECOD, useNA = "ifany"))

cat("\nRemaining missing DSDECOD:\n")
print(ds %>% filter(is.na(DSDECOD)) %>% count(DSTERM, sort = TRUE))

# --- Save output ---
saveRDS(ds, "ds_domain.rds")
write.csv(ds, "ds_domain.csv", row.names = FALSE)

message("DS domain created successfully with ", nrow(ds), " records.")

sink("ds_domain_log.txt")
cat("=== DS Domain Creation Log ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Input: pharmaverseraw::ds_raw\n")
cat("Output: ds_domain.rds / ds_domain.csv\n\n")
cat("Variables in output:\n")
print(names(ds))
cat("\nDimensions:", dim(ds), "\n")
cat("\nFirst 10 records:\n")
print(head(ds, 10))
cat("\nDSDECOD distribution:\n")
print(table(ds$DSDECOD, useNA = "ifany"))
cat("\nMissing DSDTC:", sum(is.na(ds$DSDTC)), "\n")
cat("Missing DSSTDTC:", sum(is.na(ds$DSSTDTC)), "\n")
cat("Missing DSSTDY:", sum(is.na(ds$DSSTDY)), "\n")
cat("\nProgram completed successfully with", nrow(ds), "records.\n")
sink()
message("Log file saved.")
