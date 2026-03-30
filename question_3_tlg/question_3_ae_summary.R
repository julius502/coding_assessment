# =============================================================================
# Question 3: TLG - Adverse Events Reporting
# Input: pharmaverseadam::adae, pharmaverseadam::adsl
# Output: AE summary table, AE severity plot, Top 10 AE plot
# =============================================================================

library(gtsummary)
library(ggplot2)
library(pharmaverseadam)
library(dplyr)

# --- Load data ---
adae <- pharmaverseadam::adae
adsl <- pharmaverseadam::adsl

# --- Filter to TEAEs only ---
teae <- adae %>% filter(TRTEMFL == "Y")

# --- Get subject counts per treatment ---
n_arm <- adsl %>%
  filter(ACTARM != "Screen Failure") %>%
  count(ACTARM, name = "N")

cat("Subject counts per arm:\n")
print(n_arm)
cat("TEAE records:", nrow(teae), "\n")

# =============================================================================
# OUTPUT 1: Summary Table using gtsummary
# =============================================================================

# Prepare data - one record per subject per SOC/AETERM
ae_subj <- teae %>%
  select(USUBJID, ACTARM, AESOC, AETERM) %>%
  distinct()

# Join to get total N per arm for percentages
ae_tbl <- ae_subj %>%
  left_join(adsl %>% select(USUBJID, ACTARM) %>% distinct(), by = "USUBJID")

# Count subjects with each AE term by treatment
ae_summary <- teae %>%
  distinct(USUBJID, ACTARM, AESOC, AETERM) %>%
  group_by(ACTARM, AESOC, AETERM) %>%
  summarise(n = n_distinct(USUBJID), .groups = "drop") %>%
  left_join(n_arm, by = "ACTARM") %>%
  mutate(pct = round(100 * n / N, 1),
         cell = paste0(n, " (", pct, "%)")) %>%
  select(AESOC, AETERM, ACTARM, cell) %>%
  tidyr::pivot_wider(names_from = ACTARM, values_from = cell, values_fill = "0 (0.0%)") %>%
  arrange(AESOC, AETERM)

# Add total column
ae_total <- teae %>%
  distinct(USUBJID, AESOC, AETERM) %>%
  group_by(AESOC, AETERM) %>%
  summarise(n_total = n_distinct(USUBJID), .groups = "drop") %>%
  mutate(
    N_total = n_distinct(adsl$USUBJID[adsl$ACTARM != "Screen Failure"]),
    pct_total = round(100 * n_total / N_total, 1),
    Total = paste0(n_total, " (", pct_total, "%)")
  ) %>%
  select(AESOC, AETERM, Total)

ae_summary <- ae_summary %>%
  left_join(ae_total, by = c("AESOC", "AETERM")) %>%
  arrange(desc(as.numeric(sub(" .*", "", Total))))

# Save as HTML using gt
library(gt)
ae_gt <- ae_summary %>%
  gt() %>%
  tab_header(
    title = "Treatment Emergent Adverse Events",
    subtitle = "Subjects with at least one TEAE by System Organ Class and Preferred Term"
  ) %>%
  cols_label(
    AESOC = "System Organ Class",
    AETERM = "Preferred Term"
  )

gtsave(ae_gt, "ae_summary_table.html")
cat("AE summary table saved.\n")

# =============================================================================
# OUTPUT 2: Plot 1 - AE Severity Distribution by Treatment
# =============================================================================

sev_data <- teae %>%
  filter(!is.na(AESEV)) %>%
  group_by(ACTARM, AESEV) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(AESEV = factor(AESEV, levels = c("MILD", "MODERATE", "SEVERE")))

p1 <- ggplot(sev_data, aes(x = ACTARM, y = count, fill = AESEV)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(
    values = c("MILD" = "#F08080", "MODERATE" = "#228B22", "SEVERE" = "#6495ED"),
    name = "Severity/Intensity"
  ) +
  labs(
    title = "AE severity distribution by treatment",
    x = "Treatment Arm",
    y = "Count of AEs"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 15, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

ggsave("ae_severity_plot.png", plot = p1, width = 8, height = 6, dpi = 150)
cat("AE severity plot saved.\n")

# =============================================================================
# OUTPUT 3: Plot 2 - Top 10 Most Frequent AEs with 95% CI
# =============================================================================

n_total <- n_distinct(adsl$USUBJID[adsl$ACTARM != "Screen Failure"])

top10 <- teae %>%
  distinct(USUBJID, AETERM) %>%
  group_by(AETERM) %>%
  summarise(n = n_distinct(USUBJID), .groups = "drop") %>%
  mutate(
    pct = n / n_total,
    # Clopper-Pearson 95% CI
    ci_low  = qbeta(0.025, n, n_total - n + 1),
    ci_high = qbeta(0.975, n + 1, n_total - n)
  ) %>%
  arrange(desc(n)) %>%
  slice(1:10) %>%
  mutate(AETERM = reorder(AETERM, pct))

p2 <- ggplot(top10, aes(x = pct, y = AETERM)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.3) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Top 10 Most Frequent Adverse Events",
    subtitle = paste0("n = ", n_total, " subjects; 95% Clopper-Pearson CIs"),
    x = "Percentage of Patients (%)",
    y = NULL
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("ae_top10_plot.png", plot = p2, width = 8, height = 6, dpi = 150)
cat("Top 10 AE plot saved.\n")

message("Question 3 complete - all outputs saved.")

# Fix deprecated geom_errorbarh
p2 <- ggplot(top10, aes(x = pct, y = AETERM)) +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), 
                width = 0.3, orientation = "y") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Top 10 Most Frequent Adverse Events",
    subtitle = paste0("n = ", n_total, " subjects; 95% Clopper-Pearson CIs"),
    x = "Percentage of Patients (%)",
    y = NULL
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("ae_top10_plot.png", plot = p2, width = 8, height = 6, dpi = 150)

# Save log
sink("ae_tlg_log.txt")
cat("=== AE TLG Creation Log ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Input: pharmaverseadam::adae, pharmaverseadam::adsl\n\n")
cat("TEAE records:", nrow(teae), "\n")
cat("Subject counts per arm:\n")
print(n_arm)
cat("\nTop 10 AEs:\n")
print(top10 %>% select(AETERM, n, pct) %>% arrange(desc(n)))
cat("\nOutputs saved:\n")
cat("- ae_summary_table.html\n")
cat("- ae_severity_plot.png\n")
cat("- ae_top10_plot.png\n")
cat("\nProgram completed successfully.\n")
sink()
message("Log saved.")
