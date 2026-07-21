##############################################################################
# 16. 全样本童年 SES 提取（独立于财富记录）
# 从各队列的协调文件中直接提取父母教育变量，不以财富/收入为前提
##############################################################################

library(dplyr)
library(haven)

source(if (file.exists("analysis/00_config.R")) "analysis/00_config.R" else "00_config.R")
out_dir <- PROJECT_DIR
df_main <- readRDS(file.path(out_dir, "clean_longitudinal_data.rds"))

# 全样本人员清单（只要有过至少一次认知评估）
all_ids <- df_main %>%
  filter(!is.na(cogtot)) %>%
  distinct(id, study, enroll_year)

cat("全样本人数:", nrow(all_ids), "\n")

cat("================================================================\n")
cat("  HRS: rafeduc, rameduc (连续年限)\n")
cat("================================================================\n")
hrs_raw <- read_dta(HRS_FILE(),
  col_select = c("hhidpn", "rafeduc", "rameduc"))
hrs_parent <- hrs_raw %>%
  mutate(id = paste0("HRS_", hhidpn)) %>%
  select(id, child_dad_edu = rafeduc, child_mom_edu = rameduc)
cat("  HRS 提取:", nrow(hrs_parent), "行\n")

cat("================================================================\n")
cat("  ELSA: radadeduage, ramomeduage (1-7 类别码)\n")
cat("================================================================\n")
elsa_raw <- read_dta(ELSA_FILE(),
  col_select = c("idauniq", "radadeduage", "ramomeduage"))
elsa_parent <- elsa_raw %>%
  mutate(id = paste0("ELSA_", idauniq)) %>%
  select(id, child_dad_edu = radadeduage, child_mom_edu = ramomeduage)
cat("  ELSA 提取:", nrow(elsa_parent), "行\n")

cat("================================================================\n")
cat("  SHARE: radadeducl, ramomeducl (0-3 类别码)\n")
cat("================================================================\n")
share_raw <- read_dta(SHARE_FILE(),
  col_select = c("mergeid", "radadeducl", "ramomeducl"))
share_parent <- share_raw %>%
  mutate(id = paste0("SHARE_", mergeid)) %>%
  select(id, child_dad_edu = radadeducl, child_mom_edu = ramomeducl)
cat("  SHARE 提取:", nrow(share_parent), "行\n")

cat("================================================================\n")
cat("  CHARLS: radadeducl, ramomeducl (0-3 类别码)\n")
cat("================================================================\n")
charls_raw <- read_dta(CHARLS_FILE(),
  col_select = c("ID", "radadeducl", "ramomeducl"))
charls_parent <- charls_raw %>%
  mutate(id = paste0("CHARLS_", ID)) %>%
  select(id, child_dad_edu = radadeducl, child_mom_edu = ramomeducl)
cat("  CHARLS 提取:", nrow(charls_parent), "行\n")

cat("================================================================\n")
cat("  合并、去标签、与全样本匹配\n")
cat("================================================================\n")

all_parent <- bind_rows(
  haven::zap_labels(hrs_parent),
  haven::zap_labels(elsa_parent),
  haven::zap_labels(share_parent),
  haven::zap_labels(charls_parent)
)

# 只保留分析人群中的 ID
fullsample_parent <- all_ids %>%
  left_join(all_parent, by = "id")

cat("全样本匹配:", nrow(fullsample_parent), "人\n")
cat("有父母教育信息:", sum(!is.na(fullsample_parent$child_dad_edu) |
                          !is.na(fullsample_parent$child_mom_edu)), "人\n")
cat("双亲均缺失:", sum(is.na(fullsample_parent$child_dad_edu) &
                       is.na(fullsample_parent$child_mom_edu)), "人\n")

# 按队列统计
cat("\n--- 按队列统计 ---\n")
fullsample_parent %>%
  group_by(study) %>%
  summarize(
    n = n(),
    has_dad = sum(!is.na(child_dad_edu)),
    has_mom = sum(!is.na(child_mom_edu)),
    has_either = sum(!is.na(child_dad_edu) | !is.na(child_mom_edu)),
    both_missing = sum(is.na(child_dad_edu) & is.na(child_mom_edu)),
    .groups = "drop"
  ) %>%
  mutate(pct_covered = round(100 * has_either / n, 1)) %>%
  print()

saveRDS(fullsample_parent, file.path(out_dir, "fullsample_childhood_ses.rds"))
cat("\n全样本童年 SES 已保存至 fullsample_childhood_ses.rds\n")
