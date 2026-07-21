##############################################################################
# 12. 逐波逆概率加权 (IPCW) 应对竞争风险 (12_build_ipw.R)
# 修正版本 (Phase 3.5 Final): 严格控制个体基线、死亡风险剔除与条件概率
##############################################################################

library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)

source(if (file.exists("analysis/00_config.R")) "analysis/00_config.R" else "00_config.R")
out_dir <- PROJECT_DIR
df <- readRDS(file.path(out_dir, "data_with_covariates.rds"))
cov_long <- readRDS(file.path(out_dir, "all_covariates_long.rds"))

cat("================================================================\n")
cat("  Step 1: 构建严格的纵向状态转移矩阵与风险集\n")
cat("================================================================\n")

# 只针对纵向样本中的个体
ids <- unique(df$id)
cov_sub <- cov_long %>% filter(id %in% ids) %>% arrange(id, wave)

# 确定每个人的分析基线波次
base_wave_df <- df %>% filter(!is.na(cogtot)) %>% 
  group_by(id) %>% 
  summarize(analysis_base_wave = min(wave)) %>%
  ungroup()

cov_sub <- cov_sub %>% left_join(base_wave_df, by="id")

# 加入 df_cog 提取必要协变量
df_cog <- df %>% select(id, wave, cogtot, age_base_c, enroll_year_c, female, edu3_f, study, 
                        base_married, base_diabe, base_hibpe, base_hearte, base_stroke, time_in_study)

cov_sub <- cov_sub %>% left_join(df_cog, by=c("id", "wave"))

# 填充不随时间改变的基线协变量
cov_sub <- cov_sub %>%
  group_by(id) %>%
  fill(female, edu3_f, study, analysis_base_wave, age_base_c, enroll_year_c,
       base_married, base_diabe, base_hibpe, base_hearte, base_stroke, .direction = "downup") %>%
  ungroup()

# 处理慢病缺失值 (NA -> 999 因子)
vars_to_missing <- c("base_diabe", "base_hibpe", "base_hearte", "base_stroke", "base_married")
for(v in vars_to_missing) {
  cov_sub[[v]][is.na(cov_sub[[v]])] <- 999
  cov_sub[[v]] <- as.factor(cov_sub[[v]])
}

# 状态定义与滞后变量 (lag)
cov_sub <- cov_sub %>%
  mutate(iwstat_num = as.numeric(haven::zap_labels(iwstat))) %>%
  mutate(status = case_when(
    iwstat_num %in% c(5, 6) ~ "3_Dead",
    iwstat_num == 0 ~ "4_Inapplicable",
    iwstat_num == 1 ~ "1_Observed",
    TRUE ~ "2_Dropped" # iwstat_num %in% c(4, 7, 9)
  )) %>%
  mutate(is_observed = ifelse(status == "1_Observed" & !is.na(cogtot), 1, 0)) %>%
  group_by(id) %>%
  arrange(wave) %>%
  mutate(
    lag_status = lag(status),
    lag_cogtot = lag(cogtot),
    lag_is_observed = lag(is_observed),
    lag_time = lag(time_in_study)
  ) %>%
  ungroup()

cat("\n================================================================\n")
cat("  Step 2: 拟合条件稳定化 IPCW 模型 (排除死亡竞争风险)\n")
cat("================================================================\n")

# 构建条件风险集：
# 1. 仅考虑每个人 analysis_base_wave 之后的波次
# 2. 如果当前 wave 死亡 (3_Dead)，将其剔除出需要被“重新加权”的风险集
# 3. 如果上一波已经死亡，自然不再进入风险集
cov_model_data <- cov_sub %>% 
  filter(wave > analysis_base_wave) %>%
  filter(status != "3_Dead" & status != "4_Inapplicable") %>%
  filter(is.na(lag_status) | lag_status != "3_Dead")

# 对于 lag_cogtot 缺失的情况（比如上一波是 2_Dropped），填补总体均值以防样本流失，
# 并加入 lag_is_observed 作为一个指示变量
mean_cog <- mean(cov_model_data$lag_cogtot, na.rm=TRUE)
cov_model_data$lag_cogtot_imp <- ifelse(is.na(cov_model_data$lag_cogtot), mean_cog, cov_model_data$lag_cogtot)
cov_model_data$lag_is_observed[is.na(cov_model_data$lag_is_observed)] <- 1
cov_model_data$wave_f <- as.factor(cov_model_data$wave)

# 分子模型 (仅含主要暴露变量和时间因子)
cat("拟合分子模型 (Numerator)...\n")
m_num <- glm(is_observed ~ wave_f * study + female + edu3_f, 
             family = binomial(), data = cov_model_data)

# 分母模型 (加入年龄、基线健康状况、之前是否观测到、之前的认知水平)
cat("拟合分母模型 (Denominator)...\n")
m_den <- glm(is_observed ~ wave_f * study + female + edu3_f + 
             age_base_c + base_diabe + base_hibpe + base_hearte + base_stroke + base_married +
             lag_is_observed + lag_cogtot_imp, 
             family = binomial(), data = cov_model_data)

cov_model_data$p_num <- predict(m_num, type = "response", newdata = cov_model_data)
cov_model_data$p_den <- predict(m_den, type = "response", newdata = cov_model_data)

# 稳定化条件权重，处理间歇性缺失
# 如果 is_observed == 1，贡献 p_num / p_den
# 如果 is_observed == 0，贡献 (1 - p_num) / (1 - p_den)
cov_model_data <- cov_model_data %>%
  mutate(
    w_t = ifelse(is_observed == 1, p_num / p_den, (1 - p_num) / (1 - p_den))
  )

# 累乘得到最终权重 (1% - 99% Truncation)
cov_model_data <- cov_model_data %>%
  mutate(w_t_trunc = pmin(pmax(w_t, quantile(w_t, 0.01, na.rm=T)), quantile(w_t, 0.99, na.rm=T))) %>%
  group_by(id) %>%
  arrange(wave) %>%
  mutate(ipcw = cumprod(w_t_trunc)) %>%
  ungroup()

ipcw_weights <- cov_model_data %>% select(id, wave, ipcw)

# 合并回最终分析集
df_ipcw <- df %>% left_join(ipcw_weights, by=c("id", "wave"))
# 基础波次和观测到的认知因为没有被剔除，给基线赋1
df_ipcw$ipcw[is.na(df_ipcw$ipcw) & !is.na(df_ipcw$cogtot) & df_ipcw$wave == df_ipcw$baseline_wave] <- 1 
# 为了防止有些虽然观测到但是前面有inapplicable等导致的NA
df_ipcw <- df_ipcw %>% 
  group_by(id) %>% arrange(wave) %>% fill(ipcw, .direction="down") %>% ungroup()
df_ipcw$ipcw[is.na(df_ipcw$ipcw) & !is.na(df_ipcw$cogtot)] <- 1

# Truncate final weight again globally to prevent highly influential outliers
q99_final <- quantile(df_ipcw$ipcw, 0.99, na.rm=TRUE)
df_ipcw$ipcw <- pmin(df_ipcw$ipcw, q99_final)

cat("\n最终参与建模的观测的 IPCW 分布:\n")
summary_ipcw <- summary(df_ipcw$ipcw[!is.na(df_ipcw$cogtot)])
print(summary_ipcw)

cat("\n--- 各队列 x 教育组权重诊断 ---\n")
df_ipcw %>% 
  filter(!is.na(cogtot)) %>%
  group_by(study, edu3_f) %>%
  summarize(
    N = n(),
    mean_W = mean(ipcw, na.rm=T),
    max_W = max(ipcw, na.rm=T),
    prop_truncated = mean(ipcw == q99_final, na.rm=T) * 100,
    ESS = (sum(ipcw, na.rm=T))^2 / sum(ipcw^2, na.rm=T),
    ESS_ratio = ESS / N,
    .groups = "drop"
  ) %>%
  print(n=Inf)

cat("\n================================================================\n")
cat("  Step 3: 保存带 IPCW 的数据集\n")
cat("================================================================\n")

saveRDS(df_ipcw, file.path(out_dir, "data_with_ipcw.rds"))

cat("\nIPCW 数据集保存完成 (data_with_ipcw.rds)。\n")
