# This is the script for making figures for Chapter 1
# JAGS models and final output files

### Libraries:
library(rjags)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(stringr)
library(knitr)
library(gt)
library(gtExtras)

### Load data and files:

### Making a table for data products:

# read in table of products
data_products <- read.csv("CHAPTER_1/Figures/Tables/Data_Products_Table.csv")
# fix column names
colnames(data_products) <- str_replace_all(colnames(data_products), "\\.", " ")
colnames(data_products)[6] <- "Temporal Resolution (for analyses)"

table <- data_products %>% 
  gt() %>%
  tab_options(table.font.names = 'Times New Roman',
              column_labels.font.weight = 'bold',
              heading.title.font.size = 14,
              heading.subtitle.font.size = 14,
              table.font.color = 'black',
              table.font.size = 12) %>%
  tab_header(title = "Table 1: Covariate Data for Univariate and Multivariate Model Selection") %>%
  opt_align_table_header("left") %>%
  #gt_theme_pff() %>%
  cols_width(everything() ~ px(90)) %>%
  opt_table_font(size = 12)

gtsave(table, "2025_05_11_T1_data_products_table.png")


# # png("2024_09_Ch1_data_products_table.png", 
# #     width = 8, height = 5, units ="in", res = 800)
# # prep table
# table <- kbl(data_products, align = "l") %>%
#   kable_classic(full_width = F, html_font = "Cambria", font_size = 10) %>%
#   row_spec(0, bold = TRUE)

# save_kable(table, "2024_09_Ch1_data_products_table.html")
# webshot::webshot("2024_09_Ch1_data_products_table.html",
#                  file = "2024_09_Ch1_data_products_table.png")
# dev.off



### HEATMAP ### ------
mnames <- c(
  "2024-10-02_modelrun_joint_a_21_b_14",
  "2024-10-04_modelrun_joint_a_18_b_14",
  "2024-10-04_modelrun_joint_a_20_b_17",
  "2024-10-04_modelrun_joint_a_20_b_15",
  "2024-10-05_modelrun_joint_a_20_b_14",
  "2024-09-25_modelrun_alpha_a_10_b_1",
  "2024-09-26_modelrun_alpha_a_18_b_1",
  "2024-09-27_modelrun_alpha_a_20_b_1",
  "2024-09-27_modelrun_alpha_a_21_b_1",
  "2024-09-11_modelrun_beta_a_1_b_13",
  "2024-09-11_modelrun_beta_a_1_b_14",
  "2024-09-12_modelrun_beta_a_1_b_15",
  "2024-09-13_modelrun_beta_a_1_b_16",
  "2024-09-13_modelrun_beta_a_1_b_17")

# for making supplement:
file_names <- list.files("CHAPTER_1/all_JAGS_models/")
mnames <- gsub("_data.RData", "", file_names)

columns <- c("R", "alpha0", "alpha[1]", "alpha[2]", "alpha[3]", "alpha[4]",
             "alpha[5]", "alpha[6]", "alpha[7]", "alpha[8]", "beta0", 
             "beta[1]", "beta[2]", "beta[3]","beta[4]",
             "beta[5]", "beta[6]", "beta[7]", "beta[8]", 
             "beta[9]", "pa0", "DIC")
rows <- gsub(".*n_", "" ,mnames)

#making empty matrices to fill
# for means:
best_models <- matrix(data = NA, nrow = length(mnames), ncol = length(columns))
# for SDs:
SDs <- best_models
a_vars <- matrix(data = NA, nrow = length(mnames), ncol = length(grep("^a", columns))-1)
b_vars <- matrix(data = NA, nrow = length(mnames), ncol = length(grep("^b", columns))-1)
colnames(best_models) <- columns
rownames(best_models) <- rows
colnames(SDs) <- columns
rownames(SDs) <- rows
colnames(a_vars) <- columns[grep("^a", columns)][-1]
colnames(b_vars) <- columns[grep("^b", columns)][-1]

for (i in 1:length(mnames)){
  # load model
  # model <- (paste0("CHAPTER_1/2024_09_JAGS_models/best_models/",
  #                  mnames[i],"_data.RData"))
  model <- (paste0("CHAPTER_1/all_JAGS_models/",
                   mnames[i], "_data.RData"))
  load(model)
  # remove burn in
  jpout <- model_info$jpout
  burnin <- 100000
  jburn <- window(jpout, start = burnin)
  # extract means
  summ <- summary(jburn)
  # fill in DIC
  best_models[i,"DIC"] <- model_info$dic[[2]]
  # put names in data frame
  for (j in names(summ$statistics[,"Mean"])){
    best_models[i, j] <- summ$statistics[j,"Mean"]
  }
  for (j in names(summ$statistics[,"SD"])){
    SDs[i, j] <- summ$statistics[j,"SD"]
  }
  # get variable names
  if (model_info$metadata$covs == 1){
    for(k in 1:ncol(dpls[[model_info$metadata$modelrun_a]])){
      a_vars[i, k] <- colnames(dpls[[model_info$metadata$modelrun_a]])[k]
    }
    # for (k in 1:length(names(model_info$metadata$data$a))){
    #   a_vars[i, k] <- names(model_info$metadata$data$a)[k]
    # }
  } else if (model_info$metadata$covs == 2) {
    for(k in 1:ncol(dmls[[model_info$metadata$modelrun_b]])){
      b_vars[i, k] <- colnames(dmls[[model_info$metadata$modelrun_b]])[k]
    }
    # for (k in 1:length(names(model_info$metadata$data$b))){
    #   b_vars[i, k] <- names(model_info$metadata$data$b)[k]
    # }
  } else if (model_info$metadata$covs == 3){
    for(k in 1:ncol(dpls[[model_info$metadata$modelrun_a]])){
      a_vars[i, k] <- colnames(dpls[[model_info$metadata$modelrun_a]])[k]
    }
    # for (k in 1:length(names(model_info$metadata$data$a))){
    #   a_vars[i, k] <- names(model_info$metadata$data$a)[k]
    # }
    for(k in 1:ncol(dmls[[model_info$metadata$modelrun_b]])){
      b_vars[i, k] <- colnames(dmls[[model_info$metadata$modelrun_b]])[k]
    # for (k in 1:length(names(model_info$metadata$data$b))){
    #   b_vars[i, k] <- names(model_info$metadata$data$b)[k]
    }
  }
}
rm(jburn, jpout, model_info, summ, burnin, i, j, k, model)

# round results to 3 digits
best_models <- as.data.frame(best_models) %>% mutate_if(is.numeric, round, digits = 3)
SDs <- as.data.frame(SDs) %>% mutate_if(is.numeric, round, digits = 3)
# calculate delDIC
best_models$delDIC <- min(best_models$DIC) - best_models$DIC
# sort
#best_models <- best_models[order(best_models$delDIC, decreasing = T),]
# make new alpha and beta intercept columns
best_models[which(!is.na(best_models$alpha0)),'alpha[1]'] <- c(best_models$alpha0[which(is.na(best_models$`alpha[1]`))])
best_models[which(!is.na(best_models$beta0)), 'beta[1]'] <- c(best_models$beta0[which(is.na(best_models$`beta[1]`))])

# matching models to correct variable parameters
# put covariate names in one table
#a_vars[which(a_vars[,1] == "int"),1] <- "a_int"
#b_vars[which(b_vars[,1] == "int"),1] <- "b_int"
a_vars[c(which(!is.na(a_vars)))] <- paste0("a_", a_vars[c(which(!is.na(a_vars)))])
b_vars[c(which(!is.na(b_vars)))] <- paste0("b_", b_vars[c(which(!is.na(b_vars)))])
vars <- as.data.frame(cbind(a_vars, b_vars))
rownames(vars) <- rows
# make matrix to fill
model_params <- matrix(data = NA, nrow = length(mnames), ncol = ncol(vars))
model_SDs <- matrix(data = NA, nrow = length(mnames), ncol = ncol(vars))
colnames(model_params) <- c(as.character(vars[which(!is.na(vars["alpha[8]"]))[1], grep("^a", names(vars))]),
                            as.character(vars[which(!is.na(vars["beta[9]"]))[1], grep("^b", names(vars))]))
colnames(model_SDs) <- c(as.character(vars[which(!is.na(vars["alpha[8]"]))[1], grep("^a", names(vars))]),
                            as.character(vars[which(!is.na(vars["beta[9]"]))[1], grep("^b", names(vars))]))

# fill in matrix values
for (i in 1:nrow(model_params)){
  for(j in which(!is.na(vars[i,]))){
    var <- vars[i,j]
    name <- names(vars)[j]
    model_params[i,var] <- best_models[i,name]
    model_SDs[i,var] <- SDs[i,name]
  }
}

# fill in intercepts
model_params[which(is.na(model_params[,"a_int"])), "a_int"] <- best_models$alpha0[which(!is.na(best_models$alpha0))]
model_params[which(is.na(model_params[,"b_int"])), "b_int"] <- best_models$beta0[which(!is.na(best_models$beta0))]
model_SDs[which(is.na(model_SDs[,"a_int"])), "a_int"] <- SDs$alpha0[which(!is.na(SDs$alpha0))]
model_SDs[which(is.na(model_SDs[,"b_int"])), "b_int"] <- SDs$beta0[which(!is.na(SDs$beta0))]

# make final table
JAGS_params <- as.data.frame(cbind(model_params, best_models$R, best_models$pa0, best_models$DIC, best_models$delDIC))#, 
JAGS_SDs <- as.data.frame(cbind(model_SDs, SDs$R, SDs$pa0))
#alpha0 = best_models$alpha0,
#beta0 = best_models$beta0))
# rename columns
JAGS_params <- JAGS_params %>%
  rename('Alpha Intercept' = a_int) %>%
  rename('Beta Intercept' = b_int)
colnames(JAGS_params) <- c(gsub("a_", "", colnames(JAGS_params)))
colnames(JAGS_params) <- c(gsub("b_", "", colnames(JAGS_params)))
colnames(JAGS_params) <- c(gsub("_", replacement = " ", colnames(JAGS_params)))

# make values for NAs to be able to print numbers in table
JAGS_params[is.na(JAGS_params)] <- 1000
# make this dark grey:
# all NAs/1000s
JAGS_params.col <- JAGS_params > 999
# set all values that satisfy the condition to color
JAGS_params.col <- gsub("TRUE", "grey45", JAGS_params.col)
# set all values that do not satisfy the condition to "black"
JAGS_params.col <- gsub("FALSE", "black", JAGS_params.col)
# convert to matrix
JAGS_params.col <- matrix(JAGS_params.col, ncol = ncol(JAGS_params))
# add row and column names to make broken-up figures:
#rownames(JAGS_params.col) <- rownames(JAGS_params)
colnames(JAGS_params.col) <- colnames(JAGS_params)

# make heatmap:
png("2024_10_07_JAGS_params_heatmap.png",
    width = 10, height = 8, units = "in", res = 300)
superheat(as.matrix(t(JAGS_params)), 
          scale = FALSE, # the scale normalizes to mean 0 SD 1
          left.label.text.size = 4,
          bottom.label.text.size = 4,
          bottom.label.text.angle = 90,
          heat.pal = c("dodgerblue", "skyblue1", "white", "pink1","firebrick1"),
          heat.pal.values = c(0, 0.45, 0.5, 0.55, 1),
          heat.lim = c(-7,7),
          legend.num.ticks = 8,
          column.title = "Models",
          row.title = "Parameters",
          column.title.size = 4,
          row.title.size = 4,
          X.text = as.matrix(t(JAGS_params)),
          X.text.col = t(JAGS_params.col),
          heat.na.col = "grey45",
          X.text.size = 3,
          extreme.values.na = TRUE,
          left.label.size = 0.3,
          bottom.label.size = 0.15,
          legend.vspace = 0.01)
dev.off()
