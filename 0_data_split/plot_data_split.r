library(ggupset)
library(ggplot2)
library(tidyverse, warn.conflicts = FALSE)
library(argparse)
library(kableExtra)
library(gridExtra)
library(grid)
library(nord)
library(ggrepel)
library(RColorBrewer)
library(patchwork)
library(plotly)
library(htmlwidgets)
# library(gghalves) and install it if not exists
if (!requireNamespace("gghalves", quietly = TRUE)) {
        install.packages("gghalves")
}

library(gghalves)

source("/home/nrlab/wang04/ulyses/models/cnn_xgboost_dataset_hyperparams.R")

################################################################################

# reain the test data from final_test_hyper
load(final_test_hyper)
all_filtered <- readRDS(xgboost_input_data_file_hyper) %>%
        filter(assay == "sd") %>%
        select(-rowname, -colname, -assay, -value) %>%
        distinct()
final_test_bamid <- final_test_primary_vec


ans <- all_filtered %>%
        mutate(split = ifelse(primary %in% final_test_bamid, "Test", "Train"))

ans2 <- ans %>%
        mutate(`all` = TRUE) %>%
        mutate(`[0, 0.03]` = case_when(
                ichorcna_tf_strat %in% c("Healthy", "[0, 0.03]") ~ TRUE,
                TRUE ~ FALSE
        )) %>%
        mutate(`(0.03, 0.1]` = case_when(
                ichorcna_tf_strat %in% c("Healthy", "(0.03, 0.1]") ~ TRUE,
                TRUE ~ FALSE
        )) %>%
        mutate(`(0.1, 1]` = case_when(
                ichorcna_tf_strat %in% c("Healthy", "(0.1, 1]") ~ TRUE,
                TRUE ~ FALSE
        ))

# pivot longer, gather `all`, `[0, 0.03]`, `(0.03, 0.1]`, `(0.1, 1]` to a column called "strat", values to a column called "strat_bool"
ans3 <- ans2 %>%
        pivot_longer(
                cols = c(`all`, `[0, 0.03]`, `(0.03, 0.1]`, `(0.1, 1]`),
                names_to = "strat",
                values_to = "strat_bool"
        )


ans4 <- ans3 %>%
        group_by(cohort, split, strat, strat_bool) %>%
        # summarise how many rows with strat_bool == TRUE
        summarise(count = n()) %>%
        ungroup() %>%
        filter(strat_bool == TRUE)

# set the "Healthy" to the first factor level in cohort
ans4$cohort <- factor(ans4$cohort)

ans4$cohort <- relevel(ans4$cohort, "Healthy")

# reverse the factor level of cohort
ans4$cohort <- fct_rev(ans4$cohort)

# set the factor levels of strat to `[0, 0.03]`, `(0.03, 0.1]`, `(0.1, 1]` and `all`
ans4$strat <- factor(ans4$strat, levels = c("[0, 0.03]", "(0.03, 0.1]", "(0.1, 1]", "all"))

# reverse the factor level of strat
ans4$strat <- fct_rev(ans4$strat)

# set the factor levels of split to "Train" and "Test"
ans4$split <- factor(ans4$split, levels = c("Train", "Test"))


# set color to distinct color, the number of colors equals to the number of cohorts in ans4
n_colors <- length(unique(ans4$cohort))
# use distinct colors for categorical data

qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# set seed to 0
set.seed(1)
manual_cols <- sample(col_vector, n_colors)
names(manual_cols) <- unique(ans4$cohort)
# change the color of Healthy to the "black"
manual_cols["Healthy"] <- "darkgrey"


ans5 <- ans4 %>%
        group_by(split, strat) %>%
        mutate(total_count = sum(count)) %>%
        # calculate thte sum of "Healthy" cohort
        mutate(total_control = sum(count[cohort == "Healthy"])) %>%
        # calculate the sum of non-healthy
        mutate(total_case = total_count - total_control) %>%
        # drop group_by
        ungroup()
ans6 <- ans5 %>%
        mutate(y = paste0(strat, "\n", "Ctrl:", total_control, "; ", "Case:", total_case))
# train
ans6_train <- ans6 %>%
        filter(split == "Train")

ans6_train_all <- ans6_train %>%
        filter(strat == "all") %>%
        pull(y) %>%
        unique()
# make ans6_train_0_0.03
ans6_train_0_0.03 <- ans6_train %>%
        filter(strat == "[0, 0.03]") %>%
        pull(y) %>%
        unique()

# make ans6_train_0.03_0.1
ans6_train_0.03_0.1 <- ans6_train %>%
        filter(strat == "(0.03, 0.1]") %>%
        pull(y) %>%
        unique()

ans6_train_0.1_1 <- ans6_train %>%
        filter(strat == "(0.1, 1]") %>%
        pull(y) %>%
        unique()

# set levels of y to  ans6_train_0_0.03, ans6_train_0.03_0.1, ans6_train_0.1_1, ans6_train_all
ans6_train$y <- factor(ans6_train$y, levels = c(ans6_train_0_0.03, ans6_train_0.03_0.1, ans6_train_0.1_1, ans6_train_all))
# reverse the factor level of y
ans6_train$y <- fct_rev(ans6_train$y)

p_train <- ans6_train %>%
        ggplot(aes(y = y, x = count, fill = cohort)) +
        geom_bar(position = "stack", stat = "identity", width = 0.5) +
        # facet_wrap(~split, ncol = 1) +
        theme_classic() +
        # put legend at the bottom
        theme(legend.position = "right") +
        labs(
                x = "Count (n)",
                y = "ichorCNA TF strata",
                fill = "Strat all"
        ) +
        # remove y axis
        theme(axis.title.y = element_blank()) +
        # manually set the color of fill
        scale_fill_manual(values = manual_cols) +
        # x axis expand to 0
        scale_x_continuous(expand = c(0, 0)) +
        # y axis expand to 0
        scale_y_discrete(expand = c(0, 0)) +
        coord_cartesian(xlim = c(0, 1300)) +
        # remove strip border
        theme(
                # strip background color to lightgrey
                strip.background = element_blank()
        ) +
        # change axis text size to 5
        theme(axis.text = element_text(size = 7)) +
        # axis title size to 7
        theme(axis.title = element_text(size = 7)) +
        # annotate the count on the bar using ggrepel
        geom_text_repel(aes(label = count), position = position_stack(vjust = 0), size = 2, max.overlaps = Inf)
# save as pdf

# test

ans6_test <- ans6 %>%
        filter(split == "Test")

ans6_test_all <- ans6_test %>%
        filter(strat == "all") %>%
        pull(y) %>%
        unique()
# make ans6_train_0_0.03
ans6_test_0_0.03 <- ans6_test %>%
        filter(strat == "[0, 0.03]") %>%
        pull(y) %>%
        unique()

# make ans6_train_0.03_0.1
ans6_test_0.03_0.1 <- ans6_test %>%
        filter(strat == "(0.03, 0.1]") %>%
        pull(y) %>%
        unique()

ans6_test_0.1_1 <- ans6_test %>%
        filter(strat == "(0.1, 1]") %>%
        pull(y) %>%
        unique()

# set levels of y to  ans6_train_0_0.03, ans6_train_0.03_0.1, ans6_train_0.1_1, ans6_train_all
ans6_test$y <- factor(ans6_test$y, levels = c(ans6_test_0_0.03, ans6_test_0.03_0.1, ans6_test_0.1_1, ans6_test_all))

# reverse the factor level of y
ans6_test$y <- fct_rev(ans6_test$y)

# plot stacked bar, y is `all`, stacked by cohort
p_test <- ggplot(ans6_test, aes(y = y, x = count, fill = cohort)) +
        geom_bar(position = "stack", stat = "identity", width = 0.5) +
        # facet_wrap(~split, ncol = 1) +
        theme_classic() +
        # put legend at the bottom
        theme(legend.position = "right") +
        labs(
                x = "Count (n)",
                y = "ichorCNA TF strata",
                fill = "Strat all"
        ) +
        # remove y axis
        theme(axis.title.y = element_blank()) +
        # manually set the color of fill
        # x axis expand to 0
        scale_x_continuous(expand = c(0, 0)) +
        # y axis expand to 0
        scale_y_discrete(expand = c(0, 0)) +
        # set x limits between 0 and 1200
        coord_cartesian(xlim = c(0, 1300)) +
        # remove strip border
        theme(
                # strip background color to lightgrey
                strip.background = element_blank()
        ) +
        # change axis text size to 5
        theme(axis.text = element_text(size = 7)) +
        # axis title size to 7
        theme(axis.title = element_text(size = 7)) +
        geom_text_repel(aes(label = count), position = position_stack(vjust = 0), size = 2, max.overlaps = Inf) +
        scale_fill_manual(values = manual_cols)

# stack p_train and p_test, use the same x axis scale
p <- (p_train + p_test) + plot_layout(guides = "collect", ncol = 1)



plotfile <- file.path("/home/nrlab/wang04/ulyses/0_data_split", "data_split_plot.pdf")

ggsave(plotfile,
        p,
        width = 180,
        height = 210,
        units = "mm"
)
message("Saved plot to ", plotfile)

ans7 <- ans6 %>%
        mutate(y = factor(y, levels = c(
                ans6_test_all, ans6_train_all,
                ans6_test_0.1_1, ans6_train_0.1_1,
                ans6_test_0.03_0.1, ans6_train_0.03_0.1,
                ans6_test_0_0.03, ans6_train_0_0.03
        )))

ans_table_all_cohorts <- ans7 %>%
        mutate(strat_label = paste0(split, ": ", strat))

# ans_table_all_cohorts$strat_label <- factor(ans_table_all_cohorts$strat_label, levels = c(
#         "Test: all", "Train: all",
#         "Test: (0.1, 1]", "Train: (0.1, 1]",
#         "Test: (0.03, 0.1]", "Train: (0.03, 0.1]",
#         "Test: [0, 0.03]", "Train: [0, 0.03]"
# ))

ans_table_all_cohorts$strat_label <- factor(ans_table_all_cohorts$strat_label, levels = c(
        "Train: [0, 0.03]", "Test: [0, 0.03]",
        "Train: (0.03, 0.1]", "Test: (0.03, 0.1]",
        "Train: (0.1, 1]", "Test: (0.1, 1]",
        "Train: all", "Test: all"
))

ans_table_all_cohorts <- ans_table_all_cohorts %>%
        select(cohort, count, strat_label) %>%
        mutate(cohort = fct_relevel(cohort, "Healthy")) %>%
        # sort the strat_label
        arrange(cohort, strat_label) %>%
        # pivot wider, make strat_label as columns, count as values
        pivot_wider(names_from = strat_label, values_from = count, values_fill = 0)
# set the Healthy to the first factor level in cohort

# add a sum row to the end of the table

summary_row <- ans_table_all_cohorts %>%
        summarise(across(where(is.numeric), sum, na.rm = TRUE))

# # Add a label for the summary row (optional)
summary_row <- summary_row %>%
        mutate(cohort = "Total")

# # Bind the summary row to the original data
ans_table_all_cohorts <- bind_rows(ans_table_all_cohorts, summary_row)

subheader <- str_extract(colnames(ans_table_all_cohorts), "cohort|Train|Test") %>%
        t() %>%
        as_tibble()
# set column names
colnames(subheader) <- colnames(ans_table_all_cohorts)
subheader$cohort <- " "
subheader <- subheader %>% mutate(across(everything(), as.character))
ans_table_all_cohorts <- ans_table_all_cohorts %>%
        mutate(across(everything(), as.character))

# remove train: or Test: string from the column names

# add subheader as 2nd row to the table

df_with_subheader <- ans_table_all_cohorts %>%
        slice(1) %>% # Get the first row
        bind_rows(subheader, .) %>% # Add the subheader row after the first row
        bind_rows(slice(ans_table_all_cohorts, -1))

colnames(df_with_subheader) <- str_replace_all(colnames(df_with_subheader), "Train: ", "")
colnames(df_with_subheader) <- str_replace_all(colnames(df_with_subheader), "Test: ", "")


cohort_split_table <- file.path("/home/nrlab/wang04/ulyses/0_data_split", "data_split_cohort_split_table.csv")
write_csv(df_with_subheader, cohort_split_table)
message("Saved table to ", cohort_split_table)


################################################################################
# plot one the same plot
################################################################################



# plot stacked bar, y is `all`, stacked by cohort
p_one <- ggplot(ans7, aes(y = y, x = count, fill = cohort)) +
        geom_bar(position = "stack", stat = "identity", width = 0.5) +
        # facet_wrap(~split, ncol = 1) +
        theme_classic() +
        # put legend at the bottom
        theme(legend.position = "bottom") +
        labs(
                x = "Count (n)",
                y = "ichorCNA TF strata",
                fill = "Strat all"
        ) +
        # remove y axis
        theme(axis.title.y = element_blank()) +
        # manually set the color of fill
        # x axis expand to 0
        scale_x_continuous(expand = c(0, 0)) +
        # y axis expand to 0
        scale_y_discrete(expand = c(0, 0)) +
        # set x limits between 0 and 1200
        coord_cartesian(xlim = c(0, 1300)) +
        # remove strip border
        theme(
                # strip background color to lightgrey
                strip.background = element_blank()
        ) +
        # change axis text size to 5
        theme(axis.text = element_text(size = 7)) +
        # axis title size to 7
        theme(axis.title = element_text(size = 7)) +
        geom_text_repel(aes(label = count), position = position_stack(vjust = 0), size = 2, max.overlaps = Inf) +
        scale_fill_manual(values = manual_cols) +
        # make the legend smaller
        theme(legend.key.size = unit(0.4, "cm")) +
        theme(legend.text = element_text(size = 5)) +
        theme(legend.title = element_text(size = 5))



plotfile <- file.path("/home/nrlab/wang04/ulyses/0_data_split", "data_split_plot_one.pdf")

ggsave(plotfile,
        p_one,
        width = 180,
        height = 100,
        units = "mm"
)
message("Saved plot to ", plotfile)


ans_table_all_cohorts <- ans7 %>%
        mutate(strat_label = paste0(split, ": ", y)) %>%
        select(cohort, count, strat_label) %>%
        # pivot wider, make strat_label as columns, count as values
        pivot_wider(names_from = strat_label, values_from = count)

cohort_split_table <- file.path("/home/nrlab/wang04/ulyses/0_data_split", "data_split_cohort_split_table.csv")
write_csv(ans_table_all_cohorts, cohort_split_table)
message("Saved table to ", cohort_split_table)


################################################################################
# plot Healthy vs Cancer
################################################################################

# reain the test data from final_test_hyper
load(final_test_hyper)
all_filtered <- readRDS(xgboost_input_data_file_hyper) %>%
        filter(assay == "sd") %>%
        select(-rowname, -colname, -assay, -value) %>%
        distinct()
final_test_bamid <- final_test_primary_vec


ans <- all_filtered %>%
        mutate(split = ifelse(primary %in% final_test_bamid, "Test", "Train")) %>%
        # change the cohort to Healthy and Cancer
        mutate(cohort = case_when(
                cohort == "Healthy" ~ "Healthy",
                TRUE ~ "Cancer"
        ))

ans2 <- ans %>%
        mutate(`all` = TRUE) %>%
        mutate(`[0, 0.03]` = case_when(
                ichorcna_tf_strat %in% c("Healthy", "[0, 0.03]") ~ TRUE,
                TRUE ~ FALSE
        )) %>%
        mutate(`(0.03, 0.1]` = case_when(
                ichorcna_tf_strat %in% c("Healthy", "(0.03, 0.1]") ~ TRUE,
                TRUE ~ FALSE
        )) %>%
        mutate(`(0.1, 1]` = case_when(
                ichorcna_tf_strat %in% c("Healthy", "(0.1, 1]") ~ TRUE,
                TRUE ~ FALSE
        ))

# pivot longer, gather `all`, `[0, 0.03]`, `(0.03, 0.1]`, `(0.1, 1]` to a column called "strat", values to a column called "strat_bool"
ans3 <- ans2 %>%
        pivot_longer(
                cols = c(`all`, `[0, 0.03]`, `(0.03, 0.1]`, `(0.1, 1]`),
                names_to = "strat",
                values_to = "strat_bool"
        )


ans4 <- ans3 %>%
        group_by(cohort, split, strat, strat_bool) %>%
        # summarise how many rows with strat_bool == TRUE
        summarise(count = n()) %>%
        ungroup() %>%
        filter(strat_bool == TRUE)

# set the "Healthy" to the first factor level in cohort
ans4$cohort <- factor(ans4$cohort)

ans4$cohort <- relevel(ans4$cohort, "Healthy")

# reverse the factor level of cohort
ans4$cohort <- fct_rev(ans4$cohort)

# set the factor levels of strat to `[0, 0.03]`, `(0.03, 0.1]`, `(0.1, 1]` and `all`
ans4$strat <- factor(ans4$strat, levels = c("[0, 0.03]", "(0.03, 0.1]", "(0.1, 1]", "all"))

# reverse the factor level of strat
ans4$strat <- fct_rev(ans4$strat)

# set the factor levels of split to "Train" and "Test"
ans4$split <- factor(ans4$split, levels = c("Train", "Test"))


# set color to distinct color, the number of colors equals to the number of cohorts in ans4
n_colors <- length(unique(ans4$cohort))
# use distinct colors for categorical data

qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# set seed to 0
set.seed(1)
manual_cols <- sample(col_vector, n_colors)
names(manual_cols) <- unique(ans4$cohort)
# change the color of Healthy to the "black"
manual_cols["Healthy"] <- "darkgrey"


ans5 <- ans4 %>%
        group_by(split, strat) %>%
        mutate(total_count = sum(count)) %>%
        # calculate thte sum of "Healthy" cohort
        mutate(total_control = sum(count[cohort == "Healthy"])) %>%
        # calculate the sum of non-healthy
        mutate(total_case = total_count - total_control) %>%
        # drop group_by
        ungroup()
ans6 <- ans5 %>%
        mutate(y = paste0(strat, "\n", "Ctrl:", total_control, "; ", "Case:", total_case))
# train
ans6_train <- ans6 %>%
        filter(split == "Train")

ans6_train_all <- ans6_train %>%
        filter(strat == "all") %>%
        pull(y) %>%
        unique()
# make ans6_train_0_0.03
ans6_train_0_0.03 <- ans6_train %>%
        filter(strat == "[0, 0.03]") %>%
        pull(y) %>%
        unique()

# make ans6_train_0.03_0.1
ans6_train_0.03_0.1 <- ans6_train %>%
        filter(strat == "(0.03, 0.1]") %>%
        pull(y) %>%
        unique()

ans6_train_0.1_1 <- ans6_train %>%
        filter(strat == "(0.1, 1]") %>%
        pull(y) %>%
        unique()

# set levels of y to  ans6_train_0_0.03, ans6_train_0.03_0.1, ans6_train_0.1_1, ans6_train_all
ans6_train$y <- factor(ans6_train$y, levels = c(ans6_train_0_0.03, ans6_train_0.03_0.1, ans6_train_0.1_1, ans6_train_all))
# reverse the factor level of y
ans6_train$y <- fct_rev(ans6_train$y)

# test

ans6_test <- ans6 %>%
        filter(split == "Test")

ans6_test_all <- ans6_test %>%
        filter(strat == "all") %>%
        pull(y) %>%
        unique()
# make ans6_train_0_0.03
ans6_test_0_0.03 <- ans6_test %>%
        filter(strat == "[0, 0.03]") %>%
        pull(y) %>%
        unique()

# make ans6_train_0.03_0.1
ans6_test_0.03_0.1 <- ans6_test %>%
        filter(strat == "(0.03, 0.1]") %>%
        pull(y) %>%
        unique()

ans6_test_0.1_1 <- ans6_test %>%
        filter(strat == "(0.1, 1]") %>%
        pull(y) %>%
        unique()

# set levels of y to  ans6_train_0_0.03, ans6_train_0.03_0.1, ans6_train_0.1_1, ans6_train_all
ans6_test$y <- factor(ans6_test$y, levels = c(ans6_test_0_0.03, ans6_test_0.03_0.1, ans6_test_0.1_1, ans6_test_all))

# reverse the factor level of y
ans6_test$y <- fct_rev(ans6_test$y)

ans7 <- ans6 %>%
        mutate(y = factor(y, levels = c(
                ans6_test_all, ans6_train_all,
                ans6_test_0.1_1, ans6_train_0.1_1,
                ans6_test_0.03_0.1, ans6_train_0.03_0.1,
                ans6_test_0_0.03, ans6_train_0_0.03
        )))

ans8 <- ans7 %>%
        mutate(strat_label = paste0(split, ": ", strat))

ans8$strat_label <- factor(ans8$strat_label, levels = c(
        "Test: all", "Train: all",
        "Test: (0.1, 1]", "Train: (0.1, 1]",
        "Test: (0.03, 0.1]", "Train: (0.03, 0.1]",
        "Test: [0, 0.03]", "Train: [0, 0.03]"
))

# set cohort levels to c("Healthy", "Caner")
ans8$cohort <- factor(ans8$cohort, levels = c("Cancer", "Healthy"))
# plot stacked bar, y is `all`, stacked by cohort
p_one_cancer_healthy <- ggplot(ans8, aes(y = strat_label, x = count, fill = cohort)) +
        geom_bar(position = "stack", stat = "identity") +
        # facet_wrap(~split, ncol = 1) +
        theme_classic() +
        # put legend at the bottom
        theme(legend.position = "right") +
        labs(
                x = "Count (n)",
                y = "ichorCNA TF strata",
                fill = "Strat all"
        ) +
        # remove y axis
        theme(axis.title.y = element_blank()) +
        # manually set the color of fill
        # x axis expand to 0
        scale_x_continuous(expand = c(0, 0)) +
        # y axis expand to 0
        scale_y_discrete(expand = c(0, 0)) +
        # set x limits between 0 and 1200
        coord_cartesian(xlim = c(0, 1300)) +
        # remove strip border
        theme(
                # strip background color to lightgrey
                strip.background = element_blank()
        ) +
        # change axis text size to 5
        theme(axis.text = element_text(size = 7)) +
        # axis title size to 7
        theme(axis.title = element_text(size = 7)) +
        # annotate the count on the bar not using ggrepel
        geom_text(aes(label = count), position = position_stack(vjust = 0.5), size = 2) +
        scale_fill_manual(values = c("Healthy" = "darkgrey", "Cancer" = "tomato3")) +
        # make the legend smaller
        theme(legend.key.size = unit(0.4, "cm")) +
        theme(legend.text = element_text(size = 5)) +
        theme(legend.title = element_text(size = 5)) +
        # set legend title to "Cohort"
        guides(fill = guide_legend(title = "Cohort"))



plotfile <- file.path("/home/nrlab/wang04/ulyses/0_data_split", "data_split_plot_one_cancer_and_healthy.pdf")

ggsave(plotfile,
        p_one_cancer_healthy,
        width = 100,
        height = 70,
        units = "mm"
)
message("Saved plot to ", plotfile)



################################################################################

p_one_cancer_healthy_two_col <- ggplot(ans8, aes(y = strat, x = count, fill = cohort)) +
        geom_bar(position = "stack", stat = "identity") +
        # facet_wrap(~split, ncol = 1) +
        theme_classic() +
        # put legend at the bottom
        labs(
                x = "Count (n)",
                y = "ichorCNA TF strata",
                fill = "Strat all"
        ) +
        # remove y axis
        theme(axis.title.y = element_blank()) +
        # manually set the color of fill
        # x axis expand to 0
        scale_x_continuous(expand = c(0, 0)) +
        # y axis expand to 0
        scale_y_discrete(expand = c(0, 0)) +
        # set x limits between 0 and 1200
        coord_cartesian(xlim = c(0, 1300)) +
        # remove strip border
        theme(
                # strip background color to lightgrey
                strip.background = element_blank()
        ) +
        # change axis text size to 5
        theme(axis.text = element_text(size = 5)) +
        # axis title size to 7
        theme(axis.title = element_text(size = 5)) +
        # annotate the count on the bar not using ggrepel
        geom_text(aes(label = count), position = position_stack(vjust = 0.5), size = 2) +
        scale_fill_manual(values = c("Healthy" = "darkgrey", "Cancer" = "tomato3"), drop = TRUE) +
        # make the legend smaller
        theme(legend.key.size = unit(0.4, "cm")) +
        theme(legend.text = element_text(size = 5)) +
        theme(legend.title = element_text(size = 5)) +
        # set legend title to "Cohort"
        guides(fill = guide_legend(title = "Cohort")) +
        theme(legend.position.inside = c(-0.3, 0.5)) +
        # change facet strip text size to 5
        theme(strip.text = element_text(size = 5)) +
        facet_wrap(~split, ncol = 2) +
        # remove legend
        theme(legend.position = "none")



plotfile <- file.path("/home/nrlab/wang04/ulyses/0_data_split", "data_split_plot_one_cancer_and_healthy_two_col.pdf")

ggsave(plotfile,
        p_one_cancer_healthy_two_col,
        width = 90,
        height = 50,
        units = "mm"
)
message("Saved plot to ", plotfile)

################################################################################
# prepare a summary table
################################################################################

ans_table_cancer_healthy <- ans8 %>%
        select(cohort, count, strat_label) %>%
        # pivot wider, make strat_label as columns, count as values
        pivot_wider(names_from = strat_label, values_from = count)

cancer_healthy_table <- file.path("/home/nrlab/wang04/ulyses/0_data_split", "data_split_cancer_healthy_table.csv")
write_csv(ans_table_cancer_healthy, cancer_healthy_table)
message("Saved table to ", cancer_healthy_table)




################################################################################
# test train and test tf
################################################################################

# add_ichorcna_tf_strat_col2 <- function(input) {
#         ans <- input |>
#                 dplyr::mutate(ichorcna_tf_strat = dplyr::case_when(
#                         ichorcna_tf <= 0.03 ~ "[0, 0.03]",
#                         ichorcna_tf <= 0.1 & ichorcna_tf > 0.03 ~ "(0.03, 0.1]",
#                         ichorcna_tf <= 1 & ichorcna_tf > 0.1 ~ "(0.1, 1]",
#                         TRUE ~ ichorcna_tf
#                 ))
#         return(ans)
# }


# test_bamid <- read_csv("/scratchc/nrlab/wang04/ulyses_results_iteration/iteration_3_merge_below_3/xgboost_model/xgb10/all/xgboost_input_independent_test.feat.sd.tf.all.csv") |>
#         pull(primary) |>
#         unique()

# load("/home/nrlab/wang04/ulyses/models/final_test_samples.rda")


train_color = "#ffa600"
test_color = "#8a508f"
coldata_tibble2 <- all_filtered %>%
        mutate(train_test = ifelse(primary %in% final_test_bamid, "Test", "Train"))

# coldata_tibble2 <- add_ichorcna_tf_strat_col2(coldata_tibble) |>
#         mutate(train_test = ifelse(primary %in% test_bamid, "Test", "Train"))

# if cohort is "Healthy", set ichorcna_tf_strat to "Healthy"
# coldata_tibble2$ichorcna_tf_strat <- ifelse(coldata_tibble2$cohort == "Healthy", "Healthy", coldata_tibble2$ichorcna_tf_strat)

# set factor levels of ichorcna_tf_strat
# coldata_tibble2$ichorcna_tf_strat <- factor(coldata_tibble2$ichorcna_tf_strat, levels = c("Healthy", "[0, 0.03]", "(0.03, 0.1]", "(0.1, 1]"))
coldata_tibble2$train_test <- factor(coldata_tibble2$train_test, levels = c("Train", "Test"))


tf_coldata_tibble <- select(coldata_tibble2, primary, ichorcna_tf, train_test, ichorcna_tf_strat)
tf_coldata_tibble <- tf_coldata_tibble |>
        group_by(ichorcna_tf_strat, train_test) |>
        mutate(count = n()) |>
        ungroup() |>
        mutate(x_label = paste0(train_test, "\n", "(", count, ")")) |>
        arrange(ichorcna_tf_strat, train_test)
x_label_levels <- tf_coldata_tibble$x_label |> unique()
tf_coldata_tibble$x_label <- factor(tf_coldata_tibble$x_label, levels = x_label_levels)


p_train_test_tf <- ggplot(tf_coldata_tibble, mapping = aes(train_test, ichorcna_tf)) +
        geom_violin(fill = "lightgrey", color = "lightgrey", width = 1) +
        # boxplot without outliers
        geom_boxplot(aes(fill = train_test), color = "black", outlier.shape = NA, width = 0.3, linewidth = 0.2) +
        geom_jitter(
                color = "black",
                fill = "#e9e6e6",
                alpha = 0.9,
                width = 0.2,
                size = 0.2,
                # shape = 21,
                stroke = 0.01
        ) +
        # boxplot without outliers
        labs(x = "Tumour Fraction Stratification", y = "Tumour Fraction") +
        scale_fill_manual(values = c("Train" = train_color, "Test" = test_color)) +
        ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif", vjust = 0.5, ref.group = "Train", size = pvalue_font_size_hyper) +
        theme_classic() +
        theme_NPJ() +
        # rotate x axis 45
        theme(axis.text = element_text(size = 5.5)) +
        # axis title size to 6
        theme(axis.title = element_text(size = 7)) +
        # strip text size to 6
        theme(strip.text = element_text(size = 5.5)) +
        # strip background rect boarder linewidth to 0.1
        theme(strip.background = element_rect(linewidth = 0, fill = "lightgrey", color = NA)) +
        # rotate x axis labels 45 degree
        # remove x axis title
        theme(axis.title.x = element_blank()) +
        # hide legend
        theme(legend.position = "none") +
        # change x axis label to x_label
        # scale_x_discrete(labels = tf_coldata_tibble$x_label) +
        facet_wrap(~ichorcna_tf_strat, nrow = 1, scales = "free", drop = TRUE)


ggsave(file.path("/home/nrlab/wang04/ulyses/0_data_split", "train_test_tf.pdf"), plot = p_train_test_tf, width = 3.6 * 3, height = 3, units = "cm")
message("Saved plot to ", file.path("/home/nrlab/wang04/ulyses/0_data_split", "train_test_tf.pdf"))


# plot the all group
all_cancer_data <- tf_coldata_tibble %>%
        filter(ichorcna_tf_strat != "Healthy") %>%
        mutate(ichorcna_tf_strat = "All")

p_train_test_tf_all <- ggplot(all_cancer_data, mapping = aes(train_test, ichorcna_tf)) +
        geom_violin(fill = "lightgrey", color = "lightgrey", width = 1) +
        # boxplot without outliers
        geom_boxplot(aes(fill = train_test), color = "black", outlier.shape = NA, width = 0.3, linewidth = 0.2) +
        geom_jitter(
                color = "black",
                fill = "#e9e6e6",
                alpha = 0.9,
                width = 0.2,
                size = 0.2,
                # shape = 21,
                stroke = 0.01
        ) +
        # boxplot without outliers
        labs(x = "Tumour Fraction Stratification", y = "Tumour Fraction") +
        scale_fill_manual(values = c("Train" = train_color, "Test" = test_color)) +
        ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif", vjust = 0.5, ref.group = "Train", size = pvalue_font_size_hyper) +
        theme_classic() +
        theme_NPJ() +
        # rotate x axis 45
        theme(axis.text = element_text(size = 5.5)) +
        # axis title size to 6
        theme(axis.title = element_text(size = 7)) +
        # strip text size to 6
        theme(strip.text = element_text(size = 5.5)) +
        # strip background rect boarder linewidth to 0.1
        theme(strip.background = element_rect(linewidth = 0, fill = "lightgrey", color = NA)) +
        # rotate x axis labels 45 degree
        # remove x axis title
        theme(axis.title.x = element_blank()) +
        # hide legend
        theme(legend.position = "none") +
        # change x axis label to x_label
        # scale_x_discrete(labels = tf_coldata_tibble$x_label) +
        facet_wrap(~ichorcna_tf_strat, nrow = 1, scales = "free", drop = TRUE)


ggsave(file.path("/home/nrlab/wang04/ulyses/0_data_split", "train_test_tf_all.pdf"), plot = p_train_test_tf_all, width = 3.6, height = 3, units = "cm")
message("Saved plot to ", file.path("/home/nrlab/wang04/ulyses/0_data_split", "train_test_tf_all.pdf"))




################################################################################
# train test PCA of different tf categories of different features. (as a grid)
################################################################################

xgb_base_path <- "/scratchc/nrlab/wang04/ulyses_results_iteration/iteration_3_merge_below_3/xgboost_model/xgb10"

sd_feat_filename <- "xgboost_input\\.feat\\.sd\\.tf\\..*\\.csv"
cnv_feat_filename <- "xgboost_input\\.feat\\.cnv\\.tf\\..*\\.csv"
length_feat_filename <- "xgboost_input\\.feat\\.length\\.tf\\..*\\.csv"
ctRatio_feat_filename <- "xgboost_input\\.feat\\.ctRatio\\.tf\\..*\\.csv"
slRatio_feat_filename <- "xgboost_input\\.feat\\.slRatio\\.tf\\..*\\.csv"
patterns <- list(
        sd_feat_filename,
        cnv_feat_filename,
        length_feat_filename,
        ctRatio_feat_filename,
        slRatio_feat_filename
)

feats <- lapply(patterns, list.files, path = xgb_base_path, full.names = TRUE, recursive = TRUE)
names(feats) <- c("sd", "cnv", "length", "ctRatio", "slRatio")


feats_read <- lapply(feats, lapply, read_csv, id = "filename")


feats_read_bind <- lapply(feats_read, bind_rows)

# write a function to add a tf_strata column
add_strata_train <- function(x) {
        x <- x %>%
                mutate(tf_strata = case_when(
                        str_detect(filename, "0-0.03") ~ "[0, 0.03]",
                        str_detect(filename, "0.03-0.1") ~ "(0.03, 0.1]",
                        str_detect(filename, "0.1-1") ~ "(0.1, 1]",
                        str_detect(filename, "all") ~ "all",
                        TRUE ~ "unknown"
                )) %>%
                select(-filename, -bicohort, -patient_id, -primary) %>%
                mutate(train_test = "Train")
}

feats_read_bind_clean <- lapply(feats_read_bind, add_strata_train)



# handling test data------------------------------------------------------------



sd_feat_filename <- "xgboost_input(_independent_test)?\\.feat\\.sd\\.tf\\..*\\.csv"
cnv_feat_filename <- "xgboost_input(_independent_test)?\\.feat\\.cnv\\.tf\\..*\\.csv"
length_feat_filename <- "xgboost_input(_independent_test)?\\.feat\\.length\\.tf\\..*\\.csv"
ctRatio_feat_filename <- "xgboost_input(_independent_test)?\\.feat\\.ctRatio\\.tf\\..*\\.csv"
slRatio_feat_filename <- "xgboost_input(_independent_test)?\\.feat\\.slRatio\\.tf\\..*\\.csv"
patterns <- list(
        sd_feat_filename,
        cnv_feat_filename,
        length_feat_filename,
        ctRatio_feat_filename,
        slRatio_feat_filename
)

feats <- lapply(patterns, list.files, path = xgb_base_path, full.names = TRUE, recursive = TRUE)
names(feats) <- c("sd", "cnv", "length", "ctRatio", "slRatio")


feats_read <- lapply(feats, lapply, read_csv, show_col_types = FALSE)

# for loop
for (i in names(feats_read)) {
        names(feats_read[[i]]) <- feats[[i]]
}

feats_read_bind <- lapply(feats_read, bind_rows, .id = "filename")

# write a function to add a tf_strata column
add_strata_train <- function(x) {
        x <- x %>%
                mutate(tf_strata = case_when(
                        str_detect(filename, "0-0.03") ~ "[0, 0.03]",
                        str_detect(filename, "0.03-0.1") ~ "(0.03, 0.1]",
                        str_detect(filename, "0.1-1") ~ "(0.1, 1]",
                        str_detect(filename, "all") ~ "all",
                        TRUE ~ "unknown"
                )) %>%
                mutate(train_test = case_when(
                        str_detect(filename, "_independent_test") ~ "Test",
                        !str_detect(filename, "_independent_test") ~ "Train",
                        TRUE ~ "unknown"
                )) %>%
                select(-filename, -patient_id)
        
        # set the col 'primary' as rownames
        # x <- column_to_rownames(x, var = "primary")

}

feats_read_bind_clean <- lapply(feats_read_bind, add_strata_train)


lapply(feats_read_bind_clean, FUN = function(x) {
        x %>%
                group_by(tf_strata, train_test) %>%
                summarise(n = n())
})



# feats_read_bind_flat <- bind_rows(feats_read_bind, .id = "feat_label")



pca_viz <- function(which_feat = "sd", tf_strata = "all",
                    inputData = feats_read_bind_clean,
                    plot_width = 9,
                    plot_height = 6.5,
                    train_color = "royalblue",
                    test_color = "tomato3",
                    plot_saving_dir = "/home/nrlab/wang04/ulyses/0_data_split") {
        # Extract the data frame corresponding to the specified feature
        ans <- inputData[[which_feat]] %>%
                filter(`tf_strata` == UQ(tf_strata)) %>%
                # for clarity purpose
                filter(primary != "EE87123")


        # Select all columns except 'train_test' for PCA analysis
        active <- ans |>
                select(-train_test, -tf_strata, -bicohort) |>
                select(-matches("ctRatio\\.19|slRatio\\.19|cnv\\.19")) |>
                column_to_rownames(var = "primary")


        # remove any cols containing 'ctRatio.19', 'slRatio.19' and 'cnv.19'

        # Extract the 'train_test' column for color grouping in PCA plot
        train_test_vec <- ans$train_test
        train_test_vec <- factor(train_test_vec, levels = c("Train", "Test"))

        # Perform PCA on the selected features
        res.pca <- prcomp(active)

        # Create a PCA plot with specific visualization settings
        p <- factoextra::fviz_pca_ind(res.pca,
                geom = c("point"),
                # col.ind = train_test_vec, # color points by 'train_test'
                palette = c("Train" = train_color, "Test" = test_color),
                alpha.ind = 0.7,
                # habillage = ans$bicohort,
                habillage = ans$train_test,
                legend.title = "Data split",
                pointshape = 21,
                title = paste0(which_feat, ": ", tf_strata),
                repel = TRUE,
                addEllipses = TRUE,
                ellipse.level = 0.95,
                ellipse.type = "confidence",
                pointsize = 0.4
        ) +
                theme_NPJ() +
                theme(
                        plot.title = element_text(size = 5, hjust = 0), # Set the plot title size to 5.5
                        legend.title = element_text(size = 5), # Set the legend title size to 5.5
                        legend.text = element_text(size = 5), # Set the legend key font size to 5.5
                        legend.key.size = unit(10, "pt"),
                        legend.position = "none"
                ) +
                theme(legend.position = "none")
        

        p_with_text <- factoextra::fviz_pca_ind(res.pca,
                geom = c("point", "text"),
                # col.ind = train_test_vec, # color points by 'train_test'
                palette = c("Train" = train_color, "Test" = test_color),
                alpha.ind = 0.7,
                # habillage = ans$bicohort,
                habillage = ans$train_test,
                legend.title = "Data split",
                pointshape = 21,
                title = paste0(which_feat, ": ", tf_strata),
                repel = FALSE,
                addEllipses = TRUE,
                ellipse.level = 0.95,
                ellipse.type = "confidence",
                pointsize = 0.4
        ) +
                theme_NPJ() +
                theme(
                        plot.title = element_text(size = 5, hjust = 0), # Set the plot title size to 5.5
                        legend.title = element_text(size = 5), # Set the legend title size to 5.5
                        legend.text = element_text(size = 5), # Set the legend key font size to 5.5
                        legend.key.size = unit(10, "pt"),
                        legend.position = "none"
                ) +
                theme(legend.position = "none")



        # Define the filename for saving the plot
        plotfile <- file.path(plot_saving_dir, paste0("train_test_PCA_", which_feat, "_", tf_strata, ".pdf"))
        plotfile_with_text <- file.path(plot_saving_dir, paste0("train_test_PCA_", which_feat, "_", tf_strata, ".with_text.pdf"))
        plotfile_html <- gsub("\\.pdf$", ".html", plotfile_with_text)

        # Save the plot to a PDF file
        ggsave(
                filename = plotfile,
                plot = p,
                width = plot_width,
                height = plot_height,
                dpi = 300,
                units = "cm"
        )
        # Inform the user that the plot has been saved
        message("Saved to ", plotfile)

        # convert plot to plotly 
        p_plotly <- ggplotly(p_with_text)
        saveWidget(p_plotly, file = plotfile_html, selfcontained = TRUE)
        message("saved to html file:", plotfile_html)

        # Save the plot to a PDF file
        ggsave(
                filename = plotfile_with_text,
                plot = p_with_text,
                width = 20,
                height = 20,
                dpi = 300,
                units = "cm"
        )

        # Inform the user that the plot has been saved
        message("Saved to ", plotfile_with_text)
        return(p)
}

which_feats <- names(feats_read_bind_clean)
which_tf_strata <- feats_read_bind_clean[[1]] |>
        pull("tf_strata") |>
        unique()

params <- expand.grid(which_feat = which_feats, tf_strata = which_tf_strata)

results <- params %>%
        pmap(~ pca_viz(
                which_feat = ..1,
                tf_strata = ..2,
                inputData = feats_read_bind_clean,
                train_color = train_color,
                test_color = test_color,
                plot_width = 18 / 5,
                plot_height = 12.4 / 4
        ))


