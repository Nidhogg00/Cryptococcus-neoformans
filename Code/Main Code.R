library(tidyverse)
library(readr)
library(fields)
library(reshape2)
library(patchwork)
library(car) # for vif
library(factoextra) # for PCA
library(MASS)
library(dplyr)
library(MASS)
library(VennDiagram)
library(DESeq2)
library(glmnet)
library(ggplot2)
library(GGally)
library(mice)
library(pheatmap)
library(xtable)
library(mpath)

# ------------------------------ #
# Load Data
# ------------------------------ #
counts_all <- read_csv("data/counts_all.csv")
Gat201_samplesheet <- read_csv("data/Gat201_samplesheet.csv")
H99_all_genes_promoter_500nt_2mer_counts <- 
    read_csv("data/H99_all_genes_promoter_500nt_2mer_counts.csv", skip = 10)
H99_all_genes_promoter_500nt_3mer_counts <- 
    read_csv("data/H99_all_genes_promoter_500nt_3mer_counts.csv", skip = 10)
H99_all_genes_promoter_500nt_4mer_counts <- 
    read_csv("data/H99_all_genes_promoter_500nt_4mer_counts.csv", skip = 10)
H99_all_genes_promoter_500nt_5mer_counts <- 
    read_csv("data/H99_all_genes_promoter_500nt_5mer_counts.csv")

# ------------------------------ #
# Exploratory Data Analysis
# ------------------------------ #

# Extract specific columns from counts_all for further analysis
rawY <- counts_all %>% 
    dplyr::select(!c(2:6)) %>%  # Remove columns 2 to 6
    dplyr::select(
        contains("Geneid"),
        contains("_0_"),
        contains("_30_"),
        contains("_120_"),
        contains("_240_")
    )

# Display summary statistics of the selected data (excluding the first column)
summary(unlist(rawY[, -1]))


# ------------------------------ #
# Pairwise Scatter Plot Matrix
# ------------------------------ #

# Select the relevant columns from the counts_all data frame
selected_data <- 
    counts_all[, c("A_R_120_1", "A_R_120_2", "a_R_120_1", "a_R_120_2")]

# Create a scatter plot matrix using the ggpairs function from the 
#GGally package
scatter_plot_matrix <- ggpairs(selected_data,
                               lower = list
                               (continuous = wrap("points", color = "blue", 
                                                  size = 1)),
                               diag = list(continuous = 
                                               wrap("densityDiag", 
                                                    fill = "lightblue")),
                               upper = list(continuous = wrap("cor", size = 4)),
                               title = "Pairwise Scatter Plot Matrix")

# Print the scatter plot matrix
print(scatter_plot_matrix)


# ------------------------------ #
# Principal Component Analysis (PCA)
# ------------------------------ #

# Retrieve the order of column names from the samplesheet
order <- Gat201_samplesheet$Title

# Ensure the first column is 'Geneid' followed by the order from the samplesheet
order <- c("Geneid", order)

# Reorder the columns of rawY according to the specified order
rawY <- rawY %>% dplyr::select(all_of(order))

# Check if the column names of rawY match the specified order
colnames(rawY) == order

# Perform PCA on the transposed data (excluding the Geneid column)
pca <- prcomp(t(rawY %>% dplyr::select(!Geneid)))

# Combine PCA projections with the sample sheet
Title_pc <- data.frame(Title = colnames(rawY)[-1], pca$x)

# Create a data frame with the proportion of variance 
#explained by each principal component
propvar_df <- tibble(
    PC = seq.int(1L, ncol(pca$x)),
    prop_var = pca$sdev^2 / sum(pca$sdev^2)
)

# Plot the proportion of variance explained by each principal component
plot_percentvar <- 
    ggplot(data = propvar_df, 
           aes(x = PC, y = prop_var)) +
    geom_col(fill = "blue") +
    scale_x_continuous("Principal Component",
                       limits = c(0.4, 10.6), 
                       breaks = 1L:10L,
                       expand = c(0, 0)) + 
    scale_y_continuous("Proportion of Variance", expand = c(0, 0)) +
    theme(panel.grid.major = element_blank())

# Print the plot showing the proportion of variance explained 
#by each principal component
plot_percentvar

# Combine the sample sheet with PCA projections
df_PC <- bind_cols(Gat201_samplesheet %>%
                       dplyr::mutate()) %>%
    left_join(Title_pc, by = "Title")

# Plot PCA results, color-coded by time and shaped by GAT201 status
ggplot(data = df_PC, aes(colour = factor(Time), shape = GAT201)) +
    geom_point(aes(x = PC1, y = PC2)) +
    theme_bw()

# Display summary statistics for the PCA, including the proportion 
#of variance explained
summary(pca)
summary(pca)$importance[2, ]

# ------------------------------ #
# PCA Analysis for Time 240
# ------------------------------ #

# Filter the data frame to include only samples with a time of 240
filtered_pcdf <- df_PC %>% 
    filter(Time == 240)

# Plot PCA results for samples at time 240, color-coded by GAT201 status
ggplot(data = filtered_pcdf,
       aes(colour = GAT201)) +
    geom_point(aes(x = PC1, y = PC2)) +
    labs(title = "PCA of Time 240 Conditions",
         x = "Principal Component 1",
         y = "Principal Component 2") +
    theme_minimal()



# ------------------------------ #
# Outlier Detection and Imputation
# ------------------------------ #

# Check for Outliers using Boxplots -------------------

# Generate a boxplot to visualize gene expression counts for each sample
plotBoxplot_rawY <- rawY %>% 
    dplyr::select(!Geneid) %>% 
    gather(key = "Sample", value = "GE") %>% 
    ggplot(aes(x = factor(Sample, levels = 
                              colnames(rawY %>% dplyr::select(!Geneid))), 
               y = GE)) +
    geom_boxplot(fill = "white", color = "blue") +   
    labs(title = "Gene Expression Boxplot") +
    xlab("") + 
    ylab("Gene Expression Count") +
    scale_y_continuous(expand = c(0.01, 0)) +
    coord_flip() +  # Flip coordinates for better visualization
    theme_bw() +  
    theme(axis.text.x = element_text(angle = 0, hjust = 1))

# Generate a log-transformed boxplot for gene expression counts at 
#0 and 30 time points
plotBoxplot_rawY_log <- rawY %>% 
    dplyr::select(!Geneid) %>%  
    log1p() %>%  # Apply log transformation (log1p handles log(0) as well)
    dplyr::select(contains("_0_"), contains("_30_")) %>%
    gather(key = "Sample", value = "GE") %>% 
    ggplot(aes(x = factor(Sample, levels = 
                              colnames(rawY %>% dplyr::select(!Geneid))), 
               y = GE)) +
    geom_boxplot(fill = "white", color = "blue") +   
    labs(title = "Gene Expression Boxplot (Log-Transformed)") +
    xlab("") + 
    ylab("Gene Expression Count") +
    scale_y_continuous(expand = c(0.01, 0)) +
    coord_flip() +  # Flip coordinates for better visualization
    theme_bw() +  
    theme(axis.text.x = element_text(angle = 0, hjust = 1))

# Analyze Gene Expression Count Frequencies -------------------

# Count the frequency of each unique gene expression count in the dataset
countFreq <- as.data.frame(table(unlist(rawY[, -1]))) %>% 
    as_tibble() %>% 
    arrange(desc(Freq))  # Sort by frequency in descending order

# Rename columns for clarity
colnames(countFreq) = c("countGE", "freq")  

# Select the top 10 most frequent gene expression counts
df_Plot <- countFreq %>% head(10)

# Generate a bar plot for the frequency of gene expression counts
plotBarplot_rawY <- ggplot(df_Plot, aes(x = countGE, y = freq)) +
    geom_bar(stat = "identity", color = 
                 c("red", rep("white", nrow(df_Plot) - 1)), fill = 
                 "lightblue") +
    labs(title = "Count Frequency") +
    xlab("Gene Expression Count") + 
    ylab("Frequency") +
    scale_y_continuous(expand = c(0.01, 0)) +
    scale_x_discrete(limits = rev(df_Plot$countGE)) +  
    # Reverse the x-axis for better readability
    coord_flip() +  # Flip coordinates for better visualization
    theme_bw() +  
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.text.y = element_text(angle = 0, hjust = 0.5))

# Combine and display the boxplots and bar plot for outlier detection
plotBoxplot_rawY | (plotBoxplot_rawY_log / plotBarplot_rawY)

# Correlation Analysis ------------------- 

# Compute the correlation matrix for the gene expression data
corY <- cor(rawY %>% dplyr::select(!Geneid))

# Reshape the correlation matrix for visualization
melted_corr <- melt(corY, varnames = c("Var1", "Var2"))

# Heatmap of the Correlation Matrix ------------------- 

# Generate a heatmap to visualize the correlations between samples
plotHeat_cor <- ggplot(melted_corr, aes(x = Var2, y = Var1, fill = value)) +  
    geom_tile() +  
    scale_fill_gradient2(low = "cyan", mid = "white", high = "red",   
                         midpoint = 0.5, limit = c(0, 1),   
                         space = "Lab", name = "Correlation\nCoefficient") +  
    xlab("") + 
    ylab("") +
    theme_bw() +  
    theme(axis.text.x = element_text(angle = 90, hjust = 1))  
# Rotate x-axis labels for better readability
plotHeat_cor

# Clustered Heatmap ------------------- 

# Generate a clustered heatmap of the correlation matrix
plotHeatClust_cor <- pheatmap(corY, 
                              cluster_rows = TRUE, 
                              cluster_cols = TRUE,   
                              fontsize = 8, 
                              cellwidth = 10, 
                              cellheight = 10)
plotHeatClust_cor



# ------------------------------ #
# Outlier Detection and Imputation
# ------------------------------ #

# Create a list to hold gene expression data for different time groups
list_count_group <- list()
list_count_group[[1]] <- rawY %>% 
    dplyr::select(contains("_0_")) %>% log1p()  
# Log-transform the data for the 0 minute group
list_count_group[[2]] <- rawY %>% 
    dplyr::select(contains("_30_")) %>% log1p()  
# Log-transform the data for the 30 minute group
list_count_group[[3]] <- rawY %>% 
    dplyr::select(contains("_120_") | contains("_240_")) %>% log1p()  
# Log-transform the data for the 120 and 240 minute groups

# Iterate over each time group to detect and handle outliers
for(don in 1:length(list_count_group)){
    count_group <- list_count_group[[don]]
    for(i in 1:ncol(count_group)){
        
        # 1. Principal Component Analysis (PCA) -------------------
        df_pca <- count_group %>% dplyr::select(-i)  
        # Exclude the current column from PCA
        pca_gr <- prcomp(count_group, scale. = TRUE)  
        # Perform PCA on the current group
        sumry <- summary(pca_gr)
        cum_prop_contrib <- sumry$importance[3,]  
        # Calculate cumulative proportion of variance explained
        
        # Determine the number of principal components to retain, 
        #ensuring cumulative variance > 90%
        threshold <- 0.90
        num_selPC <- which(cum_prop_contrib >= threshold)[1]
        pca_gr_sel <- pca_gr$x %>% as_tibble() %>%
            dplyr::select(all_of(1:num_selPC))
        
        # Prepare data for linear regression using selected principal components
        df_lm <- count_group %>% dplyr::select(all_of(i)) %>%
            bind_cols(pca_gr_sel) %>%
            set_names(c("Y", paste0("PC", 1:num_selPC)))
        reslm <- lm(Y ~ ., data = df_lm)
        
        # 2. Outlier Detection using DFFITS Criterion -------------------
        n = nrow(df_lm)  # Number of observations
        p = ncol(df_lm) - 1  # Number of predictors
        threshold <- 2 * sqrt((p + 1) / n)  # Calculate DFFITS threshold
        val_diffits <- dffits(reslm)  # Compute DFFITS values
        
        # Identify and replace outliers with NA
        idNA <- which(abs(val_diffits) > threshold)
        list_count_group[[don]][idNA, i] <- NA
    }
}

# Combine the corrected data groups back together and 
#reverse the log transformation
Y_NA <- bind_cols(list_count_group); Y_NA <- exp(Y_NA) - 1
Y_NA <- bind_cols(rawY %>% dplyr::select(Geneid), Y_NA)  
# Add back the Geneid column
summary(unlist(Y_NA[, -1]))  # Summary of the corrected data

# Identify rows with NA values (i.e., where outliers were detected)
rowNA <- Y_NA %>% dplyr::filter_all(any_vars(is.na(.)))
View(rowNA)  # View the rows with NA values

# Save the data with outliers replaced by NA
save(Y_NA, file = "data/Y_NA.rda")

# Visualize the corrected gene expression data using a boxplot
plotBoxplot_Y_NA <- Y_NA %>% dplyr::select(!Geneid) %>% 
    gather(key = "Sample", value = "GE") %>% 
    ggplot(aes(x = factor(Sample, levels = colnames
                          (rawY %>% dplyr::select(!Geneid))), y = GE)) +
    geom_boxplot(fill = "white", color = "blue") +   
    labs(title = "Gene Expression Boxplot (Corrected)") +
    xlab("") + 
    ylab("Gene Expression Corrected Count") +
    scale_y_continuous(expand = c(0.01, 0)) +
    coord_flip() +
    theme_bw() +  
    theme(axis.text.x = element_text(angle = 0, hjust = 1))

# Compare the original and corrected boxplots
plotBoxplot_rawY | plotBoxplot_Y_NA

# Outlier Imputation using Multiple Imputation by 
#Chained Equations (MICE) -------------------

# Perform multiple imputation on the data with NA values
mice_Y <- mice(Y_NA %>% dplyr::select(-Geneid),
               m = 5, maxit = 30, method = "pmm")
Y_mice <- complete(mice_Y) %>% as_tibble() %>% 
    bind_cols(Y_NA %>% dplyr::select(Geneid)) %>% 
    dplyr::select(Geneid, everything())

# Save the imputed data
save(Y_mice, file = "data/Y_mice.rda")

# Visualize the imputed gene expression data using a boxplot
plotBoxplot_Y <- Y_mice %>% dplyr::select(!Geneid) %>% 
    gather(key = "Sample", value = "GE") %>% 
    ggplot(aes(x = factor(Sample, levels = 
                              colnames(rawY %>% dplyr::select(!Geneid))), 
               y = GE)) +
    geom_boxplot(fill = "white", color = "blue") +   
    labs(title = "Gene Expression Boxplot (Imputed)") +
    xlab("") + 
    ylab("Gene Expression Corrected Count") +
    scale_y_continuous(expand = c(0.01, 0)) +
    coord_flip() +
    theme_bw() +  
    theme(axis.text.x = element_text(angle = 0, hjust = 1))

# Compare the original and imputed boxplots
plotBoxplot_rawY | plotBoxplot_Y

# Summarize the original and imputed data
summary(unlist(Y_mice[, -1]))
summary(unlist(rawY[, -1]))

# ------------------------------ #
# Data Saving
# ------------------------------ #

# Prepare the sample metadata for saving
X <- Gat201_samplesheet %>% 
    dplyr::select(Title, Group, Strain, GAT201, Condition, Time, BioRep) %>% 
    dplyr::mutate(across(everything(), ~as.factor(.x)))
save(X, file = "data/X.rda")

# Analyze and filter motif counts for 2-mers, 3-mers, 4-mers, and 5-mers
rawM_2mer <- H99_all_genes_promoter_500nt_2mer_counts 
summary(unlist(rawM_2mer[, -1]))
rawM_2mer %>% dplyr::select(-Gene) %>% rowSums()
rawM_2mer %>% dplyr::filter(rowSums(rawM_2mer[,-1]) != 499)
ID_Gene_2mer <- rawM_2mer %>% 
    dplyr::filter(rowSums(rawM_2mer[,-1]) > 400) %>% dplyr::select(Gene)

rawM_3mer <- H99_all_genes_promoter_500nt_3mer_counts 
rawM_3mer %>% dplyr::select(-Gene) %>% rowSums()
rawM_3mer %>% dplyr::filter(rowSums(rawM_3mer[,-1]) != 498)
ID_Gene_3mer <- rawM_3mer %>% 
    dplyr::filter(rowSums(rawM_3mer[,-1]) > 400) %>% dplyr::select(Gene) 

rawM_4mer <- H99_all_genes_promoter_500nt_4mer_counts 
rawM_4mer %>% dplyr::select(-Gene) %>% rowSums()
rawM_4mer %>% dplyr::filter(rowSums(rawM_4mer[,-1]) != 497)
ID_Gene_4mer <- rawM_4mer %>% 
    dplyr::filter(rowSums(rawM_4mer[,-1]) > 400) %>% dplyr::select(Gene)

# Combine gene IDs across the different motif lengths and check consistency
list_GeneID_Motifs <- bind_cols(ID_Gene_2mer, ID_Gene_3mer, ID_Gene_4mer)
resCheck <- apply(list_GeneID_Motifs, 1, function(x) all(x == x[1]))
all(resCheck == TRUE)

# Filter the gene expression data to include only genes with motifs detected
Y <- Y_mice %>% dplyr::filter(Geneid %in% ID_Gene_2mer$Gene)

# Reorder columns according to the sample metadata
order <- Gat201_samplesheet$Title
order <- c("Geneid", order)
Y <- Y %>% dplyr::select(all_of(order))

# Ensure column names match the desired order
colnames(Y) == order

# Save the filtered and ordered gene expression data
save(Y, file = "data/Y.rda")

# Save motif count data for 2-mers, 3-mers, 4-mers, and 5-mers
M_2mers <- H99_all_genes_promoter_500nt_2mer_counts %>% 
    dplyr::filter(rowSums(rawM_2mer[,-1]) > 400) %>%
    dplyr::filter(Gene %in% Y$Geneid)
save(M_2mers, file = "data/M_2mers.rda")

M_3mers <- H99_all_genes_promoter_500nt_3mer_counts %>% 
    dplyr::filter(rowSums(rawM_2mer[,-1]) > 400) %>%
    dplyr::filter(Gene %in% Y$Geneid)
save(M_3mers, file = "data/M_3mers.rda")

M_4mers <- H99_all_genes_promoter_500nt_4mer_counts %>% 
    dplyr::filter(rowSums(rawM_2mer[,-1]) > 400) %>%
    dplyr::filter(Gene %in% Y$Geneid)
save(M_4mers, file = "data/M_4mers.rda")

M_5mers <- H99_all_genes_promoter_500nt_5mer_counts %>% 
    dplyr::filter(rowSums(rawM_2mer[,-1]) > 400) %>%
    dplyr::filter(Gene %in% Y$Geneid)
save(M_5mers, file = "data/M_5mers.rda")


# ------------------------------ #
# Significance of GAT201 Using Negative Binomial Regression
# ------------------------------ #

# Set a random seed for reproducibility
set.seed(5)

# Combine 'Condition' and 'Time' into a single column called 'Condition_Time'
newX <- X %>% unite("Condition_Time", Condition, Time, sep = "_", 
                    remove = FALSE)

# Create a design matrix for regression, including 'Condition_Time' and 'GAT201'
design_matrix <- model.matrix(~ Condition_Time + GAT201 - 1, data = newX)

# Check the column names of the design matrix
colnames(design_matrix)

# Calculate offset values as the median 
#of log-transformed gene expression counts
offset_values <- apply(log1p(Y[,-1]), 2, median)

# Fit a negative binomial regression model for each gene
df_NBreg <- NULL
for(g in 1:nrow(Y)){
    model <- glm.nb(unlist(Y[g, -1]) ~ design_matrix + offset(offset_values))
    df_NBreg <- rbind(df_NBreg, 
                      summary(model)$coefficients["design_matrixGAT201WT",])
}

# Adjust p-values using the Benjamini-Hochberg (BH) method
df_pval_NBreg <- df_NBreg %>% as_tibble() %>% 
    mutate(Geneid = Y$Geneid,
           pval_adj = p.adjust(df_NBreg[,4], method = "BH")
    )

# Count the number of genes with adjusted p-values less than 0.01
sum(df_pval_NBreg$pval_adj < 0.01)

# Count the number of genes with raw p-values less than 0.01
sum(df_pval_NBreg$`Pr(>|z|)` < 0.01)

# Identify significant genes with adjusted p-values less than 0.01
gr_NBreg <- df_pval_NBreg %>% filter(pval_adj < 0.01) %>% pull(Geneid)

# For Particular Gene Analysis -------------------

# Example: Analyze gene at row 762
g = 762
Y[g,]$Geneid  # Display the Geneid

# Fit and summarize the negative binomial regression model for this gene
model <- glm.nb(unlist(Y[g, -1]) ~ design_matrix + offset(offset_values))
summary(model)

# Example: Analyze gene at row 2686
g = 2686
Y[g,]$Geneid  # Display the Geneid

# Fit and summarize the negative binomial regression model for this gene
model <- glm.nb(unlist(Y[g, -1]) ~ design_matrix + offset(offset_values))
summary(model)

# ------------------------------ #
# Intersection Analysis of GAT201-Affected Genes
# ------------------------------ #

# Separate data into two groups: Group A and Group B/M
x <- Y %>% dplyr::select(-Geneid) %>% 
    dplyr::select(contains("A_"), contains("a_")) %>% 
    log1p()
y <- Y %>% dplyr::select(-Geneid) %>% 
    dplyr::select(contains("B_"), contains("M_")) %>%  
    log1p()

# Using Wilcoxon Test --------------------

# Perform Wilcoxon rank-sum test for each gene
stat <- pval <- NULL
for(i in 1:nrow(Y)) {
    wt <- wilcox.test(unlist(x[i,]), unlist(y[i,]), exact = FALSE)
    stat <- c(stat, wt$statistic)
    pval <- c(pval, wt$p.value)
}

# Adjust p-values using the BH method
pval_adj <- p.adjust(pval, method = "BH")

# Identify differentially expressed genes with adjusted p-values less than 0.05
diff_exp <- rep(0, nrow(Y))
diff_exp[pval_adj < 0.05] <- 1

# Store results in a data frame
df_pval_wilcoxtest <- data.frame(Geneid = Y$Geneid, 
                                 dexp = diff_exp,
                                 stat = stat, 
                                 pval = pval, 
                                 pval_adj = pval_adj) %>% as_tibble()

# Count the number of significant genes with adjusted p-values less than 0.1
sum(df_pval_wilcoxtest$pval_adj < 0.1)

# Identify significant genes with adjusted p-values less than 0.05
gr_wilcoxtest <- df_pval_wilcoxtest %>% filter(pval_adj < 0.05) %>% pull(Geneid)

# Negative Binomial with DESeq2 --------------------

# Copy counts_all data frame and remove 'Geneid'
counts_all_noid <- counts_all %>% 
    column_to_rownames(var = "Geneid") %>%
    dplyr::select(-Chr, -Start, -End, -Strand, -Length)

# Remove rows with all zero counts
counts_1 <- counts_all_noid[rowSums(counts_all_noid) != 0,] 

# Prepare the sample metadata for DESeq2 analysis
Gat201_samplesheet <- Gat201_samplesheet %>%
    dplyr::select(-SampleID, -Media, -Timepoint) %>%
    dplyr::mutate(
        GAT201 = factor(GAT201),
        Condition_Time = factor(paste(Condition, Time, sep = "_"))
    )

# Check the alignment of column names between counts data and sample metadata
colnames(counts_1) == Gat201_samplesheet$Title

# Create a DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
    countData = counts_1,
    colData = Gat201_samplesheet,
    design = ~ Condition_Time + GAT201
)

# Perform differential expression analysis using DESeq2
dds <- DESeq(dds)
res <- results(dds)
head(res)

# Count the number of significant genes with adjusted p-values less than 0.05
sum(res$padj < 0.05)

# Identify significant genes with adjusted p-values less than 0.05
gr_DESeq <- counts_1 %>% dplyr::filter(res$padj < 0.05) %>% rownames()

# ------------------------------ #
# Venn Diagram and Overlap Analysis of Significant Genes
# ------------------------------ #

# Define the significant gene sets from different methods
gr_diff_wilcoxtest <- df_pval_wilcoxtest %>% 
    filter(pval_adj < 0.1) %>% pull(Geneid)
gr_diff_NBreg <- df_pval_NBreg %>% filter(pval_adj < 0.01) %>% pull(Geneid)
gr_diff_DESeq <- counts_1 %>% filter(res$padj < 0.05) %>% rownames()

# Calculate the number of significant genes in each method
lg_diff <- c(length(gr_diff_wilcoxtest), length(gr_diff_NBreg), 
             length(gr_diff_DESeq))
lg_diff

# Create a Venn diagram to visualize the overlap of significant 
#genes across methods
vdiff <- venn.diagram(
    x = list("Wilcox" = gr_diff_wilcoxtest, "GLM" = gr_diff_NBreg, 
             "DESeq2" = gr_diff_DESeq), 
    euler = TRUE, 
    filename = "report/figure/vdiag_diff.png", 
    imagetype = "png"
)

# Find the intersection of significant genes across all three methods
geneid_diff <- intersect(intersect(gr_diff_wilcoxtest, gr_diff_NBreg),
                         gr_diff_DESeq)

# Define the non-significant gene sets from different methods
gr_equal_wilcoxtest <- df_pval_wilcoxtest %>% 
    filter(pval_adj > 0.1) %>% pull(Geneid)
gr_equal_NBreg <- df_pval_NBreg %>% filter(pval_adj > 0.01) %>% pull(Geneid)
gr_equal_DESeq <- counts_1 %>% filter(res$padj > 0.05) %>% rownames()

# Calculate the number of non-significant genes in each method
lg_equal <- c(length(gr_equal_wilcoxtest), length(gr_equal_NBreg), 
              length(gr_equal_DESeq))
lg_equal

# Create a Venn diagram to visualize the overlap of 
#non-significant genes across methods
vequal <- venn.diagram(
    x = list("Wilcox" = gr_equal_wilcoxtest, "GLM" = gr_equal_NBreg, 
             "DESeq2" = gr_equal_DESeq), 
    euler = TRUE,  
    filename = "report/figure/vdiag_equal.png", 
    imagetype = "png"
)

# Find the intersection of non-significant genes across all three methods
geneid_equal <- intersect(intersect(gr_equal_wilcoxtest, 
                                    gr_equal_NBreg), gr_equal_DESeq)

# Create a data frame to label genes as significant (1) or non-significant (0)
geneid_sel <- data.frame(
    Geneid = c(geneid_diff, geneid_equal),
    Label = c(rep(1, length(geneid_diff)), rep(0, length(geneid_equal)))
) 

# Save the selected gene IDs and their labels to a CSV file
write.csv(geneid_sel, "data/geneid_sel.csv", row.names = FALSE)


# ------------------------------ #
# Combination of K-mers / Motifs
# ------------------------------ #

# Function to replace DNA bases with their complementary bases
Replace_bases <- function(string) {
    string <- gsub("A", "t", string)
    string <- gsub("T", "a", string)
    string <- gsub("C", "g", string)
    string <- gsub("G", "c", string)
    string <- toupper(string)
    reversed_string <- paste(rev(strsplit(string, "")[[1]]), collapse = "")
    return(reversed_string)
}

# Function to combine motifs by summing counts of complementary sequences
Combine_motifs <- function(df_mers) {
    df_mers_replaced <- df_mers
    
    # Apply the base replacement and sum complementary sequences
    for (i in 1:ncol(df_mers)) {
        original_col <- colnames(df_mers)[i]
        replaced_col <- Replace_bases(original_col)
        df_mers_replaced[, i] <- df_mers[, i] + 
            df_mers[, match(replaced_col, colnames(df_mers))]
    }
    
    # Update column names to indicate original and replaced motifs
    new_colnames <- paste(colnames(df_mers), "_", 
                          sapply(colnames(df_mers), Replace_bases), sep = "")
    colnames(df_mers_replaced) <- new_colnames
    
    # Handle columns where the motif and its 
    #complement are identical (e.g., TA_TA)
    for (colname in colnames(df_mers_replaced)) {
        parts <- strsplit(colname, "_")[[1]]
        if (parts[1] == parts[2]) {
            df_mers_replaced[[colname]] <- df_mers_replaced[[colname]] / 2
        }
    }
    
    # Remove duplicated columns
    df_mers_replaced <- df_mers_replaced[, !duplicated(t(df_mers_replaced))]
    return(df_mers_replaced)
}

# Example: Generate a random matrix of k-mer counts for testing
df_mers <- matrix(rpois(64 * 5, 5), ncol = 64)
bases <- c("A", "T", "C", "G")
all_combinations <- expand.grid(bases, bases, bases)
colnames(df_mers) <- apply(all_combinations, 1, paste, collapse = "")

# Combine motifs using the function
Combine_motifs(df_mers)

# Combine motifs for different k-mer lengths (2, 3, 4, 5)
M_2mers_comb <- bind_cols(Geneid = M_2mers$Gene, 
                          Combine_motifs(M_2mers %>% dplyr::select(-Gene)))
M_3mers_comb <- bind_cols(Geneid = M_3mers$Gene, 
                          Combine_motifs(M_3mers %>% dplyr::select(-Gene)))
M_4mers_comb <- bind_cols(Geneid = M_4mers$Gene, 
                          Combine_motifs(M_4mers %>% dplyr::select(-Gene)))
M_5mers_comb <- bind_cols(Geneid = M_5mers$Gene, 
                          Combine_motifs(M_5mers %>% dplyr::select(-Gene)))

# Combine selected genes with gene expression and k-mer/motif counts
df_countGE <- geneid_sel %>% left_join(Y, by = "Geneid")
df_2mers <- geneid_sel %>% left_join(M_2mers_comb, by = "Geneid")
df_3mers <- geneid_sel %>% left_join(M_3mers_comb, by = "Geneid")
df_4mers <- geneid_sel %>% left_join(M_4mers_comb, by = "Geneid")
df_5mers <- geneid_sel %>% left_join(M_5mers_comb, by = "Geneid")

# Prepare data for logistic regression analysis
df_reg_2mers <- gather(df_countGE, key = "Sample", 
                       value = "Count", -Geneid, -Label) %>% 
    arrange(desc(Label), Geneid) %>% 
    left_join(df_2mers, by = c("Geneid", "Label"))

df_reg_3mers <- gather(df_countGE, key = "Sample", 
                       value = "Count", -Geneid, -Label) %>% 
    arrange(desc(Label), Geneid) %>% 
    left_join(df_3mers, by = c("Geneid", "Label"))

df_reg_4mers <- gather(df_countGE, key = "Sample", 
                       value = "Count", -Geneid, -Label) %>% 
    arrange(desc(Label), Geneid) %>% 
    left_join(df_4mers, by = c("Geneid", "Label"))

df_reg_5mers <- gather(df_countGE, key = "Sample", 
                       value = "Count", -Geneid, -Label) %>% 
    arrange(desc(Label), Geneid) %>% 
    left_join(df_5mers, by = c("Geneid", "Label"))

# Combine all k-mer data into a single data frame
df_reg_allmers <- gather(df_countGE, key = "Sample", 
                         value = "Count", -Geneid, -Label) %>% 
    arrange(desc(Label), Geneid) %>% 
    left_join(df_2mers, by = c("Geneid", "Label")) %>% 
    left_join(df_3mers, by = c("Geneid", "Label")) %>% 
    left_join(df_4mers, by = c("Geneid", "Label")) %>% 
    left_join(df_5mers, by = c("Geneid", "Label"))

# Save the data frames for further analysis
save(df_reg_2mers, 
     df_reg_3mers, 
     df_reg_4mers, 
     df_reg_5mers, 
     df_reg_allmers, 
     file = "data/df_reg_kmers.rda")

# ------------------------------ #
# Motif-Based Logistic Regression of GAT201's Influence
# ------------------------------ #

# Function to perform logistic regression using LASSO (L1 regularization)
coeflogi <- function(df_reg_kmers) {
    x_kmers <- df_reg_kmers %>% dplyr::select(-c(1:4)) %>% as.matrix()
    y_kmers <- df_reg_kmers %>% dplyr::select(Label) %>% pull()
    
    # Perform cross-validated logistic regression with LASSO
    cv_fit <- cv.glmnet(x_kmers, y_kmers, alpha = 1, family = "binomial")
    
    # Extract coefficients at the optimal lambda value
    coef(cv_fit, s = "lambda.1se")
    
    # Fit the final LASSO model using the optimal lambda
    lasso_model <- glmnet(x_kmers, y_kmers, alpha = 1, 
                          lambda = cv_fit$lambda.1se, family = "binomial")
    
    # Return the coefficients of the final model
    coefficient <- coef(lasso_model)
    return(coefficient)
}

# Perform logistic regression analysis on 2-mers, 3-mers, 4-mers, and 5-mers
set.seed(5)
coeflogi(df_reg_2mers)

set.seed(5)
coeflogi(df_reg_3mers)

set.seed(5)
coeflogi(df_reg_4mers)

set.seed(5)
coeflogi(df_reg_5mers)


# ------------------------------ #
# Negative Binomial Regression Analysis of GAT201 Influence on Gene 
#Expression Using Motif Features
# ------------------------------ #

# Function to perform logistic regression using LASSO for a specific label (1 or 0)
coefmer <- function(df_reg_kmers, l){
    # Filter data for specific samples and time points (A_ or a_ at 240 minutes)
    df_reg_WT <- df_reg_kmers %>% 
        dplyr::filter((str_detect(Sample, "A_") | str_detect(Sample, "a_")) & 
                          str_detect(Sample, "_240_"))
    
    # Prepare the feature matrix (x_kmers) and response variable (y_kmers)
    x_kmers <- df_reg_WT %>% 
        dplyr::filter(Label == l) %>% 
        dplyr::select(-c(1:4)) %>% as.matrix()
    
    y_kmers <- df_reg_WT %>% 
        dplyr::filter(Label == l) %>% 
        dplyr::select(Count) %>% pull()
    
    # Estimate the dispersion parameter (theta) for the negative binomial model
    fit_nb <- glm.nb(y_kmers ~ x_kmers)
    
    # Perform cross-validated LASSO regression using the estimated theta
    cv_fit <- cv.glmnet(x_kmers, y_kmers, alpha = 1, 
                        family = negative.binomial(theta = fit_nb$theta))
    
    # Extract coefficients using the optimal lambda value
    lasso_model <- glmnet(x_kmers, y_kmers, alpha = 1, 
                          lambda = cv_fit$lambda.1se, 
                          family = negative.binomial(theta = fit_nb$theta))
    coefficient <- coef(lasso_model)
    return(coefficient)
}

# Apply the logistic regression function to 2-mers, 3-mers data
set.seed(6)
coefmer(df_reg_2mers, 1)
set.seed(5)
coefmer(df_reg_2mers, 0)

set.seed(6)
coefmer(df_reg_3mers, 1)
set.seed(6)
coefmer(df_reg_3mers, 0)

# Estimate theta (dispersion parameter) for specific data subsets
df_reg_WT <- df_reg_2mers %>% 
    dplyr::filter((str_detect(Sample, "A_") | str_detect(Sample, "a_")) & 
                      str_detect(Sample, "_240_"))

x_2mers <- df_reg_WT %>% 
    dplyr::filter(Label == 1) %>% 
    dplyr::select(-c(1:4)) %>% as.matrix()

y_2mers <- df_reg_WT %>% 
    dplyr::filter(Label == 1) %>% 
    dplyr::select(Count) %>% pull()

# Estimate theta for Label 1
fit_nb <- glm.nb(y_2mers ~ x_2mers)
theta1 <- fit_nb$theta

# Repeat for Label 0
x_2mers <- df_reg_WT %>% 
    dplyr::filter(Label == 0) %>% 
    dplyr::select(-c(1:4)) %>% as.matrix()

y_2mers <- df_reg_WT %>% 
    dplyr::filter(Label == 0) %>% 
    dplyr::select(Count) %>% pull()

fit_nb <- glm.nb(y_2mers ~ x_2mers)
theta0 <- fit_nb$theta

# Function to perform logistic regression using pre-estimated theta values
coefmer_h <- function(df_reg_kmers, l, a){
    df_reg_WT <- df_reg_kmers %>% 
        dplyr::filter((str_detect(Sample, "A_") | str_detect(Sample, "a_")) & 
                          str_detect(Sample, "_240_"))
    
    x_kmers <- df_reg_WT %>% 
        dplyr::filter(Label == l) %>% 
        dplyr::select(-c(1:4)) %>% as.matrix()
    
    y_kmers <- df_reg_WT %>% 
        dplyr::filter(Label == l) %>% 
        dplyr::select(Count) %>% pull()
    
    # Perform cross-validated LASSO regression using the provided theta value
    cv_fit <- cv.glmnet(x_kmers, y_kmers, alpha = 1, 
                        family = negative.binomial(theta = a))
    
    # Extract coefficients using the optimal lambda value
    lasso_model <- glmnet(x_kmers, y_kmers, alpha = 1, 
                          lambda = cv_fit$lambda.1se, 
                          family = negative.binomial(theta = a))
    coefficient <- coef(lasso_model)
    return(coefficient)
}

# Apply the function to 4-mers and 5-mers data using pre-estimated theta values
set.seed(6)
coefmer_h(df_reg_4mers, 1, theta1)
set.seed(6)
coefmer_h(df_reg_4mers, 0, theta0)

set.seed(6)
m1 <- coefmer_h(df_reg_5mers, 1, theta1)
set.seed(6)
m2 <- coefmer_h(df_reg_5mers, 0, theta0)

# Convert the coefficients to dense matrices or data frames
dense_m1 <- as.matrix(m1)
dense_m2 <- as.matrix(m2)

df_m1 <- as.data.frame(dense_m1)
df_m2 <- as.data.frame(dense_m2)

# Combine the two sets of coefficients into a single data frame
combined_df <- cbind(df_m1, df_m2)

# Convert the combined data frame to a LaTeX table using the xtable package
xtable_combined <- xtable(combined_df)

# Print the LaTeX table with options for longtable format and no 
#floating environment
print(xtable_combined, 
      type = "latex", 
      tabular.environment = "longtable", 
      floating = FALSE)


# ------------------------------ #
# Multiple Motif Prediction
# ------------------------------ #

# Significant + WT --------------------

# Filter the data for specific samples (WT at 240 minutes) and Label 
#1 (significant)
df_reg_allmers_s <- df_reg_allmers %>% 
    dplyr::filter((str_detect(Sample, "A_") | str_detect(Sample, "a_")) & 
                      str_detect(Sample, "_240_")) %>% 
    dplyr::filter(Label == 1)

# Select the k-mers matrix for analysis
kmers_matrix <- df_reg_allmers_s %>% 
    dplyr::select(-c(1:4)) %>% as.matrix()

# Standardize the matrix (mean = 0, variance = 1)
kmers_matrix_standardized <- scale(kmers_matrix)

# Perform PCA analysis on the standardized matrix
pca_model <- prcomp(kmers_matrix_standardized, center = TRUE, scale. = TRUE)

# View the contribution of each principal component
summary(pca_model)

# Calculate the contribution rate of each principal component
pca_contribution <- pca_model$sdev^2 / sum(pca_model$sdev^2)

# Create a data frame to display the contribution of each principal component
contribution_df <- data.frame(Principal_Component = 1:length(pca_contribution),
                              Contribution = pca_contribution) %>% 
    dplyr::arrange(desc(Contribution))  # Sort by contribution rate

# Extract the top principal components based on their contribution
num_components <- 300  
pca_features <- pca_model$x[, 1:num_components]

# Select the top components with the highest contribution rates
top_components <- contribution_df$Principal_Component[1:num_components]

# Extract the contribution rates of these top components
top_contributions <- pca_contribution[top_components]

# Calculate the total contribution rate of the selected top components
total_contribution <- sum(top_contributions)

# Plot the cumulative explained variance by principal components
explained_variance_ratio <- pca_model$sdev^2 / sum(pca_model$sdev^2)
cumulative_explained_variance <- cumsum(explained_variance_ratio)

df <- data.frame(Principal_Component = 1:length(cumulative_explained_variance),
                 Cumulative_Explained_Variance = cumulative_explained_variance)

ggplot(df, aes(x = Principal_Component, y = Cumulative_Explained_Variance)) +
    geom_line(color = "blue") +
    geom_point(color = "blue") +
    geom_hline(yintercept = 0.9556, linetype = "dashed", color = "red") + 
    # Add a horizontal line at 0.9556
    geom_vline(xintercept = which(cumulative_explained_variance >= 0.9556)[1], 
               linetype = "dashed", color = "red") +  
    # Add a vertical line for the corresponding component number
    labs(title = "Cumulative Explained Variance by Principal Components",
         x = "Number of Principal Components",
         y = "Cumulative Explained Variance") +
    theme_minimal()

# Extract the features (k-mers) corresponding to the selected 
#principal components
pca_features <- pca_model$x[, top_components]

# Use the original count data for the LASSO regression
counts <- df_reg_allmers_s$Count

# Perform LASSO regression using the PCA features and counts
cv_fit <- cv.glmnet(pca_features, counts, alpha = 1, 
                    family = negative.binomial(theta = theta1))

# Plot the cross-validated model to visualize lambda selection
plot(cv_fit)

# Extract coefficients at the optimal lambda value
lasso_model <- glmnet(pca_features, counts, alpha = 1, 
                      lambda = cv_fit$lambda.1se, 
                      family = negative.binomial(theta = theta1))

# Extract and view the model coefficients
m3 <- coef(lasso_model)

# Extract the loadings (rotations) from the PCA model
loadings <- pca_model$rotation

all_important_kmers1 <- c()

# Identify and merge the most important k-mers for each 
#selected principal component
for (pc in top_components) {
    cat("Main k-mer features for PC", pc, ":\n")
    important_kmers <- sort(abs(loadings[, pc]), decreasing = TRUE)
    top_kmers <- names(important_kmers)[1:10]  
    # Select the top 10 most important k-mers
    all_important_kmers1 <- unique(c(all_important_kmers1, top_kmers))  
    # Merge and remove duplicates
    cat(top_kmers, "\n")
}

# Output all unique important k-mer features
cat("All unique important k-mer features:\n")
print(all_important_kmers1)



# Non-Significant + WT --------------------

# Filter the data for specific samples (WT at 240 minutes) and 
#Label 0 (non-significant)
df_reg_allmers_s <- df_reg_allmers %>% 
    dplyr::filter((str_detect(Sample, "A_") | str_detect(Sample, "a_")) & 
                      str_detect(Sample, "_240_")) %>% 
    dplyr::filter(Label == 0)

# Select the k-mers matrix for analysis
kmers_matrix <- df_reg_allmers_s %>% 
    dplyr::select(-c(1:4)) %>% as.matrix()

# Standardize the matrix (mean = 0, variance = 1)
kmers_matrix_standardized <- scale(kmers_matrix)

# Perform PCA analysis on the standardized matrix
pca_model <- prcomp(kmers_matrix_standardized, center = TRUE, scale. = TRUE)

# View the contribution of each principal component
summary(pca_model)

# Calculate the contribution rate of each principal component
pca_contribution <- pca_model$sdev^2 / sum(pca_model$sdev^2)

# Create a data frame to display the contribution of each principal component
contribution_df <- data.frame(Principal_Component = 1:length(pca_contribution),
                              Contribution = pca_contribution) %>% 
    dplyr::arrange(desc(Contribution))  # Sort by contribution rate

# Extract the top 300 principal components based on their contribution
num_components <- 300  
pca_features <- pca_model$x[, 1:num_components]

# Select the top 300 components with the highest contribution rates
top_components <- contribution_df$Principal_Component[1:num_components]

# Extract the contribution rates of these top components
top_contributions <- pca_contribution[top_components]

# Calculate the total contribution rate of the selected top components
total_contribution <- sum(top_contributions)

# Extract the features (k-mers) corresponding to the selected principal 
#components
pca_features <- pca_model$x[, top_components]

# Use the original count data for the LASSO regression
counts <- df_reg_allmers_s$Count

# Perform LASSO regression using the PCA features and counts
cv_fit <- cv.glmnet(pca_features, counts, alpha = 1, 
                    family = negative.binomial(theta = theta0))

# Plot the cross-validated model to visualize lambda selection
plot(cv_fit)

# Extract coefficients at the optimal lambda value
lasso_model <- glmnet(pca_features, counts, alpha = 1, 
                      lambda = cv_fit$lambda.1se, 
                      family = negative.binomial(theta = theta0))

# Extract and view the model coefficients
m4 <- coef(lasso_model)

# Extract the loadings (rotations) from the PCA model
loadings <- pca_model$rotation

all_important_kmers2 <- c()

# Identify and merge the most important k-mers for each 
#selected principal component
for (pc in top_components) {
    cat("Main k-mer features for PC", pc, ":\n")
    important_kmers <- sort(abs(loadings[, pc]), decreasing = TRUE)
    top_kmers <- names(important_kmers)[1:10]  
    # Select the top 10 most important k-mers
    all_important_kmers2 <- unique(c(all_important_kmers2, top_kmers))  
    # Merge and remove duplicates
    cat(top_kmers, "\n")
}

# Output all unique important k-mer features
cat("All unique important k-mer features:\n")
print(all_important_kmers2)

# ------------------------------ #
# Generate LaTeX Tables for Significant and Non-Significant k-mers
# ------------------------------ #

# Convert the coefficients of the significant (m3) and non-significant 
#(m4) models to dense matrices
dense_m3 <- as.matrix(m3)
dense_m4 <- as.matrix(m4)

# Convert the matrices to data frames
df_m3 <- as.data.frame(dense_m3)
df_m4 <- as.data.frame(dense_m4)

# Combine the two data frames, m3 and m4, into one table
combined_df <- cbind(df_m3, df_m4)

# Convert the combined data frame into a LaTeX table using the xtable package
xtable_combined <- xtable(combined_df)

# Print the LaTeX table with options for longtable format and 
#no floating environment
print(xtable_combined, 
      type = "latex", 
      tabular.environment = "longtable", 
      floating = FALSE)

# Create a matrix for the significant k-mers (all_important_kmers1) 
#and convert it to LaTeX
kmers_matrix <- matrix(all_important_kmers1, ncol = 5, byrow = TRUE)

# Generate the LaTeX code for the k-mers table
kmers_table <- xtable(kmers_matrix)

# Print the LaTeX table for significant k-mers
print(kmers_table, 
      type = "latex", 
      tabular.environment = "longtable", 
      floating = FALSE)

# Create a matrix for the non-significant k-mers (all_important_kmers2) 
#and convert it to LaTeX
kmers_matrix <- matrix(all_important_kmers2, ncol = 5, byrow = TRUE)

# Generate the LaTeX code for the k-mers table
kmers_table <- xtable(kmers_matrix)

# Print the LaTeX table for non-significant k-mers
print(kmers_table, 
      type = "latex", 
      tabular.environment = "longtable", 
      floating = FALSE)