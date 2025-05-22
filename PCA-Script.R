setwd("C:/Users/khosravi/OneDrive/Desktop/MetaAnalysis_khosravi/PCA/PCA-mergedStudies/mergedStudies")


####################
#packages
####################
install.packages("ggplot2")
install.packages("plotly")
install.packages("dplyr")
install.packages("pheatmap")
install.packages("corrplot")  
library(corrplot)
library(ggplot2)
library(plotly)
library(dplyr)
library(pheatmap)
####################
#GSE95400 preparation
####################
GSE95400 <- read.delim("GSE95400_series_matrix.txt", sep = "\t", header = TRUE)

rownames(GSE95400) <- GSE95400$ID_REF
GSE95400 <- GSE95400[,-1]
View(GSE95400)
#transpose
GSE95400_t <- t(GSE95400)
View(GSE95400_t)
#make the info table
GSE95400_t_info <- data.frame(
  SampleId=c("GSM2509848", "GSM2509849","GSM2509850", "GSM2509851"),
  Group = c(rep("Treated",time=2),rep("Control",time=2)),
  Study = c(rep("GSE95400",time=4))
 
)
View(GSE95400_t_info)



######################
#GSE115341 preparation
######################
GSE115341 <- read.delim("GSE115341_series_matrix.txt", sep = "\t", header = TRUE)
rownames(GSE115341) <- GSE115341$ID_REF
GSE115341 <- GSE115341[,-1]

# Transpose
GSE115341_t <- t(GSE115341)
View(GSE115341_t)
#make the info table 
GSE115341_t_info <- data.frame(
  SampleId = c("GSM3175736", "GSM3175737", "GSM3175738", "GSM3175739", "GSM3175740", "GSM3175741"),
  Group = c(rep("Control", time = 3), rep("Treated", time = 3)),
  Study = c(rep("GSE115341", time = 6))
)

#################
#PCA Preparation
#################
common_genes <- intersect(rownames(GSE95400), rownames(GSE115341))

# Subset both datasets to common genes
GSE95400_common <- GSE95400[common_genes, ]
GSE115341_common <- GSE115341[common_genes, ]

# Transpose for PCA (samples as rows)
GSE95400_t_common <- t(GSE95400_common)
GSE115341_t_common <- t(GSE115341_common)

# Combine data
combined_data <- rbind(GSE95400_t_common, GSE115341_t_common)

# Combine info
combined_info <- rbind(GSE95400_t_info, GSE115341_t_info)


####################
# CLEANING: remove NA or zero-variance genes
####################
clean_data <- combined_data[, colSums(is.na(combined_data)) == 0]
variances <- apply(clean_data, 2, var)
clean_data <- clean_data[, variances > 0]

#do PCA #error
pca_res <- prcomp(clean_data, center = TRUE, scale. = TRUE)

pca_df <- as.data.frame(pca_res$x)
pca_df$SampleId <- rownames(pca_df)

# add grup and study info to pca by sampleId
pca_df <- left_join(pca_df, combined_info, by = "SampleId")

####################
# 2d_plot
####################
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, shape = Study, label = SampleId)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, size = 3) +
  labs(title = "PCA 2D Plot", x = "PC1", y = "PC2") +
  theme_minimal()

print(p)

####################
# 3d_plot
####################
p3 <- plot_ly(pca_df, x = ~PC1, y = ~PC2, z = ~PC3,
              color = ~Group,
              symbol = ~Study,
              text = ~paste("Sample:", SampleId,
                            "<br>Group:", Group,
                            "<br>Study:", Study),
              type = "scatter3d", mode = "markers") %>%
  layout(title = "PCA 3D Plot",
         scene = list(
           xaxis = list(title = 'PC1'),
           yaxis = list(title = 'PC2'),
           zaxis = list(title = 'PC3')
         ))

p3

####################
# Correlation Heatmap
####################
#corelation Matrix
cor_matrix <- cor(t(clean_data)) 

rownames(combined_info) <- combined_info$SampleId

# heatmap
pheatmap(cor_matrix,
         main = "Sample Correlation Heatmap",
         annotation_row = combined_info[, c("Group", "Study")],
         annotation_col = combined_info[, c("Group", "Study")],
         show_rownames = TRUE,
         show_colnames = TRUE,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "complete",
         fontsize_row = 8,
         fontsize_col = 8)

#correlation plot
corrplot(cor_matrix,
         method = "circle",      
         type = "upper",         
         order = "hclust",      
         tl.col = "black",       
         tl.srt = 45,            
         addCoef.col = "black") 
