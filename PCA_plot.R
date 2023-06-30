library(dplyr)
library(ggfortify)
library(textshape)
library(mgsub)
data = read.csv("path/to/C_Pv_trinity_inoc_ComBat_seq_vst.txt")
sample_info = read.csv("path/to/C_info.csv")
sample_info = sample_info %>% column_to_rownames("Sample")
sample_info = sample_info[colnames(data[-1]),]
sample_info$Time_point = mgsub(sample_info$Time_point, c("1", "2", "3"), c("24h", "48h", "72h"))
sample_info$Treatment = mgsub(sample_info$Treatment, c("0","1", "2", "3"), c("Control", "Primed", "Primed + Treated", "Treated"))
pca_matrix <- data %>% column_to_rownames("X") %>% as.matrix() %>% t()
sample_pca <- prcomp(pca_matrix)
p <- autoplot(sample_pca, colour = sample_info$Time_point, shape = sample_info$Treatment)
x_lab <- p$labels$x
y_lab <- p$labels$y
a <- ggplot(as.data.frame(sample_pca$x), aes(x = PC1, y = PC2, colour = factor(sample_info$Time_point), shape = factor(sample_info$Treatment))) +
geom_point(size=3) + xlab(x_lab) + ylab(y_lab) + theme_light()
a + scale_shape_discrete(name ="Treatment",labels = c("C","P", "PT", "T")) + scale_colour_discrete(name = "Time point", labels = c("24h", "48h", "72h"))