#This script is for figure 1 of the manuscript. 

# Load required package
if (!require("VennDiagram")) install.packages("VennDiagram")
library(VennDiagram)

# Define gene counts
dc3000_unique <- 1428
p25c2_unique <- 1197
shared_genes <- 3995

# Create Venn diagram
venn.plot <- draw.pairwise.venn(
  area1 = dc3000_unique + shared_genes,
  area2 = p25c2_unique + shared_genes,
  cross.area = shared_genes,
  category = c("DC3000", "p25.c2"),
  fill = c("#0072B2", "#E69F00"), 
  alpha = 0.6,
  cat.pos = c(310, 50),
  cat.dist = 0.05,
  scaled = TRUE,
  lty = "blank",
  cex = 1.5,
  cat.cex = 1.4,
  fontfamily = "sans",
  cat.fontfamily = "sans"
)

# Save as PDF
pdf("venn_diagram_halfcolumn.pdf", width = 3.5, height = 3.5, family = "Arial")
grid.draw(venn.plot)
dev.off()
grid.draw(venn.plot)
dev.off()