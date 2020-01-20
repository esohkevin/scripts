#!/usr/bin/Rscript
#install.packages("colorspace")
library("colorspace")

#------ Examples
hcl_palettes(plot = T)
q4 <- qualitative_hcl(4, "Dark 3")
q4
sequential_hcl(4)
hclplot(qualitative_hcl(9, c=50, l=70))
diverge_hcl(9)
diverging_hcl(7, palette = "Tropic", h2 = 0, register = "mytropic")
demoplot(diverging_hcl(11, "mytropic"), type = "map")
m <- rnorm(500, mean = 4, sd = 1)
n <- 500
pcol <- diverge_hcl(n, h = 170, c = 39, l = 60, fixup = F, alpha = 0.4, rev = F)
p <- hist(m, col = pcol)
p <- density(m)
plot(p, col = pcol)
polygon(p, col = pcol)

#------ Exercise 1
setwd("~/esohdata/GWAS/popstruct/maf/")

library(data.table)

maf <- fread("freqs_camgwas.txt")
colnames(maf) <- c("CHR","BP","RAF","AAF")
head(maf)
attach(maf)
n <- 50
pcol1 <- sequential_hcl(n, h = c(0, 360 * (n - 1)/n), c = 50, l = c(30, 90), alpha = 0.9)
pcol2 <- sequential_hcl(n, h = c(0, 360 * (n - 1)/n), c = 50, l = c(30, 90), rev = TRUE, alpha = 0.9)
p1 <- hist(AAF, breaks = 100)
p2 <- hist(AAF, breaks = 100)
plot(p1, col = pcol1 )
plot(p2, col = pcol2, add = T)
d <- density(AAF)
plot(d, col = pcol1)
polygon(d, col = pcol1)
dev.off()

a <- "a.b"

a#------- Exercise 2
setwd("~/esohdata/GWAS/popstruct/admixture/")
dir()

adm_plot <- function(q = q_estimates_file, k = number_of_populations) 
  {
   for (pop in k) {
    n <- k
    admcol <- qualitative_hcl(n, h = 95, c = 150, l = 45)
    adm_base <- basename(q)
    out_name <- gsub(".Q", ".png", adm_base)
    adm <- fread(adm_base)
    tdl <- as.matrix(adm)
    png(out_name, height = 13, width = 20, units = "cm", res = 100, points = 14)
    bp <- barplot(t(tdl), col = admcol, space = 0, ylab = "Ancestral proportions", 
                  xlab = "Individual #", border = NA)
    dev.off()
  }
}

adm_plot(q = "adm-data.2.Q", k = 2)
dev.off()
