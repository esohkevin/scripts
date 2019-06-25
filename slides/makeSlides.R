#!/usr/bin/Rscript

# Set working directory
setwd("/home/esoh/Git/scripts/slides")

# Load packages
require("knitr")
require("markdown")

# Create Slides
knit("test.Rmd")
system("pandoc -s -t slidy test.md -o test.html")

