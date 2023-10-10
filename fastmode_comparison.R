library(ggplot2)

setwd("/Users/dominic/Documents/semipartm_thesis")
source("test-sptm.R")

test_sptm(nword=500,ndoc=500,spar=0.9,m=1,niter=5, fastmode = TRUE)
w_fastmode <- readRDS("/Users/dominic/Documents/semipartm_thesis/Notebooks/Data/500doc_500word_0.9spar_m1.RDS")

test_sptm(nword=500,ndoc=500,spar=0.9,m=1,niter=5, fastmode = FALSE)
wo_fastmode <- readRDS("/Users/dominic/Documents/semipartm_thesis/Notebooks/Data/500doc_500word_0.9spar_m1.RDS")

comp_times <- data.frame(
  FastMode = c("On","Off"),
  Compute = round(c(as.numeric(w_fastmode$comp.time), as.numeric(wo_fastmode$comp.time)),2)
)

saveRDS(comp_times,"Notebooks/Data/Fastmode_Comparison.RDS")

ggplot(data=comp_times, aes(x=FastMode, y=Compute)) +
  geom_bar(stat="identity", fill="steelblue") +
  geom_text(aes(label=paste0(Compute, " mins")), vjust=1.6, color="white", size=3.5) +
  theme_minimal()
