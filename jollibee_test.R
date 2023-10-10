#setwd("C:/Users/DAYTADB/Documents/Dominic Thesis Concerns/semipartm_thesis")
setwd("C:/Users/domda/Documents/semipartm_thesis")

################ Load codes for specialized functions
source("Source/fac-2.R") #matrix operations
source("Source/gen_dtm.R") #generate data
source("Source/lsa.R") #latent semantic analysis
source("Source/plsa.R") #probabilistic LSA
source("Source/lda.R") #latent dirichlet allocation
source("Source/semipartm.R") #semiparametric topic modelling

library(tidytext)
library(dplyr)
library(lattice)
library(ggplot2)

################ Read data and transform to tidy matrix
feedbackData<-read.csv("Notebooks/Data/feedback_text.csv")
feedbackText<-tibble(feedbackData[,c(3,14,13)])
colnames(feedbackText)<-c("feedback_no","line_numer","feedback_text")

feedback_tokens<-unnest_tokens(feedbackText,word,feedback_text)

# get word counts
corpusTokensCount<-count(feedback_tokens,word, sort = TRUE)
corpusMostFrequentTokens<-head(corpusTokensCount,800)
corpusLeastFrequentTokens<-tail(corpusTokensCount,800)


corpusTokensCount %>%
  count(word, sort = TRUE) %>%
  filter(n > 600) %>%
  mutate(word = reorder(word, n)) %>%
  ggplot(aes(word, n)) +
  geom_col() +
  xlab(NULL) +
  coord_flip()


# cast to document-term-matrix
feedbackTokensCounts<-feedback_tokens %>%
  group_by(feedback_no,word) %>%
  dplyr::count(feedback_no,word,sort=TRUE)
feedbackDTM<-cast_dtm(feedbackTokensCounts,feedback_no,word,n)
feedbackMatrix<-as.matrix(feedbackDTM)


################ Apply topic models
lsa<-LSA(Y=feedbackMatrix, transform="none")
plsa<-PLSA(Y=feedbackMatrix,
           burn_iter=100,max_iter=1000,set_seed=NULL, ntopics=10, tol=1e-3,chi=10)
lda<-LDA2(feedbackMatrix,set_seed=NULL,ntopics=10)
sptm.cv<-SPTM(feedbackMatrix,dtm$Z,epsilon="cv",ntopics=10,reg_type="splines",
              tol=1e-3,max_iter=1000,burn_iter=100,set_seed=923)