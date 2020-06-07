library(markovchain)
library(tidyverse)

dat <- read_csv("xu_mouse_dat_processed.csv")

#Efficiency Histogram
dat %>%
  ggplot(aes(x = avg_log2_fold_change)) +
  geom_histogram()

#Encode efficiency groups as factors
dat <- dat %>%
  mutate(group = factor(group))

dat %>%
  ggplot(aes(x = avg_log2_fold_change)) +
  geom_density(aes(fill = group), alpha = 0.5)

dat_train <- dat %>%
  sample_frac(0.8)

dat_test <- dat %>%
  setdiff(dat_train)

dat <- dat_train

#Focus on group 3, high efficiency
group_3 <- dat %>%
  filter(group == 3)

#Create matrix of nucleotides from grna
nucleotide_mat <- group_3 %>%
  pull(grna) %>%
  strsplit("") %>%
  {do.call(rbind, .)}

#Obtain position specific transition matrices
t_matrices <- list()

for(i in 1:(ncol(nucleotide_mat) - 1)){
  t_matrices[[i]] <- createSequenceMatrix(nucleotide_mat[, i:(i + 1)],
                                          possibleStates = c("A", "C", "G", "T"),
                                          toRowProbs = TRUE)
}

#extract_position_matrices()
#Given a matrix composed of sequences in rows,
#extract the position-specific transition or frequency matrices
#and return in a list. Transition versus frequency matrix determined by
#output_prob. (T)RUE returns transition, (F)ALSE returns frequency
extract_position_matrices <- function(seq_matrix, output_prob, possible_vals = c("A", "C", "G", "T")){
  
  if(class(seq_matrix) != "matrix") stop("Data must be in matrix format")
  
  t_mat_list <- list()
  
  for(i in 1:(ncol(seq_matrix) - 1)){
    
    t_mat_list[[i]] <- createSequenceMatrix(seq_matrix[, i:(i + 1)],
                                            possibleStates = possible_vals,
                                            toRowProbs = output_prob)
  }
  
  
  return(t_mat_list)
  
}

score_matrices <- t_matrices %>%
  map(~ log(.x/0.25))

#Create matrix of nucleotides from grna - test set
nucleotide_mat <- dat_test %>%
  pull(grna) %>%
  strsplit("") %>%
  {do.call(rbind, .)}

#Create list of matrices for each row
#indicating which transitions are present
freq_matrices <- apply(nucleotide_mat, 1,
                       function(x) extract_position_matrices(t(as.matrix(x)),
                                                             output_prob = FALSE,
                                                             possible_vals = c("A", "C", "G", "T")))

dat_scores <- map(freq_matrices,
                  function(x) map2(x, score_matrices, ~ .x * .y))

log_likelihood <- dat_scores %>%
  map(~ sum(unlist(.x))) %>%
  unlist()

cor(log_likelihood, dat_test$avg_log2_fold_change, method = "spearman")
