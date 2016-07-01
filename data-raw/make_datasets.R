#!/usr/bin/env Rscript
opt <- list(p=9,
            nsubj=1000,
            rho=0.5,
            output_folder = "../data/",
            seed=123)

set.seed(opt$seed)

##*************************************************************************##
## Null dataset
##*************************************************************************##
p_null <- opt$p
nsubj_null <- opt$nsubj

n_null <- 0
l_null <- 1
a_null <- 15
b_null <- 15
m_null <- rnorm(p_null, n_null, sqrt(l_null))
s_null <- rgamma(p_null, a_null, b_null)
S_null <- diag(s_null^2)
W_null <- cov2cor(S_null)
WChol_null <- t(chol(W_null))

counts_null      <- exp(mvtnorm::rmvnorm(nsubj_null, m_null, S_null))
totals_null      <- rowSums(counts_null)
compositions_null <- counts_null/totals_null

colnames(compositions_null) <- paste0("f_n_", seq_len(p_null))
colnames(counts_null)       <- paste0("f_n_", seq_len(p_null))

##*************************************************************************##
## Hard null
##*************************************************************************##
p_hard_null     <- opt$p
nsubj_hard_null <- opt$nsubj
p_dominant      <- 2

n_hard_null <- 3
l_hard_null <- 2
m_hard_null <- sort(rnorm(p_hard_null, n_hard_null, sqrt(l_hard_null)))
a_hard_null <- 1.5
b_hard_null <- 0.5
s_hard_null <- sort(rgamma(p_hard_null, a_hard_null, b_hard_null))
S_hard_null <- diag(s_hard_null^2)
W_hard_null <- cov2cor(S_hard_null)
WChol_hard_null <- t(chol(W_hard_null))

counts_hard_null <- exp(mvtnorm::rmvnorm(nsubj_hard_null, m_hard_null,
                                         S_hard_null))
totals_hard_null <- rowSums(counts_hard_null)
compositions_hard_null <- counts_hard_null/totals_hard_null

colnames(compositions_hard_null) <- paste0("f_hn_", seq_len(p_hard_null))
colnames(counts_hard_null)       <- paste0("f_hn_", seq_len(p_hard_null))

##*************************************************************************##
## Positive spike-in dataset
##*************************************************************************##
p_pos_spike     <- opt$p
nsubj_pos_spike <- opt$nsubj

n_pos_spike <- 2.5
l_pos_spike <- l_null
m_pos_spike <- rnorm(p_pos_spike, n_pos_spike, sqrt(l_pos_spike))
a_pos_spike <- a_null
b_pos_spike <- b_null
s_pos_spike <- rgamma(p_pos_spike, a_pos_spike, b_pos_spike)

W_pos_spike <- diag(p_pos_spike)
i <- p_pos_spike
k <- p_pos_spike - 1
W_pos_spike[i, k] <- as.numeric(opt$rho)
W_pos_spike[k, i] <- W_pos_spike[i, k]
S_pos_spike <- diag(s_pos_spike) %*% W_pos_spike %*% diag(s_pos_spike)

counts_pos_spike      <- exp(mvtnorm::rmvnorm(nsubj_pos_spike, m_pos_spike,
                                              S_pos_spike))
totals_pos_spike      <- rowSums(counts_pos_spike)
compositions_pos_spike <- counts_pos_spike/totals_pos_spike

colnames(compositions_pos_spike) <- paste0("f_ps_", seq_len(p_pos_spike))
colnames(counts_pos_spike)       <- paste0("f_ps_", seq_len(p_pos_spike))

##*************************************************************************##
## Negative spike-in dataset
##*************************************************************************##
p_neg_spike     <- opt$p
nsubj_neg_spike <- opt$nsubj

n_neg_spike <- 2.5
l_neg_spike <- l_null
m_neg_spike <- rnorm(p_neg_spike, n_neg_spike, sqrt(l_neg_spike))
a_neg_spike <- a_null
b_neg_spike <- b_null
s_neg_spike <- rgamma(p_neg_spike, a_neg_spike, b_neg_spike)

W_neg_spike <- diag(p_neg_spike)
i <- p_neg_spike
k <- p_neg_spike - 1
W_neg_spike[i, k] <- -as.numeric(opt$rho)
W_neg_spike[k, i] <- W_neg_spike[i, k]
S_neg_spike <- diag(s_neg_spike) %*% W_neg_spike %*% diag(s_neg_spike)

counts_neg_spike      <- exp(mvtnorm::rmvnorm(nsubj_neg_spike, m_neg_spike, S_neg_spike))
totals_neg_spike      <- rowSums(counts_neg_spike)
compositions_neg_spike <- counts_neg_spike/totals_neg_spike

colnames(compositions_neg_spike) <- paste0("f_ns_", seq_len(p_neg_spike))
colnames(counts_neg_spike)       <- paste0("f_ns_", seq_len(p_neg_spike))

devtools::use_data(counts_null,      compositions_null,
                   counts_hard_null, compositions_hard_null,
                   counts_pos_spike, compositions_pos_spike,
                   counts_neg_spike, compositions_neg_spike)
