
# Nucleotide parallelism --------------------------------------------------

plot_nuc_survival <- function(df, ylim, ybreaks, ncol){
  ggplot(df, aes(x=m, y=value, color=name, linetype = name)) + 
    geom_step(direction = "mid") +
    scale_x_continuous(limits = c(0.5, 3.5), 
                       breaks = c(1, 2, 3), 
                       minor_breaks = NULL) +
    scale_y_log10(breaks = {{ ybreaks }},
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
    labs(x=TeX("Nucleotide multiplicitiy, \\textit{m}"), 
         y=TeX("Total mutations $\\geq \\textit{m}$"),
         color = NULL) +
    scale_color_manual(values = c("blue", "red")) +
    scale_linetype(guide = 'none') +
    facet_wrap(~label, ncol = {{ ncol }}) +
    annotation_logticks(sides = "l", color = "grey40") +
    coord_cartesian(ylim = c(1, {{ ylim }})) +
    theme_bw() +
    theme(
      strip.placement = 'outside',
      strip.background = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom"
    )
}

# Gene parallelism --------------------------------------------------------

# Prepare different data sources ------------------------------------------

make_gene_length_df <- function(df, strainID){
  # takes a combined dataframe for all HAMBI species total length of gene-wise
  # length of nonsynonymous coding sites and returns a reduced dataframe for the
  # focal species strainID
  df %>%
    # filter on strainID
    filter(strainID == !!strainID) %>% 
    dplyr::select(locus_tag, l_i = ns_length) %>% 
    # remove NAs
    drop_na() %>% 
    # convert to data.frame
    column_to_rownames(var = "locus_tag") %>%
    data.frame()
}

make_raw_mut_df <- function(df, gene_lengths) {
  # formats and prepares a dataframe of observed pass gene hits per replicate
  # and includes every CDS for the genome even if there were no hits to the gene
  # (ie some rowSums are zero). Rownames are locus_tag and columns are the
  # individual replicates where at least 1 mutation was detected
  df %>%
    # first check to remove any noncoding mutations
    filter(locus_tag %in% rownames(gene_lengths)) %>%
    summarize(n = n(),
              .by = c("replicate", "locus_tag")) %>%
    pivot_wider(
      names_from = "replicate",
      values_from = "n",
      values_fill = 0
    ) %>%
    left_join(dplyr::select(
      rownames_to_column(gene_lengths, var = "locus_tag"),
      locus_tag
    ),
    .,
    by = join_by(locus_tag)) %>%
    mutate(across(-locus_tag, ~ replace_na(.x, 0))) %>%
    column_to_rownames(var = "locus_tag") %>%
    data.frame()
}

make_multiplicity_df <- function(df, gene_lengths) {
  # calculate gene multiplicity for each gene
  # m_i = n_i * mean(l_i)/l_i
  df %>%
    # first check to remove any noncoding mutations
    filter(locus_tag %in% rownames(gene_lengths)) %>%
    summarize(
      n_i = n(),
      n_replicate = n_distinct(replicate),
      .by = c("locus_tag")
    ) %>%
    select(locus_tag, n_i, n_replicate) %>%
    left_join(rownames_to_column(gene_lengths, var = "locus_tag"),
              .,
              by = join_by(locus_tag)) %>%
    replace_na(list(n_i = 0, n_replicate = 0)) %>%
    mutate(m_i = n_i * mean(l_i) / l_i) %>%
    column_to_rownames(var = "locus_tag") %>%
    data.frame()
}

prepare_all_input <- function(df, gene_table){
  # prepares a list with the required dataframes and summary values to compute the
  # gscore or to simulate the G score
  
  strainID        <- unique(df$strainID)
  gene_lengths    <- make_gene_length_df(gene_table, strainID)
  raw_mut_df      <- make_raw_mut_df(df, gene_lengths)
  multiplicity_df <- make_multiplicity_df(df, gene_lengths)
  
  return(lst(
    strainID = strainID,
    # number of replicates with at least on gene hit
    rep_num = ncol(raw_mut_df),
    # total gene hits in each replicate
    mut_reps = colSums(raw_mut_df),
    # total gene hits across all genes and replicates
    mut_tot = sum(colSums(raw_mut_df)),
    gene_lengths = gene_lengths,
    multiplicity_df = multiplicity_df,
    raw_mut_df = raw_mut_df
  ))
}

# G scores ----------------------------------------------------------------

G_score <- function(input_list){
  # number of non-syn sites in gene i
  l_i <- input_list$gene_lengths$l_i
  # compute the sum of non-syn sites in genome
  Ltot <- sum(l_i)
  # compute the average non-syn sites per gene
  Lavg <-  mean(l_i)
  # compute the number of genes in the genome
  n_genes <- length(l_i)
  # number of mutations in gene i across all replicates
  n_i <- input_list$multiplicity_df$n_i
  # total number of non-syn hits across genes and replicates
  Ntot <- input_list$mut_tot
  # compute the expected multiplicity value
  m_exp <- Ntot / n_genes
  # compute gene multiplicity for gene i
  m_i <- n_i * Lavg / l_i
  # compute the G score for each gene
  g_score <- n_i * log(m_i / m_exp)
  # net increase in log-likelihood compared to the null model of m_i / m_exp =1
  # normalized to the total mutations
  delta_l <- sum(g_score, na.rm=T) / Ntot
  return(delta_l)
}

G_score_prep <- function(input_list){
  # Special G score function preparing constants for use in G_score_resampled
  # so we don't have to unnecessarily recalculate all the constants for each
  # resample

  # number of non-syn sites in gene i
  l_i <- input_list$gene_lengths$l_i
  # compute the sum of non-syn sites in genome
  Ltot <- sum(l_i)
  # compute the average non-syn sites per gene
  Lavg <-  mean(l_i)
  # Null probability of non-syn hit on a gene based only on the
  # number of nonsynonymous sites per gene relative to the total number of
  # nonsynonymous sites in the genome (gene_probs = l_i/Ltot).
  gene_probs <- l_i/Ltot
  # number of genes in genome. In the simulation this only includes coding
  # sequences and so we don't need to filter out any genes
  n_genes <- length(l_i)
  # total number of non-syn hits across genes and replicates
  Ntot <- input_list$mut_tot
  # compute the expected multiplicity value
  m_exp <- Ntot / n_genes
  
  lst(l_i=l_i, Ltot=Ltot, Lavg=Lavg, n_genes=n_genes, 
      Ntot=Ntot, m_exp=m_exp, gene_probs=gene_probs)
}  

G_score_resampled <- function(df, input_list){
  # remaining G score function that uses constants from G_score_prep while only
  # calculating n_i dynamically
  
  # number of mutations in gene i summed across all replicates
  n_i <- rowSums(df)
  # compute gene multiplicity for gene i
  m_i = n_i * input_list$Lavg / input_list$l_i
  # compute the G score for each gene
  g_score = n_i * log(m_i / input_list$m_exp)
  # net increase in log-likelihood compared to the null model of m_i / m_exp =1
  # normalized to the total mutations
  delta_l = sum(g_score, na.rm=T) / input_list$Ntot
}

# G score permutation -----------------------------------------------------

mut_hit_sim <- function(mut_reps, gene_probs, resamples){
  # randomly sample (n=resamples) the observed number of nonsynonymous mutations
  # from each replicate (vector mut_reps) with probability based only on the
  # number of nonsynonymous sites per gene relative to the total number of
  # nonsynonymous sites in the genome (gene_probs = l_i/Ltot). Returns a list of
  # length mut_reps where for each replicate in mut_reps there is a n=resamples
  # column and n=gene_number row matrix
  map(mut_reps, \(x) rmultinom(n = resamples, size = x, prob = gene_probs))
}

convert_nameless_mat2_df <- function(m, pop_num, gene_names){
  # take output from mut_hit_sim function (list of sample x gene matrices) and
  # reformat it so that it returns a dataframe with locus_tag for each gene as
  # rownames and sampled replicates are columns
  names(m) <- LETTERS[seq_along(1: pop_num)]
  data.frame(m, row.names = gene_names)
}

gscore_p_value <- function(input_list, emp_g_score, samples){
  G_score_prepped <- G_score_prep(input_list)
  # randomly sample hits based on gene length and split to a convenient form
  multinomhitsamp_split <- mapply(asplit,
                                  MARGIN = 2,
                                  x = mut_hit_sim(input_list$mut_reps, 
                                                  G_score_prepped$gene_probs, 
                                                  samples))
  # for each random sample convert to a more useful dataframe form
  g_reps <- seq_along(1:samples) %>%
    map(\(x) convert_nameless_mat2_df(multinomhitsamp_split[x,], 
                                      input_list$rep_num, 
                                      rownames(input_list$gene_lengths))) %>%
    # calculate the g_score for each sample
    map_dbl(\(x) G_score_resampled(x, G_score_prepped))
  
  return(length(g_reps[g_reps > emp_g_score]) / length(g_reps))
}


gscore_p_value2 <- function(input_list, emp_g_score, samples){
  # also computes the observed distribution of hits per locus tag
  # See Fig S3 from this paper: https://doi.org/10.1038/s41559-020-1128-3
  G_score_prepped <- G_score_prep(input_list)
  # randomly sample hits based on gene length and split to a convenient form
  multinomhitsamp_split <- mapply(asplit,
                                  MARGIN = 2,
                                  x = mut_hit_sim(input_list$mut_reps, 
                                                  G_score_prepped$gene_probs, 
                                                  samples))
  # for each random sample convert to a more useful dataframe form
  multinomhitsamp_split_df <- seq_along(1:samples) %>%
    map(\(x) convert_nameless_mat2_df(multinomhitsamp_split[x,], 
                                      input_list$rep_num, 
                                      rownames(input_list$gene_lengths)))
  # get the observed distribution of "hits" i.e. non-syn mutations
  observed_hit_distrubution <- input_list$raw_mut_df %>% 
    rownames_to_column(var = "locus_tag") %>% 
    pivot_longer(-locus_tag) %>%
    filter(value > 0) %>% 
    summarize(n_reps = n_distinct(name), .by = c(locus_tag)) %>% 
    summarize(n = n(), .by = c(n_reps)) %>% 
    mutate(type = "Observed")
  
  null_hit_distribution <- multinomhitsamp_split_df %>% 
    map_dfr(\(x) pivot_longer(rownames_to_column(x, var = "locus_tag"), -locus_tag), .id = "resample") %>%
    filter(value > 0) %>% 
    summarize(n_reps = n_distinct(name), .by = c(resample, locus_tag)) %>% 
    summarize(n = n()/samples, .by = c(n_reps)) %>% 
    mutate(type = "Null")
  
  # calculate the g_score for each sample
  g_reps <- map_dbl(multinomhitsamp_split_df, \(x) G_score_resampled(x, G_score_prepped))
  
  # calculate p values
  g_score_pvalue <- length(g_reps[g_reps > emp_g_score]) / length(g_reps)
  
  return(lst(g_score_pvalue = g_score_pvalue, 
             hit_distribution = bind_rows(observed_hit_distrubution, null_hit_distribution)))
}

# Calculating gene-wise significance --------------------------------------

from_parallelism_statistics <- function(input_list) {
  # number of non-syn sites in gene i
  l_i <- input_list$gene_lengths$l_i
  # compute the sum of non-syn sites in genome
  Ltot <- sum(l_i)
  # number of mutations in gene i across all replicates
  n_i <- input_list$multiplicity_df$n_i
  # total number of non-syn hits across genes and replicates
  Ntot <- input_list$mut_tot
  # compute the fraction of total non-syn sites occupied by each gene
  ps <- l_i / Ltot
  # compute the expected number of mutations in each gene based on gene length
  n_exp <- Ntot * l_i / Ltot
  
  lst(
    "gene_names" = rownames(input_list$gene_lengths),
    "Ls" = l_i,
    "ns" = n_i,
    "Ltot" = Ltot,
    "Ntot" = Ntot,
    "ps" = ps,
    "n_exp" = n_exp
  )
}

calculate_poisson_log_survival <- function(n_exp, ns, ethresh, filter = FALSE){
  # survival function = 1-CDF
  survivals <- (1 - ppois(ns - 0.1, n_exp))
  # for the case if there are so many mutations in some gene that its survival is 0
  if (filter) {
    survivals[survivals==0] <- min(survivals[survivals>0])/1e10
  }
  logpvalues <- rep(0, each = length(survivals))
  logpvalues[survivals > ethresh] = -log(survivals[survivals > ethresh])
  logpvalues[survivals <= ethresh] = -(ns * log(ns / n_exp +
                                                  (ns == 0)) + ns - n_exp)[survivals <= ethresh]
  return(logpvalues)
}

calculate_unnormalized_survival_from_vector <- function(xs, min_x=NULL, max_x=NULL, min_p=1e-10){
  if (is.null(min_x)){
    min_x <- min(xs)-1
  }
  if (is.null(max_x)){
    max_x <- max(xs)+1
  }
  unique_xs <- unique(xs)
  unique_xs <- sort(append(unique_xs, c(min_x, max_x)))
  num_observations <- map_dbl(unique_xs, \(x) sum(xs >= x))
  
  # So that we can plot CDF, SF on log scale
  num_observations[0] <- num_observations[0] - min_p
  num_observations[1] <- num_observations[1] - min_p
  num_observations[-1] <- num_observations[-1] + min_p 
  
  return(list(unique_xs, num_observations))
}

map_pval_and_probs = function(observed_p, logpvalues, probabilities, ns_grid, nmin){
  sum(map2_dbl(logpvalues, probabilities, 
               \(x, y) sum((x >= observed_p) * (ns_grid >= nmin) * y)))
}

format_survivalcurves_2plot <- function(observed_ps, survivals, observed_pvalue_survival, 
                                        pstar, num_significant){
  df_step <- tibble(
    obs_p = observed_ps,
    Survival_Expected = survivals,
    Survival_Observed = observed_pvalue_survival
  ) %>%
    pivot_longer(-obs_p)  %>%
    mutate(name = factor(
      name,
      levels = c("Survival_Observed", "Survival_Expected"),
      labels = c("Survival Observed", "Survival Expected")
    ))
  
  df_seg <- tibble(
    x = c(min(observed_ps), pstar),
    y = c(num_significant, 0),
    xend = c(pstar, pstar),
    yend = c(num_significant, num_significant)
  )
  
  df_ann <- tibble(x = pstar, y = num_significant)
  
  # return
  lst(df_step = df_step,
      df_seg = df_seg,
      df_ann = df_ann
  )
}


make_gene_table <- function(pooled_pvalues, strainID, obs_gene_multiplicity, raw_hits, pstar){
  enframe(pooled_pvalues, name = "locus_tag", value = "neg_log10P") %>% 
    mutate(strainID = strainID, pstar = pstar) %>% 
    left_join(., rownames_to_column(obs_gene_multiplicity, var = "locus_tag"), join_by(locus_tag)) %>%
    left_join(., rownames_to_column(raw_hits, var = "locus_tag"), join_by(locus_tag)) %>%
    dplyr::rename(nonsyn_sites = l_i, 
                  observed_hits_n_i = n_i,
                  gene_multiplicity_m_i = m_i)
}

NullGeneLogpSurvivalFunction <- function(mylist, nmin=3){
  
  # Get formatted input
  fps <- from_parallelism_statistics(mylist)
  
  # Calculate log-P value for each gene
  logpvalues_par <- calculate_poisson_log_survival(fps$n_exp, fps$ns, 
                                                   ethresh=1e-30, filter = TRUE)
  names(logpvalues_par) <- fps$gene_names
  
  # Collect P values for genes with at least nmin (default = 3) mutations across
  # all populations this limits the number of low P values driven primarily by
  # gene length
  pooled_pvalues <- sort(logpvalues_par[fps$ns >= nmin])
  
  # Un-normalized poisson survival vector
  tsrv <- calculate_unnormalized_survival_from_vector(pooled_pvalues) 
  observed_ps <- pluck(tsrv, 1)
  observed_pvalue_survival <- pluck(tsrv, 2)
  
  # get required inputs for the survival analysis to estimate significance
  ns_grid <- rep(0:399)
  alpha <- 0.05
  
  logpvalues       <- map(fps$n_exp, \(x) calculate_poisson_log_survival(x, ns = ns_grid, ethresh = 1e-20, filter = FALSE))
  logprobabilities <- map(fps$n_exp, \(x) ns_grid * log(x) - lgamma(ns_grid + 1) - x)
  probabilities    <- map(logprobabilities, exp)
  survivals        <- map_dbl(observed_ps, \(x) map_pval_and_probs(x, logpvalues, probabilities, ns_grid, nmin))
  threshold_idx    <- which(survivals/observed_pvalue_survival <= alpha)[1]
  pstar            <- observed_ps[threshold_idx]
  num_significant  <- observed_pvalue_survival[threshold_idx]
  
  # Return pertinent results
  df2plot <- format_survivalcurves_2plot(observed_ps, survivals, observed_pvalue_survival, 
                                         pstar, num_significant)
  gene_table <- make_gene_table(pooled_pvalues, mylist$strainID, mylist$multiplicity_df, mylist$raw_mut_df, pstar)
  
  return(lst("df2plot" = df2plot, 
             "gene_table" = gene_table,
             "observed_ps" = observed_ps,
             "pooled_pvalues" = pooled_pvalues,
             "observed_pvalue_survival" = observed_pvalue_survival,
             "survivals" = survivals,
             "pstar" = pstar,
             "num_significant" = num_significant))
  
}

plotcurve <- function(input_list){
  ggplot() +
    geom_segment(
      data = input_list$df_seg,
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend
      ),
      linetype = "dashed",
      linewidth = 0.25
    ) +
    geom_step(data = input_list$df_step, aes(x = obs_p, y = value, color = name)) +
    annotate("point", x = input_list$df_ann$x, y = input_list$df_ann$y) +
    # annotate(
    #   "text",
    #   x = input_list$df_ann$x + 35,
    #   y = input_list$df_ann$y - 3,
    #   label = TeX("Critical P*\nthreshold = $1.6 \\times 10^{-4}$")
    # ) +
    labs(x = TeX("$-log_{10}P$"), y = "Number of genes") +
    scale_color_manual(values = c("blue", "grey")) +
    scale_y_log10(
      breaks = c(1, 10, 100),
      labels = scales::trans_format("log10", scales::math_format(10 ^.x))
    ) +
    annotation_logticks(sides = "l", color = "grey40") +
    coord_cartesian(ylim = c(1, 1000)) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      panel.grid = element_blank(),
    )
}

run_full_gene_parallelism_pipeline <- function(grouped, gene_lens, resamples){
  # grouped is a tibble that has been grouped at the desired hierarchy level
  # from the experiment. It should also contain a variable called "groupid" that
  # is numeric variable unique to each group. The tibble gene_lens contains the
  # number of non-synonymous sites found in each locus tag for all 23 HAMBI
  # species in a variable called "ns_length." the variable "resamples is the
  # number of resamplings you want to perform for the significance testing of
  # the G score
  
  # split into list of tibbles
  grouped_split <- group_split(grouped)
  
  # prepare ecah grouped input for gene parallelism analysis
  grouped_split_prepped <- grouped_split %>% 
    map(\(x) prepare_all_input(x, gene_lens), .progress="inputprep") 
  
  # Calculate empirical G scores for non-syn mutations
  emp_g_scores <- map(grouped_split_prepped, \(x) G_score(x), .progress="gscore")
  
  # Calculate significance of empirical G scores
  #g_score_pvals_lump <- map2_dbl(grouped_split_prepped, emp_g_scores, \(x, y) gscore_p_value(x, y, resamples))
  g_score_pvals_lump2 <- map2(grouped_split_prepped, emp_g_scores, \(x, y) gscore_p_value2(x, y, resamples), .progress="gscore_pval")
  
  # Estimate gene-wise p values and compare to null dependending only on total hits and gene length
  nglpsf_out <- map(grouped_split_prepped, \(x) NullGeneLogpSurvivalFunction(x, nmin = 1), .progress="LogPSurvival")
  
  # format G scores into a tibble
  output_gscore <-
    unnest(enframe(emp_g_scores, name = "groupid", value = "observed_g_score"), cols = observed_g_score) %>% 
    left_join(unnest(enframe(map(g_score_pvals_lump2, \(x) pluck(x, 1)), name = "groupid", value = "pvalue"), cols = "pvalue")) %>%
    left_join(group_keys(grouped))
  
  output_simulated_reps <- g_score_pvals_lump2 %>% 
    map(\(x) pluck(x, 2)) %>% 
    map_dfr(\(x) x, .id = "groupid") %>% 
    mutate(groupid = as.integer(groupid)) %>% 
    left_join(group_keys(grouped))
  
  output_gene_table <- map_dfr(nglpsf_out, \(x) pluck(x, "gene_table"), .id = "groupid") %>%
    mutate(groupid = as.integer(groupid)) %>% 
    left_join(group_keys(grouped))
  
  getdfl2plot <- map(nglpsf_out, \(x) pluck(x, "df2plot"))
  output_df2plot <- map_dfr(getdfl2plot, \(x) pluck(x, "df_step"), .id = "groupid") %>%
    left_join(map_dfr(getdfl2plot, \(x) pluck(x, "df_seg"), .id = "groupid")) %>%
    mutate(groupid = as.integer(groupid)) %>% 
    left_join(group_keys(grouped))

  return(lst(output_gscore = output_gscore, 
             output_simulated_reps = output_simulated_reps, 
             output_gene_table = output_gene_table,
             output_df2plot = output_df2plot))
}

plot_sims <- function(input_list, page){
  input_list$output_simulated_reps %>% 
    mutate(label = paste(strainID, "prey:", prey_history, "\npredator:", predator_history)) %>% 
    mutate(n_reps = factor(n_reps), 
           labels = label) %>% 
    ggplot() +
    geom_col(aes(x = n_reps, y = n, fill = type),
             position = position_dodge2(preserve = "single")) +
    labs(x = "Number of replicates with a hit in the same gene", y = "Number of genes with a hit") +
    scale_y_sqrt() + 
    scale_fill_brewer(palette = "Set1") +
    facet_wrap_paginate(vars(label), scales = "free_y",
                        labeller = labeller(label_wrap_gen()),
                        ncol = 4, nrow = 4, page = {{ page }}) +
    theme_bw() +
    theme(
      strip.placement = 'outside',
      strip.background = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom"
    )
}

# Functions for Jaccard index resampling ----------------------------------

prep4jaccard_resample <- function(gene_table_filt, genetable, strainID, ...){
  mygroup <- gene_table_filt %>% 
    filter(strainID == {{ strainID }}) %>%
    group_by(...) %>% 
    mutate(groupid = cur_group_id())
  
  mysplit <- mygroup %>% 
    group_split() %>% 
    map(\(x) pull(x, locus_tag))
  
  genetable_filt <- genetable %>% 
    filter(strainID == {{ strainID }}) %>% 
    dplyr::select(locus_tag, ns_length) %>% 
    # remove NAs
    drop_na() %>% 
    arrange(locus_tag)
  
  locus_tag <- pull(genetable_filt, locus_tag)
  l_i <- pull(genetable_filt, ns_length)
  gene_probs <- l_i/sum(l_i)
  
  return(
    lst(mygroup, mysplit, locus_tag, gene_probs)
  )
}

jaccard <- function(a, b) {
  intersection  <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  return (intersection/union)
}

pairsamplejaccard <- function(x, y, locus_tag, gene_probs){
  jaccard(
    sample(x = locus_tag, size = length(x), prob = gene_probs, replace = FALSE),
    sample(x = locus_tag, size = length(y), prob = gene_probs, replace = FALSE)
  )
}

resample_jaccard <- function(nresamples, prepped4jaccard_resample) {
  
  # defining inputs
  input_gene_lists <- prepped4jaccard_resample$mysplit
  genome_locus_tags <- prepped4jaccard_resample$locus_tag
  genome_gene_probs <- prepped4jaccard_resample$gene_probs

  # some housekeeping
  myindxgroups <- prepped4jaccard_resample$mygroup %>% 
    ungroup() %>% 
    distinct(groupid, prey_history, predator_history) %>% 
    mutate(groupname = paste0("bact: ", prey_history, ", predator: ", predator_history)) %>% 
    dplyr::select(groupid, groupname)
  
  # generate unique combinations
  combos <- combn(1:length(input_gene_lists), 2)
  
  # convert combo ids to variable names for pairs
  myindxgroupnames <- as.data.frame(t(combos)) %>% 
    `colnames<-`(c("groupid1", "groupid2")) %>%
    rownames_to_column(var = "treat_pair") %>%
    mutate(treat_pair = as.integer(treat_pair)) %>%
    left_join(rename(myindxgroups, groupid1 = groupid, groupname1 = groupname),
              by = join_by(groupid1)) %>%
    left_join(rename(myindxgroups, groupid2 = groupid, groupname2 = groupname),
              by = join_by(groupid2))
  
  # First inspect observed jaccard indexes for all pairs
  jc_obs <-
    map_dbl(asplit(combos, MARGIN = 2), 
            \(x) jaccard(input_gene_lists[[x[1]]], 
                         input_gene_lists[[x[2]]])
    )
  
  # get only combos with a Jaccard index > 0 (i.e. any overlapping elements)
  jc_obs_gt0 <- jc_obs[which(jc_obs > 0)]
  
  if (length(jc_obs_gt0) > 1) {
    ind2iterate <- asplit(combos[, which(jc_obs_gt0 > 0)], MARGIN = 2)
  } else {
    ind2iterate <- 1
  }
  
  # resampled jaccard index
  jc_resampled <- replicate(nresamples,
                            map_dbl(
                              ind2iterate,
                              \(x) pairsamplejaccard(input_gene_lists[[x[1]]], 
                                                     input_gene_lists[[x[2]]], 
                                                     genome_locus_tags, 
                                                     genome_gene_probs)
                            )
  )
  
  if (length(jc_obs_gt0) > 1) {
    # Vector of p values for each combination with multiple comparisons
    ps <- tibble(treat_pair = which(jc_obs > 0), p_values = rowSums(jc_resampled >= jc_obs_gt0) / rep(nresamples, length(jc_obs_gt0)))
  } else {
    # single p value for the single sample case
    ps <- tibble(treat_pair = which(jc_obs > 0), p_values = length(jc_resampled[jc_resampled >= jc_obs_gt0]) / nresamples)
  }
  left_join(ps, myindxgroupnames)
}

# hypergeometric test -----------------------------------------------------

enricher <- function(n_background, n_par, n_cog_background, n_cog_par, over=TRUE){
  # n_background: total genes with COG annotation 
  # n_par: Set of target genes of interest. In this case the total number of COG
  # annotated genes that have significant parallelism
  # n_cog_background: set of genes from background annotated with the COG of interest
  # n_cog_par: number of genes in a COG of interest from the total number of COG annotated
  # genes w/ parallelism
  if (over) {
    phyper(n_cog_par-1, n_cog_background, n_background-n_cog_background, n_par, lower.tail= FALSE)
  } else {
    phyper(n_cog_par, n_cog_background, n_background-n_cog_background, n_par, lower.tail= TRUE)
  }
}
