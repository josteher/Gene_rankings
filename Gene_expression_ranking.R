library(nichenetr)
library(tidyverse)
library(RaceID)




# Helper function to perform permutations on gene expression levels
permtest <- function(df,
                     n_group1 = sc@cpart %in% my_clusters %>% names(sc@cpart)[.] %>% .[grep("WT|wt|Wt",.) ] %>% length(),
                     n_group2 = sc@cpart %in% my_clusters %>% names(sc@cpart)[.] %>% .[grep("KO|ko|Ko",.) ] %>% length(),
                     my_clusters = c(18,19,23)
){
  set.seed(42)
  # Create vector with genotypes encoded by 1 or -1  depending on the occurence in the data, 
  r_g1 <- n_group1/(n_group1 + n_group2)
  r_g2 <- n_group2/(n_group1 + n_group2)
  
  # Sample according to the occurence in the selected groups of the dataset
  p <- apply(replicate(1000, df$value  * sample( c( rep(1,  round( length(df$value) * r_g1)   ),
                                                    rep(-1, round( length(df$value) * r_g2)   )
  )   ) * (1/df$total_rel_tx) * 100), 
  2, 
  function(x){ x[order(abs(x), decreasing = T)] %>% 
      cumsum() }) %>% 
    # Calculate min an max of sampled data
    apply(., 
          2,
          function(x){ c(min=min(x), max=max(x), range=max(x) - min(x) )  }) %>% 
    # Calculate mean and sd of min- and max-values
    apply(.,
          1,
          function(x){ c(mean=mean(x), sd=sd(x))  } )
  return(p)
}


# Gene expressionranking function
gene_ranks <- function(sc_object,
                       genes,
                       my_clusters = c(18,19,23),
                       n_group1 = sc@cpart %in% my_clusters %>% names(sc@cpart)[.] %>% .[grep("WT|wt|Wt",.) ] %>% length(),
                       n_group2 = sc@cpart %in% my_clusters %>% names(sc@cpart)[.] %>% .[grep("KO|ko|Ko",.) ] %>% length(),
                       min_expr = 0,
                       rank_cutoff = 500,
                       perm_rank_cutoff = 150,
                       perm_min_cells = 100,
                       p_val_cutoff = 0.05
                       ){
  

  # Metadata 
  plotmap_df <- sc_object@cpart %>% 
    as.data.frame() %>%
    rownames_to_column() %>%
    # Bind columns since cells in sc_object@cpart and tsne should be ordered the same way
    bind_cols(sc_object@umap) %>%
    magrittr::set_colnames(c("cellid", "cluster", "X", "Y")) %>%
    mutate("is_medoid"= cellid %in% sc_object@medoids, 
           cluster = factor(cluster, levels = 1:length(unique(cluster)) ) ) %>% 
    mutate(genotype = substr(cellid, 10, 11))
  
  
  # Gene expression dataframe
  genes_df <- cbind(plotmap_df, 
                    (sc_object@ndata[genes, rownames(plotmap_df), drop = F]) %>% 
                      as.matrix() %>% t() ) %>% 
    data.table::data.table() %>%
    data.table::melt(., 
                     id.vars = c('cellid', 'cluster', 'X', 'Y', 'is_medoid', 'genotype'), 
                     measure.vars = genes)

  
  # Rank genes by expression and compare by genotype
  generank <- genes_df %>% 
    filter(cluster %in% my_clusters) %>% filter(value > min_expr) %>%
    group_by(variable) %>% 
    mutate(rank = rank(-value)) %>% 
    mutate(genotype_pm = ifelse(genotype == "WT", 1, -1 )) %>% 
    mutate(new_rank = genotype_pm * rank) %>% 
    # Important to sort according to rank before cumsum
    arrange(rank) %>% 
    mutate(cum_genotype_pm = cumsum(genotype_pm)) %>% 
    
    mutate(#n_per_gene = n(),
      total_rel_tx = sum(value)) %>% 
    # Calculate contribution of transcript per gene per cell in percent(%)
    mutate(cum_transcript_pm = cumsum( 100 * (genotype_pm * value)/total_rel_tx  )) 
  
  
  # Summarize gene rankings
  generank_summary <- generank %>% 
    filter(rank < rank_cutoff) %>% 
    # Which is max value in top 200, which is max value
    summarise(max_val = max(cum_transcript_pm),
              min_val = min(cum_transcript_pm),
              range_val = max_val - min_val,
              abs_max = max(cum_transcript_pm %>% abs() ),
              n_cells = n(),
              x = (min_val + max_val)/2,
              y = log2(abs(max_val - min_val) + 0.000001)
    ) 
  
  
  # Perform permutation test
  perm_stats <- generank %>% 
    group_by(variable) %>% 
    mutate(n_cells = n()) %>% 
    filter(n_cells > perm_min_cells) %>% 
    filter(rank < perm_rank_cutoff) %>%
    ungroup() %>% 
    mutate(variable = as.character(variable) ) %>% 
    split(f = as.factor(.$variable)) %>% 
    map(. %>% permtest(., n_group1 = n_group1, n_group2 = n_group2, my_clusters = my_clusters) )
  
  
  # Combine metrics
  perm_stats_c <- perm_stats %>% 
    map(~ as.data.frame(.) %>% rownames_to_column("statistic") ) %>% 
    map( ~pivot_longer(data = ., cols = c("min", "max", "range"), names_to = "metric") ) %>% 
    bind_rows( .id = "variable")
  
  
  # Join with gene expression rank data
  perm_res <- left_join(generank_summary, 
                        perm_stats_c %>% pivot_wider(names_from = c(statistic, metric), values_from = value), 
                        on = "variable")
  
  
  # Calculate p value and extract significant genes
  sig_genes <- perm_res %>% 
    .[complete.cases(.), ] %>% 
    mutate(z_max_val =  (max_val - mean_max)/(sd_max),
           z_min_val =  (min_val - mean_min)/(sd_min),
           pvalue_max = 2*pnorm(-abs(z_max_val)),
           pvalue_min = 2*pnorm(-abs(z_min_val)),
           padj_min = p.adjust(pvalue_min, method = "BH"), 
           padj_max = p.adjust(pvalue_max, method = 'BH') ) %>% 
    filter(padj_min < p_val_cutoff | padj_max < p_val_cutoff  )
  
  
  return(list("GeneRankings" = generank,
              "GeneRankSummary"  = generank_summary,
              "SignificantGenes" = sig_genes) )
  
}



# Rankings for entire NFkB target genes
boston <- readr::read_tsv("Boston_University_NfKB_targets.tsv", col_names = F)$X4 %>% 
  convert_human_to_mouse_symbols() %>% .[!is.na(.)]
mayanla_nfkb_pred <- jsonlite::fromJSON("maayanlab.cloud_Harmonizome_api_1.0_gene_set_NF-kappaB.json")$associations$gene$symbol %>% 
  convert_human_to_mouse_symbols() %>% .[!is.na(.)]

genes <- c(mayanla_nfkb_pred, boston,   'Irf1', "Casp8", "Casp1", "Ifih1", "Ddx58", "Dhx33", "H2afy", "Cebpa", "Spib", "Cebpb") %>% 
  unique() %>% .[. %in% rownames(sc@ndata)]


# Test run
test <- gene_ranks(sc, genes)

# test <- gene_ranks(sc, sc@genes)



# Plots
test$GeneRankings %>% 
  filter(variable %in% c("Pglyrp1", "Irf8", "E2f3", "Pycard", "Brca2", "Cebpa", "Myc", "Irf1", "Casp8", "Ikbkb", "Cdk2", "Rnf31",
                         "H2afy", "Map3k14", "Nfkbia", "Pim1", "Ier3", "Bcl2l1", "Cdkn1a", "Apoe", "Naip2") ) %>% 
  ggplot(. , aes(y = -rank , x = cum_transcript_pm , color = genotype)) +
  geom_point(size = 0.1)+  
  scale_color_manual(name = "genotype", values = c('darkred', 'royalblue')) + 
  facet_wrap(~variable, scales = 'free') + 
  geom_vline(aes(xintercept=0), colour="black", lwd=0.2)+
  theme_minimal() + 
  theme(strip.text.x = element_text(size = 14, face = c('italic') ))



test$GeneRankSummary %>%  
  filter(n_cells > 100) %>% 
  # filter(x < -5) %>%
  ggplot(aes (x = x, y = y , label = variable)) +
  geom_point(alpha = 0.5, color ="grey") + 
  
  # geom_text(size = 2.5) +
  # Label only certain genes
  geom_text_repel(data = . %>% filter(variable %in% ( test$SignificantGenes %>% select(variable) %>% flatten_chr() ) ),
                  aes(x = x, y = y, label = variable),
                  size = 1.5, max.overlaps = 30) +
  # 
  # geom_text_repel(data = . %>% filter(y > 4),
  #                 aes(x = x, y = y, label = variable),
  #                 size = 1.5, max.overlaps = 30) +
  
  geom_vline(aes(xintercept=0), colour="black", lwd=0.2) +
  geom_hline(aes(yintercept=0), colour="black", lwd=0.2) +
  xlab("(min_val + max_val)/2") +
  ylab("log2(abs(max_val - min_val) + 0.000001)") +
  theme_minimal()




