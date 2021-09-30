##############

# Convert the metaphlan percent abundances to clr
transform.clr <- function(x){
  ab.x <- abundances(x)
  
  ab.x <- ab.x/100
  
  #colnames(ab.x) <- colnames(x)
  
  if (any(ab.x == 0)) {
    v <- as.vector(ab.x)
    minval <- min(v[v > 0])/2
    ab.x <- ab.x + minval
  }
  
  # Pick samples x taxa abundance matrix
  d <- t(apply(ab.x, 2, function(x) {
    log(x) - mean(log(x))
  }))
  
  if (nrow(d) == ncol(ab.x)) {
    rownames(d) <- colnames(ab.x)
    colnames(d) <- rownames(ab.x)
  } else {
    colnames(d) <- colnames(ab.x)
    rownames(d) <- rownames(ab.x)
  }
  
  ab.x <- t(d)
}

###
theme_vega <- function(){
  theme_minimal() 
}
############
readMetaphlan <- function(path="data_raw/all_humann/metaphlan_bugs",
                          pattern = "*_bugs_list.tsv",
                          type=c("relative_abundance","estimated_number_of_reads_from_the_clade")){
  
  
  mpa_files <- list.files(path=path,
                          pattern = "*_bugs_list.tsv",
                          full.names = T,
                          recursive = FALSE)
  
  #map(control_files, fread, select = c("species_id","species_abund"))
  #out <- map(control_files, fread, select = c("species_id","species_abund"))
  file_nam = gsub(path,"", mpa_files)
  file_nam = gsub("/","", file_nam)
  ans <- map2(mpa_files,
              file_nam,
              ~fread(.x,select = c("#clade_name",type[1])) %>%
                mutate(id = .y))
  #ans
  
  df <- bind_rows(ans)
  colnames(df)[1] <- "cladeName"
  
  #head(df) # check file
  df_wide <- df %>%
    reshape2::dcast(cladeName  ~ id, value.var = {{ type[1] }})
  
  colnames(df_wide) <- gsub("mg_Vega_", "VEGA", colnames(df_wide))
  
  colnames(df_wide) <- sub("\\_.*", "", colnames(df_wide))
  
  # now process the rownames to keep the lowest taxonomic classification.
  xnames <- df_wide$cladeName
  shortnames <- gsub(paste0(".+\\", "|"), "", xnames)
  rownames(df_wide) <- shortnames
  
  #######
  x2 = strsplit(xnames, split="|", fixed=TRUE)
  taxmat = matrix(NA, ncol=max(sapply(x2, length)), nrow=length(x2))
  colnames(taxmat) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")[1:ncol(taxmat)]
  rownames(taxmat) = rownames(df_wide)
  for (i in 1:nrow(taxmat)){
    taxmat[i, 1:length(x2[[i]])] <- x2[[i]]
  }
  taxmat <- gsub("[a-z]__", "", taxmat)
  taxmat <- phyloseq::tax_table(taxmat)
  
  # convert NA to zero
  df_wide[is.na(df_wide)] <- 0
  
  
  otutab <- phyloseq::otu_table(df_wide[,-1], taxa_are_rows=TRUE)
  res <- phyloseq::phyloseq(taxmat, otutab)
  
  return(res)
}


## GMM plots

plotGMMLevel1.old <- function(df,
                              level1= NULL){
  colnames(df) <- gsub(" ", "_", colnames(df))
  px <- df %>%
    filter(Level_1 %in% level1) %>%
    group_by(Group, Level_1, variable) %>%
    summarise(cpm_sum = sum(value)) %>%
    ggplot(aes(Group,cpm_sum)) +
    geom_boxplot(aes(fill=Group, color=Group), width=0.4, alpha = 0.25) +
    geom_point(aes(x = Group, color = Group),
               position = position_jitter(width = 0.15),
               size = 1, alpha = 0.25
    ) +
    stat_compare_means(comparisons = comps, aes(label = ..p.signif..)) +
    #labs(title = module_plot) +
    theme_vega() +
    rotate_x_text(angle = 45) +
    scale_fill_manual(values = cols_grps) +
    scale_color_manual(values = cols_grps) +
    ylab("Median Counts (RPKs)") +
    facet_wrap(~Level_1, scales = "free_y") +
    theme(legend.position = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) # p-labels
  return(px)
  
}

plotGMMLevel1 <- function(df,
                          level1= NULL){
  colnames(df) <- gsub(" ", "_", colnames(df))
  
  
  
  px <- df %>%
    filter(Level_1 %in% level1) %>%
    group_by(Group, Level_1, variable) %>%
    summarise(cpm_sum = sum(value)) %>% 
    ungroup()
  
  # internal for plotGMMLevel1
  do.test.plot <- function(x){
    
    stats.test <- x %>%
      group_by(Level_1) %>% 
      wilcox_test(cpm_sum ~ Group) %>%
      adjust_pvalue(method = "BH") %>%
      add_significance("p.adj") %>% 
      add_xy_position(x = "Group")
    return(stats.test)
  }
  
  wilcox.test.df <- do.test.plot(px) 
  
  plot.gm <- ggboxplot(px, x = "Group", y = "cpm_sum",
                       facet.by = "Level_1", color="Group",
                       palette = cols_grps, scales = "free_y",
                       add = "jitter",width=0.4,
                       jitter = 0.2, add.params = list(alpha=0.25)) + 
    theme_vega() + 
    stat_pvalue_manual(wilcox.test.df,
                       label = "p.adj.signif",
                       hide.ns = TRUE, 
                       tip.length = 0) +
    ylab("Median Counts (RPKs)") +
    xlab("") +
    theme(legend.position = "none",
          axis.title.x = element_blank()) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
    rotate_x_text(angle = 45) + 
    theme(legend.position = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) # p-labels
  return(plot.gm)
  
}




plotGMMLevel2 <- function(df,
                          level2= NULL){
  colnames(df) <- gsub(" ", "_", colnames(df))
  px <- df %>%
    filter(Level_2 %in% level2) %>%
    group_by(Group, Level_2, variable) %>%
    summarise(cpm_sum = sum(value)) %>% 
    ungroup()
  
  ## 
  wilx.test <- px %>%
    group_by(level2) %>% 
    wilcox_test(cpm_sum ~ Group) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>% 
    add_xy_position(x = "Group")
  
  plot.gm <- ggboxplot(px, x = "Group", y = "cpm_sum",
                       facet.by = "Level_2", color="Group",
                       palette = cols_grps, scales = "free_y",
                       add = "jitter", width=0.4,
                       jitter = 0.2, add.params = list(alpha=0.25)) + 
    theme_vega() + 
    stat_pvalue_manual(wilx.test,
                       label = "p.adj.signif",
                       hide.ns = TRUE, 
                       tip.length = 0) +
    ylab("Median Counts (RPKs)") +
    xlab("") +
    theme(legend.position = "none",
          axis.title.x = element_blank()) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
    rotate_x_text(angle = 45) + 
    theme(legend.position = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) # p-labels
  
  
  return(px)
  
}

##################################################################
####### GMM utilities
# main gmm df
getFilteredGMM <- function(gmmdf,
                           filter_level= NULL,
                           for_var = NULL){
  
  if(is.null(filter_level) | is.null(for_var)) {
    stop("check variable inputs")
  }
  
  filtered.df <- gmmdf %>%
    dplyr::filter(!!sym(filter_level) %in% for_var) %>%
    dplyr::group_by(!!sym(filter_level), Group, variable,Taxon) %>%
    # sum at all filter_level
    summarise(abundance = sum(value))
  return(filtered.df)
}


##################
plotGMMAbundance <- function(gm.df){
  
  main.df <- gm.df %>% 
    group_by(Group, variable) %>%
    summarise(sum.value = sum(abundance)) %>% 
    ungroup() 
  
  test.df <- main.df %>% 
    wilcox_test(sum.value ~ Group) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>% 
    add_xy_position(x = "Group")
  
  
  plot.gm <- ggboxplot(main.df, x = "Group", 
                       y = "sum.value",
                       color="Group",
                       palette = cols_grps, scales = "free_y",
                       add = "jitter", width=0.4,
                       jitter = 0.2, add.params = list(alpha=0.25)) + 
    theme_vega() + 
    stat_pvalue_manual(test.df,
                       label = "p.adj.signif",
                       hide.ns = TRUE, 
                       tip.length = 0) +
    ylab("Median Counts (RPKs)") +
    xlab("") +
    theme(legend.position = "none",
          axis.title.x = element_blank()) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
    rotate_x_text(angle = 45) + 
    theme(legend.position = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) # p-label
  
  return(plot.gm)
}


###################### 
# gmmdf <- aad.df
getTopContributors <- function(gmmdf,
                               top.bugs=11,
                               groups_com = unique(gmmdf$Group)){
  lst <- gmmdf %>%
    group_by(Group, Taxon) %>%
    summarise(mean=mean(abundance))
  
  #mutate(frequency = abundance / sum(abundance))
  #summarise(mean=mean(abundance))
  
  # get top for each group
  
  top.contb <- NULL
  for(g in groups_com) {
    
    lstw <- lst %>%
      dplyr::filter(Group == g) %>%
      dplyr::arrange(desc(mean))
    
    top.feature <- lstw$Taxon[1:top.bugs]
    top.contb[[g]] <- top.feature
    
  }
  
  top.ctb.uniq <- unique(unlist(top.contb, use.names = FALSE))
  return(top.ctb.uniq)
}

# sub.df <- aad.df

#############################
plotTopContributors <- function(sub.df,
                                metadata=meta_tab,
                                group_colour = cols_grps,
                                vec_colors = c("#f8fafd", "#1d3557"),
                                top.bugs=11){
  top.ctb <- getTopContributors(sub.df, groups_com = unique(sub.df$Group),top.bugs=top.bugs)
  sub.df <- sub.df %>%
    mutate(Top_taxa= factor(ifelse(Taxon %in% top.ctb, Taxon, "Other"))) %>%
    left_join(meta_tab, by= c("variable" ="sampleID")) %>%
    #filter(!is.na(Group)) %>%
    group_by(Group.x, Top_taxa) %>% 
    dplyr::summarise(sum.val = sum(abundance)) %>% 
    ungroup() %>% 
    group_by(Group.x) %>% 
    mutate(frequency = sum.val / sum(sum.val)) %>% 
    #arrange(desc(frequency)) %>% 
    #mutate(frequency = abundance / sum(abundance)) %>%
    filter(!Top_taxa %in% c("Other", "unclassified"))
  
  sub.df$Top_taxa <- gsub("_", " ",sub.df$Top_taxa)
  
  p <- ggplot(sub.df,
              aes(x = Group.x, y = Top_taxa)) +
    geom_tile(aes(fill= frequency), color="grey90") +
    scale_fill_gradientn("Relative \ncontribution",
                         colours = vec_colors,
                         #trans = "log10",
                         na.value = "white") +
    theme_vega() +
    theme(axis.title = element_blank(),
          axis.text.y = element_text(face = "italic")) +
    rotate_x_text(angle = 45)
  
  p
}

#####

plotSignificantContributors <- function(df, 
                                        fdr.method = "bonferroni",
                                        metadata=meta_tab,
                                        group_colour = cols_grps,
                                        vec_colors = c("#f8fafd", "#1d3557"),
                                        padj= 0.05,
                                        abund_cutoff = 10,
                                        prevalence.cutoff = 10/100) {
  
  wide.df <- df %>% 
    pivot_wider(names_from = variable, id_cols = Taxon, values_from = abundance) %>% 
    as.data.frame() %>% 
    column_to_rownames("Taxon")
  
  keep <- core_members(wide.df, abund_cutoff, prevalence.cutoff) 
  
  df <- df %>% filter(Taxon %in% keep)
  
  top.ctb <- gmmContributionKWTest(df, 
                                   metadata=meta_tab,
                                   fdr.method = "bonferroni") %>% 
    filter(p.adj <= padj) 
  
  sig.tax <- top.ctb %>% pull(Taxon)
  
  if (is.null(sig.tax)){
    stop("No significant contributor")
  }
  
  sub.df <- df %>%
    #mutate(Sig_taxa = factor(ifelse(Taxon %in% sig.tax, Taxon, "Other"))) %>%
    left_join(meta_tab, by= c("variable" ="sampleID")) %>%
    #filter(!is.na(Group)) %>%
    group_by(Group.x, Taxon) %>% 
    dplyr::summarise(sum.val = sum(abundance)) %>% 
    ungroup() %>% 
    group_by(Group.x) %>% 
    mutate(frequency = sum.val / sum(sum.val)) %>% 
    #arrange(desc(frequency)) %>% 
    #mutate(frequency = abundance / sum(abundance)) %>%
    filter(Taxon!="Other")
  
  sub.df <- sub.df %>% filter(Taxon %in% sig.tax)
  sub.df$Taxon <- gsub("_", " ",sub.df$Taxon)
  
  p <- ggplot(sub.df,
              aes(x = Group.x, y = Taxon)) +
    geom_tile(aes(fill= frequency), color="grey90") +
    scale_fill_gradientn("Relative \ncontribution",
                         colours = vec_colors,
                         #trans = "log10",
                         na.value = "white") +
    theme_vega() +
    theme(axis.title = element_blank(),
          axis.text.y = element_text(face = "italic")) +
    rotate_x_text(angle = 45)
  
  p
  
}


gmmContributionKWTest <- function(df, 
                                  metadata=meta_tab,
                                  fdr.method = "bonferroni"){
  sdh <- sdh.test <- NULL
  
  sdh <- df %>%
    #mutate(Top_taxa= factor(ifelse(Taxon %in% top.ctb, Taxon, "Other"))) %>%
    left_join(meta_tab, by= c("variable" ="sampleID")) %>% 
    group_by(Group.x, variable) %>% 
    mutate(frequency = abundance / sum(abundance)) %>% 
    ungroup()
  
  sdh.test <- sdh %>% 
    group_by(Taxon) %>% 
    kruskal_test(frequency ~ Group.x) %>%
    adjust_pvalue(method = fdr.method) %>%
    add_significance("p.adj")
  return(sdh.test)
}


########
getPrev <- function(x,
                    return_rank = rank_names(x),
                    return_taxa = taxa_names(x),
                    sort = TRUE, ...) {
  
  # global vars
  Taxa <- NULL
  
  if (class(x) != "phyloseq") {
    stop("Input is not an object of phyloseq class")
  }
  
  rank_nm <- .check_ranks(x, return_rank)
  rank_nm <- c("Taxa", rank_nm)
  
  tax_tib <- getTaxaTibble(x,
                           column_id = "Taxa",
                           select_cols = rank_names(x),
                           select_rows = taxa_names(x)
  ) %>%
    dplyr::select(!!!syms(rank_nm))
  
  prev_tbl <- prevalence(x) %>%
    as.data.frame() %>%
    rownames_to_column("Taxa") %>%
    rename(prevalence = ".") %>%
    left_join(tax_tib, by = "Taxa")
  
  # colnames(prev_tbl) <- c("Taxa", "prevalence", rank_nm)
  
  taxa_rt <- .check_taxa(x, return_taxa)
  
  prev_tbl <- prev_tbl %>%
    filter(Taxa %in% taxa_rt)
  
  if (sort) {
    prev_tbl <- prev_tbl %>%
      arrange(desc(prevalence))
    return(prev_tbl)
  } else {
    return(prev_tbl)
  }
}


#' @importFrom phyloseq rank_names
.check_ranks <- function(x, return_rank) {
  if (!is.null(return_rank) || any(return_rank %in% rank_names(x))) {
    return(return_rank)
  } else if (is.null(return_rank) || is.na(return_rank)) {
    return_rank <- rank_names(x)
    return(return_rank)
  } else if (!is.null(return_rank) && !any(return_rank %in% rank_names(x))) {
    stop("Please provide valid taxonomic rank names in 'return_rank' ")
  }
}

#' @importFrom phyloseq taxa_names
.check_taxa <- function(x, return_taxa) {
  if (!is.null(return_taxa) || any(return_taxa %in% taxa_names(x))) {
    return(return_taxa)
  } else if (is.null(return_taxa)) {
    return(taxa_names(x))
  } else if (!is.null(return_taxa) && !any(return_taxa %in% taxa_names(x))) {
    stop("Please provide valid names in 'return_taxa' ")
  }
}




#############################

plotdbRDA <- function(x,
                      rda_obj,
                      scale_taxa = 1,
                      scale_vars = 1,
                      group,
                      group_colors = brewer.pal(12,"Paired"),
                      taxa_arrow_color = "steelblue",
                      sample_arrow_color = "brown3") {
  
  
  
  dbRDA.mpa.sum <- summary(rda_obj)
  meta.vecs.mpa <- cbind.data.frame(center = 0, dbRDA.mpa.sum$biplot[,c(1:2)])
  
  # get axis vars
  var.explain1 <- round(eigenvals(rda_obj)[1]/sum(eigenvals(rda_obj))*100,2)
  var.explain2 <- round(eigenvals(rda_obj)[2]/sum(eigenvals(rda_obj))*100,2)
  
  # get samples scores
  site.scores.mpa <- as.data.frame(scores(rda_obj)$sites)
  
  mpa.site.scores <- site.scores.mpa %>%
    tibble::rownames_to_column(var = "sID")
  # get samples
  merge.sam <- biomeUtils::getSampleTibble(x, column_id = "sID")
  
  mpa.rda <- mpa.site.scores %>%
    dplyr::inner_join(merge.sam, by="sID")
  
  meta.vecs.mpa <- meta.vecs.mpa %>% 
    tibble::rownames_to_column(var="variable") %>% 
    dplyr::mutate(variable.plot = gsub("non_westernizedyes", "Non-Western",variable),
                  variable.plot = gsub("country", "",variable.plot))
  
  tax.tib <- biomeUtils::getTaxaTibble(ps.merge)
  
  # get species scores
  spec.scores.mpa <- cbind.data.frame(center = 0, scores(rda_obj)$species)
  spec.scores.mpa <- spec.scores.mpa %>%
    dplyr::filter(abs(RDA1) > 0.25 | abs(RDA2) >0.25) %>%
    tibble::rownames_to_column("FeatureID") %>%
    dplyr::inner_join(tax.tib, by="FeatureID")
  
  scale.arrows.sp <- scale_taxa
  scale.arrows.vars <- scale_vars
  
  rda.plot <- ggplot() +
    # Add site scores as points
    # geom_point(data = site.scores, aes(x = RDA1, y = RDA2)) +
    geom_point(data = mpa.rda,
               aes(x = RDA1, y = RDA2,
                   fill = non_westernized),
               size=1,
               shape=21,
               alpha=0.2) + 
    scale_fill_manual("Non-Westernaized",values=c("brown3", "grey50")) +
    geom_hline(yintercept = 0, color="grey70")+
    geom_vline(xintercept = 0,color="grey70") +
    geom_segment(data = spec.scores.mpa,
                 aes(x = center, y = center,
                     xend = scale.arrows.sp*RDA1, yend = scale.arrows.sp*RDA2),
                 alpha = 0.5, color = "black",
                 arrow = arrow(angle = 20, length = unit(.1, "inches"), type = "open"),
                 show.legend = F) +
    geom_label_repel(data = spec.scores.mpa,
                     aes(x = scale.arrows.sp*RDA1,
                         y = scale.arrows.sp*RDA2,
                         label = FeatureID),
                     label.padding = .3,
                     alpha = 0.7,
                     fontface="italic",
                     size=3,
                     max.overlaps = 10000) +
    # Add environmental vectors in blue
    geom_segment(data = meta.vecs.mpa,
                 aes(x = center,
                     y = center,
                     xend = scale.arrows.vars*RDA1,
                     yend = scale.arrows.vars*RDA2),
                 arrow = arrow(angle = 20,
                               length = unit(.1, "inches"),
                               type = "open"),
                 show.legend = F,
                 alpha = .5,
                 color = "blue") +
    geom_label_repel(data = meta.vecs.mpa,
                     aes(x = scale.arrows.vars*RDA1,
                         y = scale.arrows.vars*RDA2,
                         label = variable.plot),
                     color = "blue",
                     segment.alpha = 0.3,
                     alpha = 0.7,
                     size=3,
                     max.overlaps = 10000) +
    xlab(paste0("dbRDA1 [",var.explain1,"%]")) +
    ylab(paste0("dbRDA2 [",var.explain2,"%]")) +
    theme_minimal()
  return(rda.plot)
}





