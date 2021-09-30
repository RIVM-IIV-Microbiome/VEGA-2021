

library(vegan)
library(tibble)
library(dplyr)
library(ggrepel)

cols_grps <- c(
  `Meat consumer` = "#ef476f",
  Pescatarian = "#f2cc8f",
  Vegan = "#4ecdc4",
  Vegetarian = "#1a936f"
)
set.seed(183922)

mpa.abund <- t(abundances(mpa.ctabund.ps, "hellinger"))
mpa.sam <- meta(mpa.ctabund.ps)


dbRDA.mpa <- rda(mpa.abund ~ Age_groups + Group + Gender + ESBL_status,
                 permutations = 999,distance = "bray", mpa.sam)

permutest(dbRDA.mpa, permutations = 999)

dbRDA.mpa.sum <- summary(dbRDA.mpa)
meta.vecs.mpa <- cbind.data.frame(center = 0, dbRDA.mpa.sum$biplot[,c(1:2)])

var.explain1 <- round(eigenvals(dbRDA.mpa)[1]/sum(eigenvals(dbRDA.mpa))*100,2)
var.explain2 <- round(eigenvals(dbRDA.mpa)[2]/sum(eigenvals(dbRDA.mpa))*100,2)

spec.scores.mpa <- cbind.data.frame(center = 0, scores(dbRDA.mpa)$species)
site.scores.mpa <- as.data.frame(scores(dbRDA.mpa)$sites)

mpa.site.scores <- site.scores.mpa %>%
  rownames_to_column(var = "MG_Sample_ID")

mpa.rda <- mpa.site.scores %>%
  inner_join(mpa.sam, by="MG_Sample_ID")

tax.mus <- getTaxaTibble(mpa.ctabund.ps)

spec.scores.mpa <- spec.scores.mpa %>%
  filter(abs(RDA1) > 0.1 | abs(RDA2) >0.1) %>%
  rownames_to_column("FeatureID") %>%
  inner_join(tax.mus, by="FeatureID")

#clean names
spec.scores.mpa$FeatureID <- gsub("s__","",spec.scores.mpa$FeatureID)

scale.arrows.sp <- 5
scale.arrows.vars <- 1
rda.plot <- ggplot() +
  # Add site scores as points
  # geom_point(data = site.scores, aes(x = RDA1, y = RDA2)) +
  geom_point(data = mpa.rda,
             aes(x = RDA1, y = RDA2,
                 fill = Group),
             size=3,
             shape=21,
             alpha=0.7) +
  scale_fill_manual("",values=cols_grps) +
  geom_hline(yintercept = 0, color="grey70")+
  geom_vline(xintercept = 0,color="grey70") +
  geom_segment(data = spec.scores.mpa,
               aes(x = center, y = center,
                   xend = scale.arrows.sp*RDA1, yend = scale.arrows.sp*RDA2),
               alpha = 0.5, color = "black",
               arrow = arrow(angle = 20, length = unit(.1, "inches"), type = "open"),
               show.legend = F) +
  geom_text_repel(data = spec.scores.mpa,
                  aes(x = scale.arrows.sp*RDA1,
                      y = scale.arrows.sp*RDA2),
                  label = spec.scores.mpa$FeatureID,
                  label.padding = .3,
                  alpha = 0.7,
                  fontface="italic",
                  size=3) +
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
  geom_text_repel(data = meta.vecs.mpa,
                  aes(x = scale.arrows.vars*RDA1,
                      y = scale.arrows.vars*RDA2),
                  label = rownames(meta.vecs.mpa),
                  color = "blue",
                  segment.alpha = 0.3,
                  alpha = 0.7,
                  size=3) +
  xlab(paste0("dbRDA1 [",var.explain1,"%]")) +
  ylab(paste0("dbRDA2 [",var.explain2,"%]")) +
  theme_vega()

rda.plot

p.cop <- plot_density(mpa.reabund.ps, "s__Prevotella_copri", log10 = TRUE) + theme_vega()

pcopri <- abundances(mpa.reabund.ps)["s__Prevotella_copri",]/100
rda2 <- mpa.rda$RDA2

plot(log10(pcopri), log10(rda2))

rda.plot | (plot(log10(pcopri), log10(rda2)))

plot_listed_taxa(mpa.ctabund.ps,
                 select.taxa = "s__Collinsella_aerofaciens",
                 group = "ESBL_status",
                 group.colors=c("steelblue", "brown3"))

comps <- make_pairs(meta(mpa.ctabund.ps)$Group)
plot_listed_taxa(mpa.ctabund.ps,
                 select.taxa = c("s__Bacteroides_dorei", "s__Prevotella_copri",
                                 "s__Ruminococcus_bromii", "s__Collinsella_aerofaciens"),
                 group = "Group",
                 group.colors=cols_grps,
                 ncol = 2,
                 nrow=2) +
  ggpubr::stat_compare_means(comparisons = comps)
