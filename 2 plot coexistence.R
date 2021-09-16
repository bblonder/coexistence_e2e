do_plots <- function(results_table, fn)
{
  # MAKE PLOTS
  g_richness = ggplot(results_table, aes(x=factor(frac),y=richness.r2,col=method)) +
    geom_boxplot() +
    theme_bw() +
    ylim(0,1) +
    ylab("R2 for predicted richness") +
    xlab("Fraction of assemblages in training")
  ggsave(g_richness, file=sprintf('outputs_statistical/g_%s_richness.pdf',fn),width=8,height=6)
  
  g_fs = ggplot(results_table, aes(x=factor(frac),y=feasible.and.stable.balanced_accuracy,col=method)) +
    geom_boxplot() +
    theme_bw() +
    ylim(0,1) +
    ylab("Balanced accuracy of feasible/stable prediction") +
    xlab("Fraction of assemblages in training")
  ggsave(g_fs, file=sprintf('outputs_statistical/g_%s_feasible_stable.pdf',fn),width=8,height=6)
  
  
  g_composition = ggplot(results_table, aes(x=factor(frac),y=composition.balanced_accuracy,col=method)) +
    geom_boxplot() +
    theme_bw() +
    ylim(0,1) +
    ylab("Balanced accuracy of composition prediction") +
    xlab("Fraction of assemblages in training")
  ggsave(g_composition, file=sprintf('outputs_statistical/g_%s_composition.pdf',fn),width=8,height=6)
  
  g_abundance = ggplot(results_table, aes(x=factor(frac),y=abundance.r2,col=method)) +
    geom_boxplot() +
    theme_bw() +
    ylim(0,1) +
    ylab("R2 of abundance prediction") +
    xlab("Fraction of assemblages in training")
  ggsave(g_abundance, file=sprintf('outputs_statistical/g_%s_abundance.pdf',fn),width=8,height=6)
}