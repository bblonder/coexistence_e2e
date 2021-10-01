plot_obs_pred <- function(obs, pred, frac, limits)
{
  z = (as.matrix(obs) - as.matrix(pred)) / mean(as.numeric(as.matrix(obs)))
  row.names(z) = 1:nrow(z)
  print(min(z))
  print(max(z))
  #z = z[1:100,]
  
  z = z %>%
    as_tibble %>%
    rowid_to_column(var="X") %>%
    gather(key="Y", value="Z", -1) %>%
    mutate(Y = as.numeric(factor(Y)))
  
  g = ggplot(z, aes(X, Y, fill= Z)) + 
    theme_minimal() +
    theme(panel.border = element_rect(colour = "black", size = 1, fill=NA)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.position='bottom') +
    geom_tile() +
    scale_fill_gradient2(low='red',mid='white',high='blue',midpoint = 0, name='(Obs-Pred)/Obs',limits=limits) +
    xlab('Assemblage') +
    ylab('Species') +
    ggtitle(bquote(beta ~ "=" ~ .(frac)))
  
  return(g)
}



limits_this = c(-50,50)
g_annualplant_heatmap = ggarrange(
  plot_obs_pred(data_abund_annual_0.001_obs, data_abund_annual_0.001_pred, 0.0014, limits=limits_this),
  plot_obs_pred(data_abund_annual_0.01_obs, data_abund_annual_0.01_pred, 0.013, limits=limits_this),
  plot_obs_pred(data_abund_annual_0.05_obs, data_abund_annual_0.05_pred, 0.054, limits=limits_this),
  nrow=1,ncol=3,
  common.legend=TRUE,
  legend='bottom')





g_annualplant_hexbin = ggarrange(
  plot_obs_pred_scatter(data_abund_annual_0.001_obs, data_abund_annual_0.001_pred, 0.0014),
  plot_obs_pred_scatter(data_abund_annual_0.01_obs, data_abund_annual_0.01_pred, 0.013),
  plot_obs_pred_scatter(data_abund_annual_0.05_obs, data_abund_annual_0.05_pred, 0.054),
  nrow=1,ncol=3,
  common.legend=TRUE,
  legend='bottom')


ggsave(ggarrange(g_annualplant_heatmap, g_annualplant_hexbin,
                 nrow=2,ncol=1,labels='auto'), 
       file='outputs_figures/g_annualplant_abundance.png',
       width=10,height=10,
       bg = 'white')
