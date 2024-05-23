# Plot stock indicators as four panels
plot_indicators <- function(dat) {
  # Catch and TAC (panel A)
  p_catch <- ggplot(data = dat, mapping = aes(x = Year, y = Catch)) +
    geom_line() +
    geom_hline(yintercept = dat$TAC, linetype = "dashed") +
    scale_y_continuous(labels = comma) +
    labs(x = NULL, y = "Catch (1,000 t)") +
    expand_limits(y = 0) +
    annotate(
      geom = "text", x = Inf, y = Inf, label = "(A)", vjust = 2, hjust = 2
    ) +
    theme(axis.text.x = element_blank())
  # SSB, USR, and LRP (panel B)
  p_biomass <- ggplot(data = dat, mapping = aes(x = Year, y = SSB_med)) +
    geom_ribbon(mapping = aes(ymin = SSB_min, ymax = SSB_max), fill = "grey") +
    geom_line() +
    geom_hline(yintercept = dat$LRP, linetype = "dashed", colour = "red") +
    geom_hline(yintercept = dat$USR, linetype = "dotted", colour = "green") +
    scale_y_continuous(labels = comma) +
    labs(x = NULL, y = "Biomass (1,000 t)") +
    expand_limits(y = 0) +
    annotate(
      geom = "text", x = Inf, y = Inf, label = "(B)", vjust = 2, hjust = 2
    ) +
    theme(axis.text.x = element_blank())
  # F and M (panel C)
  p_mortality <- ggplot(data = dat, mapping = aes(x = Year, y = F_med)) +
    geom_ribbon(mapping = aes(ymin = F_min, ymax = F_max), fill = "grey") +
    geom_line() +
    geom_hline(yintercept = dat$F_lim, linetype = "dashed", colour = "red") +
    geom_hline(yintercept = dat$M_med, linetype = "dotted") +
    labs(y = "Mortality (/yr)") +
    expand_limits(y = 0) +
    annotate(
      geom = "text", x = Inf, y = Inf, label = "(C)", vjust = 2, hjust = 2
    )
  # Recruitment (panel D)
  p_recruitment <- ggplot(data = dat, mapping = aes(x = Year, y = R_med)) +
    geom_ribbon(mapping = aes(ymin = R_min, ymax = R_max), fill = "grey") +
    geom_line() +
    scale_y_continuous(labels = comma) +
    labs(y = "Recruitment (million)") + # 1,000
    expand_limits(y = 0) +
    annotate(
      geom = "text", x = Inf, y = Inf, label = "(D)", vjust = 2, hjust = 2
    )
  # Four-panel plot
  p <- p_catch + p_biomass + p_mortality + p_recruitment
  p
}
