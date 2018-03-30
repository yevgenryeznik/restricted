#' A simulation example
#'
#' Reproduces simulations from the paper Ryeznik Y, Sverdlov O (2018) "A Comparative Study of Restricted Randomization
#' Procedures for Multi-Arm Trials with Equal or Unequal Treatment Allocation Ratios",
#' Statistics in Medicine, Submitted.
#'
#' @param nsim number of simulations (by default equals to 10000)
#'
#' @return a list of simulations results per scenario
#'
#' @examples
#'    # 1000 simulations:
#'
#'    run_example(nsim = 1000)
#'
#' @export

run_example <- function(nsim = 10000) {

  # number of subjects
  nsbj <- 200

  # fixed allocation ratio
  w <- list(
    c( 1,  1,  1,  1),
    c( 2,  1,  1,  2),
    c( 4,  3,  2,  1),
    c(37, 21, 21, 21)
  )

  # number of treatments
  ntrt <- dim(w)[2]

  # randomization procedure(s)
  procedure <- c("CRD",  rep("PBD", 5), rep("BUD", 4), rep("MWUD", 4), rep("DL", 4), rep("DBCD", 5), rep("MaxEnt", 5))
  parameter <- c(NA, 1, 2, 3, 4, 5, 2, 3, 4, 5, 2, 4, 6, 8, 2, 4, 6, 8, 1, 2, 4, 5, 10, 0.05, 0.1, 0.25, 0.5, 1)

  # responses
  # flat responses (N(0, 1)), given treatment id's
  flat_response <- map(seq_len(nsim), ~ matrix(rnorm(4*nsbj), ncol = 4, byrow = TRUE))

  # monotone responses (N(mu[i], 1)), given treatment id's (mu = c(0, 0.1, 0.3, 0.7))
  mu <- c(0, 0.1, 0.3, 0.7)
  monotone_response <- map(seq_len(nsim), ~ matrix(rnorm(4*nsbj, rep(mu, nsbj)), ncol = 4, byrow = TRUE))

  # create scenarios for simulations
  scenario <- cross2(
    w,
    list(
      procedure,
      parameter
    ) %>%
      transpose()
  )
  output <- vector("list", length(scenario))

  for (id in seq_along(scenario)) {
    cat(paste0("Scenario ", id, " out of ", length(scenario)," is being simulated\n"))
    allocation_ratio <- scenario[[id]][[1]]
    rnd_procedure <- scenario[[id]][[2]][[1]]
    rnd_parameter <- scenario[[id]][[2]][[2]]
    out <- restricted:::.restricted_multiple_simulations_cpp(allocation_ratio, nsbj, rnd_procedure, rnd_parameter, nsim)

    target <- paste("w =", str_c(map_chr(allocation_ratio, ~ paste(.)), collapse = ":"))
    procedure_ <- if_else(rnd_procedure == "CRD", "CRD", paste0(rnd_procedure, " (", rnd_parameter, ")"))

    allocation <- map(out, ~ {
      simulation_out <- .
      simulation_out$allocation %>%
        as_data_frame() %>%
        set_names(map_chr(seq_along(w), ~ paste0("N[", ., "]"))) %>%
        add_column(subject = simulation_out$subject, .before = 1) %>%
        add_column(target = target, .before = 1) %>%
        add_column(procedure = procedure_, .before = 1)
    }) %>%
      bind_rows()

    output[[id]] <- list(
      # operational characteristics
      op = map(out, ~ {
        data_frame(
          procedure = procedure_,
          target = target,
          subject = .$subject,
          FI = .$FI,
          imbalance = .$imbalance
        )
      }) %>%
        bind_rows() %>%
        group_by(procedure, target, subject) %>%
        summarise(
          AFI = mean(FI),
          MI = max(imbalance),
          MPM = mean(imbalance)
        ) %>%
        mutate(CMPM = map_dbl(subject, ~ mean(MPM[seq_len(.)]))) %>%
        # add an "average standard deviation" of the allocation proportions
        inner_join(
          allocation  %>%
            gather(allocation, value, -procedure, -target, -subject) %>%
            group_by(procedure, target, allocation, subject) %>%
            mutate(value = value/subject) %>%
            summarise(var_value = var(value)) %>%
            group_by(procedure, subject) %>%
            summarise(ASD = sum(var_value)) %>%
            mutate(ASD = sqrt(subject*ASD)),
          by = c("procedure", "subject")
        ) %>%
        ungroup(),

      # allocation probabilities
      probability = map(out, ~ {
        simulation_out <- .
        simulation_out$probability %>%
          as_data_frame() %>%
          set_names(map_chr(seq_along(w), ~ paste0("pi[",. , "]"))) %>%
          add_column(subject = simulation_out$subject, .before = 1) %>%
          add_column(target = target, .before = 1) %>%
          add_column(procedure = procedure_, .before = 1)
      }) %>%
        bind_rows() %>%
        gather(probability, value, -procedure, -target, -subject) %>%
        group_by(procedure, target, subject, probability) %>%
        summarise(value = mean(value)) %>%
        spread(probability, value) %>%
        ungroup(),

      # allocation proportions distributions
      allocation = allocation %>%
        gather(proportion, value, -procedure, -target, -subject) %>%
        mutate(value = value/subject) %>%
        mutate(proportion = gsub("N", "rho", proportion)) %>%
        group_by(procedure, target, subject, proportion) %>%
        do(distribution = .$value) %>%
        spread(proportion, distribution),

      # treatment assignments
      obseravtion = map(out, ~ {
        data_frame(
          procedure = procedure_,
          target = target,
          subject = .$subject,
          treatment = .$treatment
        )
      }) %>%
        bind_rows() %>%
        group_by(procedure, target, subject) %>%
        do(treatment = .$treatment) %>%
        ungroup()
    )

    # generate responses and perform ANOVA test
    treatment <- out %>% map(~ .$treatment)

    # flat1 is an original flat response (M1 model)
    flat1 <- restricted:::.anova_test_per_subject_cpp (treatment, flat_response, 4L, 0.05, FALSE) %>%
      as_data_frame() %>%
      set_names(c("subject", "flat1"))

    # flat2 is an original flat response with time drift added (M2 model)
    flat2 <- restricted:::.anova_test_per_subject_cpp(treatment, flat_response, 4L, 0.05, TRUE) %>%
      as_data_frame() %>%
      set_names(c("subject", "flat2"))

    # monotone1 is an original monotone response (M1 model)
    monotone1 <- restricted:::.anova_test_per_subject_cpp(treatment, monotone_response, 4L, 0.05, FALSE) %>%
      as_data_frame() %>%
      set_names(c("subject", "monotone1"))

    # monotone2 is an original monotone response with time drift added (M2 model)
    monotone2 <- restricted:::.anova_test_per_subject_cpp(treatment, monotone_response, 4L, 0.05, TRUE) %>%
      as_data_frame() %>%
      set_names(c("subject", "monotone2"))

    # type I error/power
    typeIerror <- inner_join(
      inner_join(
        flat1,
        flat2,
        by = c("subject")
      ),
      inner_join(
        monotone1,
        monotone2,
        by = c("subject")
      ),
      by = c("subject")
    ) %>%
      add_column(target = target, .before = 1) %>%
      add_column(procedure = procedure_, .before = 1)

    output[[id]] <- append(output[[id]], list(typeIerror = typeIerror))
  }

  op <- output %>% map(~ .$op) %>%
    bind_rows()

  probability  <- output %>% map(~ .$probability) %>%
    bind_rows()

  allocation  <- output %>% map(~ .$allocation) %>%
    bind_rows()

  typeIerror <- output %>% map(~ .$typeIerror) %>%
    bind_rows()

  if (!dir.exists("./data")) {
    dir.create("./data")
  }
  save(op, probability, allocation, typeIerror,
       file = "./data/simulation_data_to_analyze.Rda")

  cat("Simulations are finished!")
}



#' A simulation example: analysis and visualization
#'
#' Visualize reproduced simulations from the paper Ryeznik Y, Sverdlov O (2018) "A Comparative Study of Restricted Randomization
#' Procedures for Multi-Arm Trials with Equal or Unequal Treatment Allocation Ratios",
#' Statistics in Medicine, Submitted.
#'
#' @param data_file a path to the *.Rda file where results of \code{run_example} are stored
#'
#'
#' @examples
#'    # by default, run_example function creates a folder "data" in a working folder and stores the results in a files
#'    # "./data/simulation_data_to_analyze.Rda":
#'
#'    visualize_example(data_file = "./data/simulation_data_to_analyze.Rda")
#'
#' @export

visualize_example <- function(data_file = "./data/simulation_data_to_analyze.Rda") {

# load data
  load(file = data_file)

# create "plots" folder in a "data" folder
  if (!dir.exists("./data/plots")) {
    if(!dir.exists("./data")) {
      dir.create("./data")
    }
    dir.create("./data/plots")
  }

# colors to use in plots
color <- c("red", "steelblue", "seagreen", "orange", "purple")

# shapes to use in plots
shape <- c(23, 21, 4, 3, 8, 24, 22)

targets <- op$target %>%
  as.vector() %>%
  unique()

designs <- map_chr(op$procedure, ~ str_split_fixed(., " ", 2)[1]) %>%
  as.vector() %>%
  unique()

response_models <- (typeIerror %>%
  gather(response_model, value, -procedure, -target, -subject))$response_model %>%
  as.vector() %>%
  unique()


## overall performance
wI <- wR <- 1
op1 <- op #%>%
  #mutate(CMPM = CMPM/subject)

ovp <- op1 %>%
  select(-MI, -MPM, -ASD) %>%
  inner_join(
    inner_join(
      # AFIs of CRD and MaxEnt (1)
      op1 %>%
        filter(procedure %in% c("CRD", "MaxEnt (1)")) %>%
        select(-MI, -MPM, -CMPM, -ASD) %>%
        spread(procedure, AFI) %>%
        rename(AFI_CRD = CRD, AFI_MaxEnt1 = `MaxEnt (1)`),
      # CMPMs of CRD and MaxEnt (1)
      op1 %>%
        filter(procedure %in% c("CRD", "MaxEnt (1)")) %>%
        select(-AFI, -MI, -MPM, -ASD) %>%
        spread(procedure, CMPM) %>%
        rename(CMPM_CRD = CRD, CMPM_MaxEnt1 = `MaxEnt (1)`),
      by = c("target", "subject")
    ),
    by = c("target", "subject")
  ) %>%
  mutate(
    design = map_chr(.$procedure, ~ str_split_fixed(., " ", 2)[1]),
    design = factor(design, levels = designs, ordered = TRUE),
    target = factor(target, levels = targets, ordered = TRUE)
  ) %>%
  # calculate rank of the design
  mutate(
    UI = CMPM/(CMPM_CRD-CMPM_MaxEnt1)-CMPM_MaxEnt1/(CMPM_CRD-CMPM_MaxEnt1),
    UR = AFI/(AFI_MaxEnt1-AFI_CRD)-AFI_CRD/(AFI_MaxEnt1-AFI_CRD),
    G = sqrt(((wI*UI)^2 +(wR*UR)^2)/(wI^2 + wR^2))
  )  %>%
  select(target, design, procedure, subject, G) %>%
  filter(!is.nan(G))

# OVPs heat mep plots
ovp_plots <- ovp %>%
  group_by(target) %>%
  do(
    ovp_plot = {
      ggplot(data = ., aes(subject, procedure)) +
        geom_tile(aes(fill=G))+
        scale_fill_gradientn(colors = rainbow(7))+
        scale_x_continuous(breaks = c(1, seq(25, 200, by = 25)))+
        xlab("Sample size")+
        ylab("Design")+
        theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 18),
              axis.title.y = element_text(family = "Helvetica", face = "bold", size = 18),
              axis.text.x  = element_text(family = "Helvetica", face = "bold", size = 14),
              axis.text.y  = element_text(family = "Helvetica", face = "bold", size = 14),
              legend.position = "right",
              legend.title = element_blank(),
              legend.text  = element_text(family = "Helvetica", face = "bold", size = 14))
    }
  ) %>%
  mutate(plot_file = paste0("./data/plots/OVP heatmap (", target, ").pdf"))

map2(ovp_plots$plot_file, ovp_plots$ovp_plot, ggsave, width = 16, height = 9, units = "in", dpi = 300)

# excel table of OVPs in the end of a trial
targets %>%
  map(~ {
  target_val <- .
  ovp %>%
    filter(subject == max(subject) & target == target_val) %>%
    select(-design, -target, -subject) %>%
    arrange(G)
}) %>%
  set_names(targets %>% map_chr(~ str_replace_all(., ":", ","))) %>%
  writexl::write_xlsx("./data/plots/ovp4.xlsx")
## =======


## excel table of Type I error/Power
response_models %>%
  map(~ {
    response_model_val <- .
    typeIerror %>%
      gather(response_model, value, -procedure, -target, -subject) %>%
      mutate(value = round(value, 2)) %>%
      filter(subject == max(subject) & response_model == response_model_val) %>%
      select(-subject, -response_model) %>%
      spread(target, value)# %>%
      #arrange(procedure)
  }) %>%
  set_names(response_models) %>%
  writexl::write_xlsx("./data/plots/alpha.xlsx")
## =======


## Operational Characteristics

id <- op1 %>%
  select(procedure) %>%
  mutate(
    design = map_chr(.$procedure, ~ str_split_fixed(., " ", 2)[1])
  ) %>%
  group_by(design) %>%
  do(id = length(unique(.$procedure))) %>%
  unnest() %>%
  select(id) %>%
  .$id %>%
  map(~ seq_len(.)) %>%
  unlist()


# MPM, AFI vs Sample Size
op_plots <- op1 %>%
  select(procedure, target, subject, AFI, MPM) %>%
  mutate(
    design = map_chr(.$procedure, ~ str_split_fixed(., " ", 2)[1]),
    design = factor(design, levels = designs, ordered = TRUE),
    target = factor(target, levels = targets, ordered = TRUE)
  ) %>%
  gather(op, value, -design, -procedure, -target, -subject) %>%
  mutate(op_txt = if_else(op == "MPM", "Momentum of probability mass", "Average forcing index")) %>%
  group_by(target, op) %>%
  do(
    op_plot = ggplot(data = ., aes(x = subject, y = value, group = procedure))+
      geom_point(aes(shape = procedure, fill = procedure, color = procedure), size = 2.2)+
      geom_line(aes(color = procedure), size = 0.1)+
      scale_shape_manual(values = shape[id])+
      scale_color_manual(values = color[id])+
      scale_fill_manual(values = color[id])+
      facet_wrap(~ design, scales = "free_y", ncol = 1)+
      xlab("Sample size")+
      ylab(unique(.$op_txt))+
      theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 18),
            axis.title.y = element_text(family = "Helvetica", face = "bold", size = 18),
            axis.text.x  = element_text(family = "Helvetica", face = "bold", size = 14),
            axis.text.y  = element_text(family = "Helvetica", face = "bold", size = 14),
            strip.text.x = element_text(family = "Helvetica", face = "bold", size = 18),
            strip.text.y = element_text(family = "Helvetica", face = "bold", size = 18),
            legend.position = "right",
            legend.title = element_blank(),
            legend.text  = element_text(family = "Helvetica", face = "bold", size = 14))
  ) %>%
  mutate(plot_file = paste0("./data/plots/op (", op, ", ", target, ").pdf"))

map2(op_plots$plot_file, op_plots$op_plot, ggsave, width = 16, height = 9, units = "in", dpi = 300)


# AFI vs CMPM

op1 %>%
  filter(subject %in% c(25, 50, 100, 200)) %>%
  select(procedure, target, subject, AFI, CMPM) %>%
  mutate(
    design = map_chr(.$procedure, ~ str_split_fixed(., " ", 2)[1]),
    target = factor(target, levels = targets),
    subject = factor(subject,
                     levels = c(25, 50, 100, 200),
                     labels = map_chr(c(25, 50, 100, 200), ~ paste0("sample size = ", .)))
  ) %>%
  ggplot(aes(x = CMPM, y = AFI, group = procedure))+
    geom_point(aes(shape = design, fill = procedure, color = procedure), size = 4.5)+
    scale_shape_manual(values = shape)+
    scale_color_manual(values = color[id])+
    scale_fill_manual(values = color[id])+
    facet_grid(target ~ subject, scales = "free")+
    xlab("Cumulative momentum of probability mass")+
    ylab("Average forcing index")+
    theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 18),
          axis.title.y = element_text(family = "Helvetica", face = "bold", size = 18),
          axis.text.x  = element_text(family = "Helvetica", face = "bold", size = 14),
          axis.text.y  = element_text(family = "Helvetica", face = "bold", size = 14),
          strip.text.x = element_text(family = "Helvetica", face = "bold", size = 18),
          strip.text.y = element_text(family = "Helvetica", face = "bold", size = 18),
          legend.position = "right",
          legend.title = element_blank(),
          legend.text  = element_text(family = "Helvetica", face = "bold", size = 14))

  ggsave("./data/plots/AFI vs CMPM.pdf", width = 16, height = 9, units = "in", dpi = 300)


# ASD vs AFI for selected subjects

  sub_op11 <- op1 %>% mutate(
    design = map_chr(.$procedure, ~ str_split_fixed(., " ", 2)[1]),
    target = factor(target, levels = targets)
  ) %>%
    filter(subject %in% c(50, 100, 200) & design %in% c("BUD", "CRD", "DBCD", "DL", "MWUD"))

  sub_op12 <- sub_op11 %>%
    filter(subject == 50 & target == targets[1])

  ggplot()+
    geom_point(data = sub_op11, aes(x = AFI, y = ASD, shape = factor(subject), color = procedure), size = 3)+
    geom_line(data = sub_op11, aes(x = AFI, y = ASD, linetype = design, color = procedure), size = 1)+
    geom_point(data = sub_op12, aes(x = AFI, y = ASD, shape = factor(subject), color = procedure), size = 3)+
    geom_text(data = sub_op12, aes(x = AFI, y = ASD, label=procedure, color = procedure),
              size = 4, angle = 30 , family = "Helvetica", check_overlap = FALSE, nudge_y = 0.22, nudge_x = 0.001)+
    labs(shape="Sample size")+
    scale_x_continuous(limits = c(0, 0.17))+
    scale_y_continuous(limits = c(0,1))+
    xlab("Average forcing index")+
    ylab("Average standard deviation of allocation proportions")+
    facet_grid(target ~ ., scales = "free_y")+
    theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 18),
          axis.title.y = element_text(family = "Helvetica", face = "bold", size = 18),
          axis.text.x  = element_text(family = "Helvetica", face = "bold", size = 14),
          axis.text.y  = element_text(family = "Helvetica", face = "bold", size = 14),
          strip.text.x = element_text(family = "Helvetica", face = "bold", size = 18),
          strip.text.y = element_text(family = "Helvetica", face = "bold", size = 18),
          legend.position = "right",
          legend.text  = element_text(family = "Helvetica", face = "bold", size = 14))

  ggsave("./data/plots/ASD vs AFI.pdf", width = 16, height = 9, units = "in", dpi = 300)

  # the same as above but separated by targets
  asd_afi_plots <- op1 %>% mutate(
    design = map_chr(.$procedure, ~ str_split_fixed(., " ", 2)[1]),
    target = factor(target, levels = targets)
  ) %>%
    filter(subject %in% c(50, 100, 200) & design %in% c("BUD", "CRD", "DBCD", "DL", "MWUD")) %>%
    group_by(target) %>%
    do(
      asd_afi_plot = ggplot(data = .)+
        geom_point(aes(x = AFI, y = ASD, shape = factor(subject), color = procedure), size = 3)+
        geom_line(aes(x = AFI, y = ASD, linetype = design, color = procedure), size = 1)+
        geom_text(data = filter(., subject == 50), aes(x = AFI, y = ASD, label=procedure, color = procedure),
                  size = 4, angle = 45 , family = "Helvetica", check_overlap = FALSE, nudge_y = 0.025, nudge_x = 0.005)+
        labs(shape="Sample size")+
        scale_x_continuous(limits = c(0, 0.17))+
        scale_y_continuous(limits = c(0,1))+
        xlab("Average forcing index")+
        ylab("Average standard deviation of allocation proportions")+
        facet_grid(. ~ target, scales = "free_y")+
        theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 18),
              axis.title.y = element_text(family = "Helvetica", face = "bold", size = 18),
              axis.text.x  = element_text(family = "Helvetica", face = "bold", size = 14),
              axis.text.y  = element_text(family = "Helvetica", face = "bold", size = 14),
              strip.text.x = element_text(family = "Helvetica", face = "bold", size = 18),
              strip.text.y = element_text(family = "Helvetica", face = "bold", size = 18),
              legend.position = "right",
              legend.text  = element_text(family = "Helvetica", face = "bold", size = 14))
    ) %>%
    mutate(plot_file = paste0("./data/plots/ASD vs AFI (", target, ").pdf"))
  map2(asd_afi_plots$plot_file, asd_afi_plots$asd_afi_plot, ggsave, width = 16, height = 9, units = "in", dpi = 300)

## =======


## Allocation proportion
  allocation_plots <- allocation %>%
    filter(subject == max(subject)) %>%
    select(-subject) %>%
    gather(proportion, value, -procedure, -target) %>%
    unnest() %>%
    group_by(target) %>%
    do(
      allocation_plot = {
        w = unique(.$target) %>%
          str_replace("w = ", "c(") %>%
          str_replace_all(":", ",") %>%
          str_c(")") %>%
          parse(text = .) %>%
          eval()

        ggplot(data = ., aes(x = procedure, y =value))+
          geom_boxplot(aes(fill = procedure))+
          scale_y_continuous(breaks = w/sum(w))+
          xlab("Design")+
          ylab("Allocation proportion")+
          facet_grid(proportion ~ ., labeller = label_parsed)+
          theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 18),
                axis.title.y = element_text(family = "Helvetica", face = "bold", size = 18),
                axis.text.x  = element_text(family = "Helvetica", face = "bold", size = 14, angle = 75,margin = margin(t = 35, b = 5)),
                axis.text.y  = element_text(family = "Helvetica", face = "bold", size = 14),
                strip.text.x = element_text(family = "Helvetica", face = "bold", size = 18),
                strip.text.y = element_text(family = "Helvetica", face = "bold", size = 22, angle = 0),
                legend.position = "none",
                legend.title = element_blank(),
                legend.text  = element_text(family = "Helvetica", face = "bold", size = 14))
    }) %>%
    mutate(plot_file = paste0("./data/plots/Allocation (", target, ").pdf"))

  map2(allocation_plots$plot_file, allocation_plots$allocation_plot, ggsave, width = 16, height = 9, units = "in", dpi = 300)


## =======

## Unconditional allocation probabilities plots: To Demonstrate an ARP property
# All designs
arp_plots <- probability %>%
  mutate(
    design = map_chr(.$procedure, ~ str_split_fixed(., " ", 2)[1]),
    design = factor(design, levels = designs, ordered = TRUE),
    target = factor(target, levels = targets, ordered = TRUE)
  ) %>%
  gather(probability, value, -design, -procedure, -target, -subject) %>%
  group_by(design, target) %>%
  do(
    arp_plot = {
      w = unique(.$target) %>%
        str_replace("w = ", "c(") %>%
        str_replace_all(":", ",") %>%
        str_c(")") %>%
        parse(text = .) %>%
        eval()
      ggplot(data = ., aes(x = subject, y = value, group = probability))+
       geom_point(aes(color = probability), size = 1.25)+
       geom_line(aes(color = probability), size = 0.5)+
       scale_y_continuous(limits = c(0, 1), breaks = w/sum(w))+
       scale_color_manual(
         values = c("darkblue", "darkgreen", "black", "red"),
         labels = c("treatment 1", "treatement 2", "treatement 3", "treatement 4")
       )+
       ggtitle(paste0("Unconditional allocation probability: ", unique(.$target), ", ", unique(.$design)))+
       xlab("number of subjects")+
       ylab("unconditional allocation probability")+
       facet_wrap(~ procedure, nrow = 1)+
       theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 12),
             axis.title.y = element_text(family = "Helvetica", face = "bold", size = 12),
             axis.text.x  = element_text(family = "Helvetica", face = "bold", size = 10),
             axis.text.y  = element_text(family = "Helvetica", face = "bold", size = 10),
             title = element_text(family = "Helvetica", face = "bold", size = 14),
             strip.text.x = element_text(family = "Helvetica", face = "bold", size = 10),
             strip.text.y = element_text(family = "Helvetica", face = "bold", size = 10),
             legend.position = "bottom",
             legend.title = element_blank(),
             legend.text  = element_text(family = "Helvetica", face = "bold", size = 16))
       }) %>%
  mutate(plot_file = paste0("./data/plots/ARP propety (", design, ", ", target, ").pdf"))

map2(arp_plots$plot_file, arp_plots$arp_plot, ggsave, width = 16, height = 9, units = "in", dpi = 300)


# Selected designs: DBCD, MaxEnt, MWUD
arp_plots <- probability %>%
  mutate(
    design = map_chr(.$procedure, ~ str_split_fixed(., " ", 2)[1]),
    design = factor(design, levels = designs, ordered = TRUE),
    target = factor(target, levels = targets, ordered = TRUE)
  ) %>%
  filter(design %in% c("DBCD", "MaxEnt", "MWUD")) %>%
  gather(probability, value, -design, -procedure, -target, -subject) %>%
  group_by(target) %>%
  do(
    arp_plot = {
      w = unique(.$target) %>%
        str_replace("w = ", "c(") %>%
        str_replace_all(":", ",") %>%
        str_c(")") %>%
        parse(text = .) %>%
        eval()
      ggplot(data = ., aes(x = subject, y = value, group = probability))+
        geom_point(aes(color = probability), size = 1.5)+
        geom_line(aes(color = probability), size = 0.25)+
        scale_y_continuous(limits = c(0, 1), breaks = w/sum(w))+
        scale_color_manual(
          values = c("darkblue", "darkgreen", "black", "red"),
          labels = c(expression(pi[1]), expression(pi[2]), expression(pi[3]), expression(pi[4]))
        )+
        xlab("Sample size")+
        ylab("Unconditional allocation probability")+
        facet_wrap(~ procedure, ncol = 5)+
        theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 18),
              axis.title.y = element_text(family = "Helvetica", face = "bold", size = 18),
              axis.text.x  = element_text(family = "Helvetica", face = "bold", size = 14),
              axis.text.y  = element_text(family = "Helvetica", face = "bold", size = 14),
              strip.text.x = element_text(family = "Helvetica", face = "bold", size = 18),
              strip.text.y = element_text(family = "Helvetica", face = "bold", size = 18),
              #legend.position = "bottom",
              legend.title = element_blank(),
              legend.text  = element_text(family = "Helvetica", face = "bold", size = 20))
    }) %>%
  mutate(plot_file = paste0("./data/plots/3 designs ARP propety (", target, ").pdf"))

map2(arp_plots$plot_file, arp_plots$arp_plot, ggsave, width = 16, height = 9, units = "in", dpi = 300)
}

