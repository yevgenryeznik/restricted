#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(restricted)

shinyServer(function(input, output) {

  simulation_output <- reactiveValues(
    op = NULL,
    ovp = NULL,
    probability = NULL,
    allocation = NULL,
    typeIerror = NULL
  )

  observeEvent(input$simulate, {
    # selected randomization procedures
    selected <- c(
      input$proc11_check,
      input$proc12_check,
      input$proc2_check,
      input$proc3_check,
      input$proc4_check,
      input$proc5_check,
      input$proc6_check,
      input$proc7_check,
      input$proc8_check
    )

    # names of randomization procedures
    procedure <- c(
      "CRD",
      "MaxEnt",
      input$proc2_name,
      input$proc3_name,
      input$proc4_name,
      input$proc5_name,
      input$proc6_name,
      input$proc7_name,
      input$proc8_name
    )

    # corresponding parameters
    parameter <- as.numeric(
      c("NA",
        1,
        input$proc2_param,
        input$proc3_param,
        input$proc4_param,
        input$proc5_param,
        input$proc6_param,
        input$proc7_param,
        input$proc8_param
      )
    )

    rnd_procedure <- procedure[selected]
    rnd_parameter <- parameter[selected]

    # target allocation ratio
    w <- as.numeric(unlist(strsplit(input$w,",")))
    simulation_output$w <- w

    # number of treatments
    ntrt <- length(w)

    # number of subjects involved
    nsbj <- as.numeric(input$nsbj)

    # number of simulation runs
    nsim <- as.numeric(input$nsim)

    # significance level to test a hypothesis
    alpha <- as.numeric(input$alpha)

    # repsonse is modelled as normal
    means <- as.numeric(unlist(strsplit(input$resp_mean,",")))
    sds <- as.numeric(unlist(strsplit(input$resp_sd,",")))
    add_time_drift <- as.logical(input$add_time_drift)

    # generate responses
    response <- map(
      seq_len(nsim),
      ~ matrix(rnorm(4*nsbj, rep(means, nsbj), rep(sds, nsbj)), ncol = ntrt, byrow = TRUE)
    )

    simulation_output$shape_values <- c(3, 17, 8, 15, 16, 25, 18, 10, 12)[selected]
    simulation_output$color_values <- c(
      "red", "darkorange", "gold3", "darkgreen",
      "blue", "darkblue", "darkviolet", "black", "deepskyblue3"
    )[selected]



    system.time({
      trials <- mcMap(function(id) {
        out <- restricted:::.restricted_multiple_simulations_cpp(w, nsbj, rnd_procedure[id], rnd_parameter[id], nsim)

        procedure_ <- if_else(rnd_procedure[id] == "CRD", "CRD", paste0(rnd_procedure[id], " (", rnd_parameter[id], ")"))
        # ===== Type I error / power =====
        treatment <- out %>% map(~ .$treatment)

        if(!add_time_drift) {
          # an original response (M1 model)
          typeIerror <- restricted:::.anova_test_per_subject_cpp(treatment, response, ntrt, alpha, FALSE) %>%
            as_data_frame() %>%
            set_names(c("subject", "value"))
        }
        else {
          # an original response with time drift added (M2 model)
          typeIerror <- restricted:::.anova_test_per_subject_cpp(treatment, response, ntrt, alpha, TRUE) %>%
            as_data_frame() %>%
            set_names(c("subject", "value"))
        }

        # type I error/power
        typeIerror <- typeIerror %>%
          filter(!is.na(value)) %>%
          add_column(procedure = procedure_, .before = 1)
        # ==========

        # simulation summary
        list(
          # operational characteristics
          op = map(out, ~ {
            data_frame(
              procedure = procedure_,
              subject = .$subject,
              FI = .$FI,
              imbalance = .$imbalance
            )
          }) %>%
            bind_rows() %>%
            group_by(procedure, subject) %>%
            summarise(
              AFI = mean(FI),
              MI = max(imbalance),
              MPM = mean(imbalance)
            ) %>%
            mutate(CMPM = map_dbl(subject, ~ mean(MPM[seq_len(.)]))) %>%
            ungroup(),

          # allocation proportions distributions
          allocation = map(out, ~ {
            simulation_out <- .
            simulation_out$allocation %>%
              as_data_frame() %>%
              set_names(map_chr(seq_along(w), ~ paste0("N[", ., "]"))) %>%
              add_column(subject = simulation_out$subject, .before = 1) %>%
              add_column(procedure = procedure_, .before = 1)
          }) %>%
            bind_rows() %>%
            gather(proportion, value, -procedure, -subject) %>%
            mutate(value = value/subject) %>%
            mutate(proportion = gsub("N", "rho", proportion)) %>%
            group_by(procedure, subject, proportion) %>%
            do(distribution = .$value) %>%
            spread(proportion, distribution),

          # allocation probabilities
          probability = map(out, ~ {
            simulation_out <- .
            simulation_out$probability %>%
              as_data_frame() %>%
              set_names(map_chr(seq_along(w), ~ paste0("pi[",. , "]"))) %>%
              add_column(subject = simulation_out$subject, .before = 1) %>%
              add_column(procedure = procedure_, .before = 1)
          }) %>%
            bind_rows() %>%
            gather(probability, value, -procedure, -subject) %>%
            group_by(procedure, subject, probability) %>%
            summarise(value = mean(value)) %>%
            spread(probability, value) %>%
            ungroup(),

          # type I error / power
          typeIerror = typeIerror
        )
        # ==========
      }, seq_along(rnd_procedure), mc.cores = detectCores())
    })


   simulation_output$op <- trials %>%
      map(~ {.$op}) %>%
      bind_rows()

   simulation_output$probability <- trials %>%
      map(~ {.$probability}) %>%
      bind_rows()


   simulation_output$allocation <- trials %>%
      map(~ {.$allocation}) %>%
      bind_rows()

   simulation_output$typeIerror <- trials %>%
     map(~ {.$typeIerror}) %>%
     bind_rows()

  designs <- unique(simulation_output$op$procedure)

   wI <- wR <- 1
   simulation_output$ovp <- simulation_output$op %>%
     select(-MI, -MPM) %>%
     inner_join(
       inner_join(
         # AFIs of CRD and MaxEnt (1)
         simulation_output$op %>%
           filter(procedure %in% c("CRD", "MaxEnt (1)")) %>%
           select(-MI, -MPM, -CMPM) %>%
           spread(procedure, AFI) %>%
           rename(AFI_CRD = CRD, AFI_MaxEnt1 = `MaxEnt (1)`),
         # CMPMs of CRD and MaxEnt (1)
         simulation_output$op %>%
           filter(procedure %in% c("CRD", "MaxEnt (1)")) %>%
           select(-AFI, -MI, -MPM) %>%
           spread(procedure, CMPM) %>%
           rename(CMPM_CRD = CRD, CMPM_MaxEnt1 = `MaxEnt (1)`),
         by = c("subject")
       ),
       by = c("subject")
     ) %>%
     mutate(
       design = map_chr(.$procedure, ~ str_split_fixed(., " ", 2)[1]),
       design = factor(design, levels = designs, ordered = TRUE)
     ) %>%
     # calculate rank of the design
     mutate(
       UI = CMPM/(CMPM_CRD-CMPM_MaxEnt1)-CMPM_MaxEnt1/(CMPM_CRD-CMPM_MaxEnt1),
       UR = AFI/(AFI_MaxEnt1-AFI_CRD)-AFI_CRD/(AFI_MaxEnt1-AFI_CRD),
       G = sqrt(((wI*UI)^2 +(wR*UR)^2)/(wI^2 + wR^2))
     )  %>%
     select(design, procedure, subject, G) %>%
     filter(!is.nan(G))

  })

  output$allocation_boxplot <- renderPlot({
    if (is.null(simulation_output$allocation)) return()
    simulation_output$allocation %>%
      unnest() %>%
      gather(variable, value, -procedure, -subject) %>%
      filter(subject == max(subject)) %>%
      ggplot(aes(x = factor(procedure), y = value, fill = factor(procedure)))+
        geom_boxplot()+
      scale_y_continuous(breaks = simulation_output$w/sum(simulation_output$w))+
        xlab("randomization procedure")+
        ylab("allocation proportion")+
      facet_grid(variable ~ ., labeller = label_parsed)+
      coord_flip()+
      theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 18),
            axis.title.y = element_text(family = "Helvetica", face = "bold", size = 18),
            axis.text.x  = element_text(family = "Helvetica", face = "bold", size = 14),
            axis.text.y  = element_text(family = "Helvetica", face = "bold", size = 14),
            strip.text.x = element_text(family = "Helvetica", face = "bold", size = 18),
            strip.text.y = element_text(family = "Helvetica", face = "bold", size = 22, angle = 0),
            legend.title = element_blank(),
            legend.text  = element_text(family = "Helvetica", face = "bold", size = 14))
  })

  output$probability_plot <- renderPlot({
    if (is.null(simulation_output$probability)) return()
    simulation_output$probability %>%
      gather(probability, value, -procedure, -subject) %>%
    ggplot(aes(x = subject, y = value, group = probability))+
      geom_line(size = 0.5)+
      geom_point(size = 1.25)+
      scale_y_continuous(limits = c(0, 1), breaks = simulation_output$w/sum(simulation_output$w))+
      xlab("sample size")+
      ylab("uncond. allocation probability")+
      facet_grid(procedure ~ probability, labeller = label_parsed)+
      theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 18),
            axis.title.y = element_text(family = "Helvetica", face = "bold", size = 18),
            axis.text.x = element_text(family = "Helvetica", face = "bold", size = 14),
            axis.text.y = element_text(family = "Helvetica", face = "bold", size = 14),
            strip.text.x = element_text(family = "Helvetica", face = "bold", size = 22),
            strip.text.y = element_text(family = "Helvetica", face = "bold", size = 16))
  })

  output$op_plot <- renderPlot({
    if (is.null(simulation_output$op)) return()
    simulation_output$op %>%
      select(-MI, -MPM) %>%
      gather(variable, value, -procedure, -subject) %>%
      mutate(variable = if_else(variable == "AFI", "Average Forcing Index", "Cumulative Momentum of Probability Mass")) %>%
    ggplot()+
      geom_line(aes(x = subject, y = value, color = procedure), size = 0.5)+
      geom_point(aes(x = subject, y = value, shape = procedure, color = procedure), size = 1.25)+
      scale_shape_manual(values = simulation_output$shape_values)+
      scale_color_manual(values = simulation_output$color_values)+
      xlab("sample size")+
      ylab("")+
      facet_wrap(~ variable, scales = "free_y", ncol = 1)+
      theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 18),
            axis.title.y = element_text(family = "Helvetica", face = "bold", size = 18),
            axis.text.x = element_text(family = "Helvetica", face = "bold", size = 14),
            axis.text.y = element_text(family = "Helvetica", face = "bold", size = 14),
            legend.title = element_blank(),
            legend.text = element_text(family = "Helvetica", face = "bold", size = 18),
            strip.text = element_text(family = "Helvetica", face = "bold", size = 18))
  })

  output$typeIerror_plot <- renderPlot({
    if (is.null(simulation_output$typeIerror)) return()
    simulation_output$typeIerror %>%
    ggplot(aes(x = subject, y = value))+
      geom_line(size = 0.5)+
      geom_point(size = 1.25)+
      xlab("sample size")+
      ylab("")+
      facet_wrap(~ procedure, scales = "free_y", ncol = 1)+
      theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 18),
            axis.title.y = element_text(family = "Helvetica", face = "bold", size = 18),
            axis.text.x = element_text(family = "Helvetica", face = "bold", size = 14),
            axis.text.y = element_text(family = "Helvetica", face = "bold", size = 14),
            legend.title = element_blank(),
            legend.text = element_text(family = "Helvetica", face = "bold", size = 18),
            strip.text = element_text(family = "Helvetica", face = "bold", size = 18))
  })

  output$ovp_table <- renderTable({
    if (is.null(simulation_output$ovp)) return()
    simulation_output$ovp %>%
      filter(subject == max(subject)) %>%
      select(procedure, G) %>%
      arrange(G)
    }, spacing = "m", digits = 3)

  output$ovp_plot <- renderPlot({
    if (is.null(simulation_output$ovp)) return()
    simulation_output$ovp %>%
      ggplot(aes(x = subject, y = factor(procedure))) +
      geom_tile(aes(fill=G))+
      coord_fixed(ratio=25)+
      scale_fill_gradientn(colors = rainbow(7))+
      #scale_x_continuous(breaks = c(1, seq(25, 200, by = 25)))+
      xlab("Sample size")+
      ylab("Design")+
      theme(axis.title.x = element_text(family = "Helvetica", face = "bold", size = 18),
            axis.title.y = element_text(family = "Helvetica", face = "bold", size = 18),
            axis.text.x  = element_text(family = "Helvetica", face = "bold", size = 14),
            axis.text.y  = element_text(family = "Helvetica", face = "bold", size = 14),
            legend.position = "right",
            legend.title = element_blank(),
            legend.text  = element_text(family = "Helvetica", face = "bold", size = 14))

  })

})

