#' Restricted randomization
#'
#' Simulates a clinical trial one time or multiple times with retricted randomizaion procedure
#'
#' @author Yevgen Ryeznik (\email{yevgen.ryeznik@gmail.com}), Oleksandr Sverdlov
#'
#' @param w vector of fixed allocation ratio ((ir)rational or integers with GCD = 1).
#' @param nsbj number of subjects to randomize.
#' @param procedure name of randomization procedure \code{c("CRD", "PBD", "BUD", "MWUD",
#'    "DBCD", "DL", "MinQD", "MaxEnt")}.
#' @param parameter tuning parameter of randomization procedure (NA for \code{"CRD"}).
#' @param nsim number of simulations (1 by default).
#'
#' @return In case of a single trial simulation (i.e. \code{nsim} equals to 1) returns
#'    a list with the following items items:
#'    \itemize{
#'       \item \code{op} -- a tibble (data frame) of operational characteristics:
#'       forcing index, imbalance.
#'       \item \code{probability} -- a tibble (data frame) of allocation probabilities at each step.
#'       \item \code{allocation} -- a tibble (data frame) of subjects' allocations at each step.
#'       \item \code{observation} -- a tibble (data frame) of treatment assignments.
#'    }
#'
#'    In case of multiple trial simulations (i.e. \code{nsim} larger than 1) returns
#'    a list with the following items items:
#'    \itemize{
#'       \item \code{op} -- a tibble (data frame) of operational characteristics:
#'       average forcing index (AFI), maximum imbalance (MI), momentum of probability
#'       mass (MPM), cumulative momentum of probability mass (CMPM), and "average
#'       standard deviation" of the allocation proportions (ASD).
#'       \item \code{probability} -- a tibble (data frame) of unconditional allocation
#'       probabilities at each step.
#'       \item \code{allocation} -- a tibble (data frame) of distrbutions of the
#'       subjects' allocations at each step.
#'       \item \code{observation} -- a tibble (data frame) of distributions of treatment
#'       assignments.
#'    }
#' @examples
#'    # a single trial simulation (without responses):
#'    #    target allocation ratio: 4:3:2:1 (i.e. 4 treatment groups)
#'    #    number of subjects involved in  trial: 100
#'    #    randomization procedure used: Drop-the-Loser (DL)
#'    #    randomization procedure parameter: 2
#'
#'    restricted(c(4, 3, 2, 1), 200, "DL", 2)
#'
#'    # 200 trial simulations without responses:
#'    #    target allocation ratio: 31:19 (i.e. 2 treatment groups)
#'    #    number of subjects involved in  trial: 50
#'    #    randomization procedure used: Maximum Entropy Constraint Balance
#'    #          Randomization (MaxEnt)
#'    #    randomization procedure parameter: 0.5
#'    #    number of simulations: 200
#'
#'    restricted(c(31, 19), 50, "MaxEnt", 0.5, nsim = 200)
#'
#' @references Ryeznik Y, Sverdlov O (2018) "A Comparative Study of Restricted Randomization Procedures
#'    for Multi-Arm Trials with Equal or Unequal Treatment Allocation Ratios",
#'    Statistics in Medicine, Submitted.
#'
#' @export
#'
restricted <- function(w, nsbj, procedure, parameter = NA, nsim = 1) {
  target <- paste("w =", str_c(map_chr(w, ~ paste(.)), collapse = ":"))

  if (!procedure %in% c("CRD", "PBD", "BUD", "MWUD", "DL", "GDL", "DBCD", "MinQD", "MaxEnt")) {
    stop('the value of input parameter procedure must be one of  c("CRD", "PBD", "BUD", "MWUD", "DL", "GDL", "DBCD", "MinQD", "MaxEnt")')
  }
  # a single trial simulation
  if (nsim == 1) {
    out <- restricted:::.restricted_one_simulation_cpp(w, nsbj, procedure, parameter)

    list(
      # operational characteristics
      op = data_frame(
        procedure = out$procedure,
        target = target,
        subject = out$subject,
        FI = out$FI,
        imbalance = out$imbalance
      ),

      # allocation probabilities
      probability = out$probability %>%
        as_data_frame() %>%
        set_names(map_chr(seq_along(w), ~ paste0("pi[",. , "]"))) %>%
        add_column(subject = out$subject, .before = 1) %>%
        add_column(target = target, .before = 1) %>%
        add_column(procedure = out$procedure, .before = 1),

      # allocation at each step
      allocation = out$allocation %>%
        as_data_frame() %>%
        set_names(map_chr(seq_along(w), ~ paste0("N[", ., "]"))) %>%
        add_column(subject = out$subject, .before = 1) %>%
        add_column(target = target, .before = 1) %>%
        add_column(procedure = out$procedure, .before = 1),

      # observations: -- treatment assignments at each step
      observation = data_frame(
        procedure = out$procedure,
        target = target,
        subject = out$subject,
        treatment = out$treatment
      )
    )
  }

  # multiple trial simulations
  else {
    out <- restricted:::.restricted_multiple_simulations_cpp(w, nsbj, procedure, parameter, nsim)

    allocation <- map(out, ~ {
      simulation_out <- .
      simulation_out$allocation %>%
        as_data_frame() %>%
        set_names(map_chr(seq_along(w), ~ paste0("N[", ., "]"))) %>%
        add_column(subject = simulation_out$subject, .before = 1) %>%
        add_column(target = target, .before = 1) %>%
        add_column(procedure = simulation_out$procedure, .before = 1)
    }) %>%
      bind_rows()

    list(
      # operational characteristics
      op = map(out, ~ {
        data_frame(
          procedure = .$procedure,
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
          add_column(procedure = simulation_out$procedure, .before = 1)
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
         procedure = .$procedure,
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
  }

}
