library(shiny)
library(tidyverse)

get_dat <- function(controls, deficits = rep(0, ncol(controls)) ) {
  if (ncol(controls) != length(deficits)) stop("nvar and defs don't match")
  require(singcar)
  
  case <- controls[1, ] + deficits
  controls <- controls[-1, ]
  
  ntask <- ncol(controls)
  ntest <- factorial(ntask)/(factorial(ntask-2)*factorial(2))
  
  
  weight_mat <- cor(controls) # Adjacency matrix consisting of correlation and effect sizes
  p_mat <- cor(controls)
  for (j in 1:(ntask-1)) { 
    for (i in (j+1):ntask){
      
      r <- cor(controls[ , i], controls[ , j])
      std_a <- (case[j] - mean(controls[ , j])) / sd(controls[ , j])
      std_b <- (case[i] - mean(controls[ , i])) / sd(controls[ , i])
      
      weight_mat[i, j] <- (std_a - std_b) / sqrt(2 - 2*r)
      p_mat[i, j] <- RSDT(case[j], case[i], controls[ , j], controls[ , i])$p.value
    }
  }
  diag(p_mat) <- NA
  sig_var <- c(p_mat[lower.tri(p_mat)])
  sig_def <- apply(rbind(case, controls), 2, function(x){
    TD(case = x[1], controls = x[-1], alternative = "two.sided", conf_int = FALSE)$p.value
  } )
  
  list(weight_mat = weight_mat, p_mat = p_mat, ntask = ntask, sig_var = sig_var, sig_def = sig_def)
}

get_sim_dat <- function(n, r, nvar, nsim, defs, typecor = "equal", bonferroni = TRUE, alpha = 0.05){
  
  ntest <- factorial(nvar)/(factorial(nvar-2)*factorial(2))
  alpha <- alpha
  alpha <- ifelse(bonferroni, alpha/ntest, alpha)
  
  if (typecor == "equal") {
    Sigma <- (1 - diag(nvar)) * r
    diag(Sigma) <- 1
  } 
  if (typecor == "simplex") {
    Sigma <- toeplitz(c(1, sapply(1:(nvar-1), function(x) r^x)))
  }
  
  sig_df <- tibble(.rows = nsim)
  for (j in 1:(nvar-1)) { 
    for (i in (j+1):nvar){
      sig_df <- sig_df %>%  mutate("{j}_vs_{i}" := 0)
    }
  }
  
  weight_mat <- array(dim = c(nvar, nvar, nsim))
  for (i in 1:nsim) {
    x <- get_dat(MASS::mvrnorm(n, rep(0, nvar),
                               Sigma = Sigma), defs)
    pmat <- x$p_mat
    sig_df[i, ] <- as.list(pmat[lower.tri(pmat)])
    weight_mat[, , i] <- x$weight_mat
    
  }
  if (nsim != 1){
    sig_var <- NULL
    sig_def <- NULL
  } else {
    sig_var <- x$sig_var
    sig_def <- x$sig_def
  }
  
  weight_mat_means <- rowMeans(weight_mat, dims = 2)
  
  pwr <- sig_df %>% pivot_longer(cols = everything(),
                                 names_to = "dissoc", values_to = "p_value") %>% 
    group_by(dissoc) %>% 
    summarise(
      power = mean(p_value < alpha)
    )
  
  list(pwr = pwr, weight_mat_means = weight_mat_means,
       sig_var = sig_var, sig_def = sig_def,
       p_mat = pmat, ntask = nvar)
}


def_net <- function(weight_mat, ntask, sig_var = NULL, sig_def = NULL, show_cor= FALSE){
  require(igraph)
  
  diag(weight_mat) <- 0
  
  if (show_cor == FALSE) {
    if (!is.null(sig_var)) line_type <- ifelse(sig_var < 0.05, 1, 2)
    gr_mode <- "undirected"
    if(!is.null(sig_var)){
      plt_sub <- "Solid lines show significant dissociations and the weights are the dissociation sizes. \n Significant deficits are denoted by red nodes"
    } else {
      plt_sub <- "Average strength of dissociation shown by edge width and weights as Z[dcc]"
    }
    
    wt <- as.matrix(weight_mat)
    wt <- c(wt[which(lower.tri(weight_mat))])
    
    wi <- abs(wt)
    wi <- wi + min(wi, na.rm = TRUE) + 1 #offset 
  } else {
    
    gr_mode <- "plus"
    if(!is.null(sig_var)){
      plt_sub <- "Solid red lines show significant dissociations and the weights are the dissociation sizes, \n dotted grey lines are the correlations. Significant deficits are denoted by red nodes"
    } else {
      plt_sub <- "Average strength of dissociation shown by edge width and weights as Z[dcc]"
    }
    wt <- as.matrix(weight_mat)
    wt <- c(wt)
    
    wi <- abs(wt)
    wi[which(upper.tri(weight_mat))] <- NA
    wi <- wi[which(wi != 0 | is.na(wi))]
    wi <- wi + min(wi, na.rm = TRUE) + 1 #offset 
    
    wt <- wt[which(wt != 0)]
    
    if (!is.null(sig_var)) {
      line_type <- rep(NA, length(wi))
      line_type[which(!is.na(wi))] <- ifelse(sig_var < 0.05, 1, 2)
      line_type[is.na(line_type)] <- 3
    }
  }
  
  ad_mat <- matrix(1, nrow = ntask, ncol = ntask)
  diag(ad_mat) <-  0
  dimnames(ad_mat) <- dimnames(weight_mat)
  
  g <- graph.adjacency(ad_mat, mode=gr_mode, diag = FALSE)
  E(g)$weight <- wt
  E(g)$width <- wi      
  
  E(g)$color[!is.na(wi)] <- "#ff4a2e"
  E(g)$color[is.na(wi)] <- 'grey'
  
  if (!is.null(sig_var)) E(g)$lty <- line_type
  
  if (!is.null(sig_def)){
    V(g)$color[sig_def < 0.05] <-  "#F47174"
    V(g)$color[sig_def > 0.05] <-  "#ffb861"
  }
  
  
  plot(g, layout = layout.circle,
       edge.label = round(E(g)$weight, 2))
  title("Network of dissociations",
        sub = plt_sub)
}



per_sig <- function(nsig){
  ggplot()+
    geom_bar(aes(x = as.factor(nsig<0.05), fill = as.factor(nsig<0.05)))+
    labs(x = "Significant dissociations") +
    theme_bw() +
    ylim(0, length(nsig))+
    theme(legend.position = "")
}

pwrplot <- function(pwr){
  ggplot(pwr, aes(x = dissoc, y = power)) +
    geom_col(position = "dodge2")+
    ggtitle("Frequency of found dissociations per test")+
    ylim(c(0, 1))+
    theme_bw()
}




ui <- fluidPage(
  
  
  fluidRow(column(1,
                  numericInput("n", label = "Number of controls", value = 25),
                  numericInput("nsim", label = "Number of iterations", value = 1),
                  numericInput("r", label = "Correlation", value = 0.25, min = 0, max = 0.99),
                  radioButtons("typecor", "Type of correlation",
                               choices = list("Equal" = "equal", "Simplex" = "simplex"), selected = "equal"),
                  numericInput("alpha", label = "Alpha", value = 0.05, min = 0, max = 1, step = 0.01),
                  checkboxInput("bonferroni", "Bonferroni correction", value = TRUE),
                  actionButton("simulate", "Simulate!")),
           column(6,
                  plotOutput("plot2")),
           column(5,
                  tableOutput("tab"))
  ),
  fluidRow(
    sidebarLayout(
      sidebarPanel(
        numericInput("nvar", "Number of tasks", value = 5),
        # place to hold dynamic inputs
        uiOutput("inputGroup"), width = 1
      ),
      mainPanel(checkboxInput("showcor", "Show correlation", value = FALSE),
                plotOutput("plot", width = "700px", height = "700px"),  width = 11)
    )
  )
  
)





server <- function(input, output, session) {
  
  #observe changes in "numInputs", and create corresponding number of inputs
  observeEvent(input$nvar, {
    output$inputGroup = renderUI({
      input_list <- lapply(1:input$nvar, function(i) {
        # for each dynamically generated input, give a different name
        inputName <- paste0("Deficit on task ", i)
        numericInput(inputName, inputName, 0)
      })
      do.call(tagList, input_list)
    })
  })
  
  defs <- eventReactive(input$simulate, {
    unlist(
      lapply(1:input$nvar, function(i) {
        inputName <- paste0("Deficit on task ", i)
        input[[inputName]]
      })
    )
  })
  
  data <- eventReactive(input$simulate, {
    get_sim_dat(input$n, input$r, input$nvar, input$nsim, defs = defs(),
                typecor = input$typecor, bonferroni = input$bonferroni, alpha = input$alpha)
  } )
  
  nosims <- eventReactive(input$simulate, {
    nosims <- input$nsim
  }) 
  
  output$tab <- renderTable({
    if (nosims() == 1){
      round(data()$p_mat, 3) 
    } else {
      NULL
    }
  }, na= "", caption = "Lower triangular shows dissociation p-values, upper triangular shows correlations")
  
  output$plot <- renderPlot({
    def_net(data()$weight_mat, data()$ntask, data()$sig_var, data()$sig_def, show_cor = input$showcor)
  }, res = 96)
  
  output$plot2 <- renderPlot({
    if (nosims() == 1){
      per_sig(data()$sig_var)  
    } else {
     pwrplot(data()$pwr) 
    }
  }, res = 96) 
  
  
}

shinyApp(ui, server)





