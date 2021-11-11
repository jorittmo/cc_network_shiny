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



def_net <- function(weight_mat, ntask, sig_var, sig_def, show_cor= FALSE){
  require(igraph)
  
  diag(weight_mat) <- 0
  
  if (show_cor == FALSE) {
    line_type <- ifelse(sig_var < 0.05, 1, 2)
    gr_mode <- "undirected"
    plt_sub <- "Solid lines show significant dissociations and the weights are the dissociation sizes. \n Significant deficits are denoted by red nodes"
    wt <- as.matrix(weight_mat)
    wt <- c(wt[which(lower.tri(weight_mat))])
    
    wi <- abs(wt)
    wi <- wi + min(wi, na.rm = TRUE) + 1 #offset 
  } else {
    
    gr_mode <- "plus"
    plt_sub <- "Solid red lines show significant dissociations and the weights are the dissociation sizes, \n dotted grey lines are the correlations. Significant deficits are denoted by red nodes"
    wt <- as.matrix(weight_mat)
    wt <- c(wt)
    
    wi <- abs(wt)
    wi[which(upper.tri(weight_mat))] <- NA
    wi <- wi[which(wi != 0 | is.na(wi))]
    wi <- wi + min(wi, na.rm = TRUE) + 1 #offset 
    
    wt <- wt[which(wt != 0)]
    line_type <- rep(NA, length(wi))
    line_type[which(!is.na(wi))] <- ifelse(sig_var < 0.05, 1, 2)
    line_type[is.na(line_type)] <- 3
  }
  
  ad_mat <- matrix(1, nrow = ntask, ncol = ntask)
  diag(ad_mat) <-  0
  dimnames(ad_mat) <- dimnames(weight_mat)
  
  g <- graph.adjacency(ad_mat, mode=gr_mode, diag = FALSE)
  E(g)$weight <- wt
  E(g)$width <- wi      
  
  E(g)$color[!is.na(wi)] <- "#ff4a2e"
  E(g)$color[is.na(wi)] <- 'grey'
  
  E(g)$lty <- line_type
  
  V(g)$color[sig_def < 0.05] <-  "#F47174"
  V(g)$color[sig_def > 0.05] <-  "#ffb861"
  
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






ui <- fluidPage(
  
  
  fluidRow(column(1,
                  numericInput("n", label = "Number of controls", value = 25),
                  numericInput("r", label = "Correlation", value = 0.25, min = 0, max = 0.99),
                  radioButtons("typecor", "Type of correlation",
                               choices = list("Equal" = "equal", "Simplex" = "simplex"), selected = "equal"),
                  actionButton("simulate", "Simulate!")),
           column(3,
                  plotOutput("plot2", width = "300px", height = "300px")),
           column(8,
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
    if (input$typecor == "equal") {
      Sigma <- (1 - diag(input$nvar)) * input$r
      diag(Sigma) <- 1
    } 
    if (input$typecor == "simplex") {
      Sigma <- toeplitz(c(1, sapply(1:(input$nvar-1), function(x) input$r^x)))
    }
    get_dat(MASS::mvrnorm(input$n, rep(0, input$nvar),
                          Sigma = Sigma), defs())
  } )
  
  
  
  output$tab <- renderTable({
    round(data()$p_mat, 3)
  }, na= "", caption = "Lower triangular shows dissociation p-values, upper triangular shows correlations")
  
  output$plot <- renderPlot({
    def_net(data()$weight_mat, data()$ntask, data()$sig_var, data()$sig_def, show_cor = input$showcor)
  }, res = 96)
  
  output$plot2 <- renderPlot({
    per_sig(data()$sig_var)
  }, res = 96)
  
  
}

shinyApp(ui, server)





