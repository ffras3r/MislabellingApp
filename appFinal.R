library(shiny)
library(DT)
library(shinydashboard)
library(ggplot2)
library(tidyverse)
library(rstan)
library(ggmosaic)
library(plotly)
library(car)
library(caret)
library(vcd)
library(shinyalert)


setwd("C:\\Users\\ffras\\OneDrive\\Desktop\\School\\self\\App")

is.binomial <- function(data){
  names <- character(0)
  for(name in colnames(data)){
    uv <- unique(data[[name]])
    if(length(uv) == 2 && all(uv %in% c(0,1))){
      names <- c(names, name)
    }
  }
  if(length(names) > 0){
    return(names)
  } else{
    return("No repsonse selected")
  }
}

standardize_vars <- function(var){
  if(!is.numeric(var) || is.integer(var)){ 
    var <- trimws(var)
    var <- as.factor(var)
  }
  return(var)
}

predict_obs <- function(response, predictors, data, cent){
  train <- sample_frac(data, cent)
  test <- anti_join(data, train)
  
  the_formula <- as.formula(paste(response, "~", paste(predictors, collapse="+")))
  final <- glm(the_formula, data = train, family=binomial)
  predictions <- predict(final, newdata = test, type = 'response')
  
  predictions <- ifelse(predictions > 0.5, 1, 0)
  true_pos <- sum(predictions == 1 & data[[response]] == 1)
  false_pos <- sum(predictions == 1 & data[[response]] == 0)
  true_neg <- sum(predictions == 0 & data[[response]] == 0)
  false_neg <- sum(predictions == 0 & data[[response]] == 1)
  
  results <- numeric(0)
  
  #accuracy
  results[1] <- (true_pos + true_neg) / (true_neg + true_pos + false_neg + false_pos)
  
  #sensitivity
  results[2] <- true_pos / (true_pos + false_neg)
  
  #specificity
  results[3] <- true_neg / (true_neg + false_pos)
  
  return(results)
}

forward_variable_selection <- function(response, predictors, theData, significance_level) {
  selected_vars <- character(0)  # Initialize an empty vector for selected variables
  while(length(predictors) > 0){
    data <- theData
    if(length(selected_vars)>0){
      prev_formula <- as.formula(paste(response, "~", paste(selected_vars, collapse = "+")))
    } else{
      prev_formula <- as.formula(paste(response, "~ 1"))
    }
    prev_model <- glm(prev_formula, data=data, family="binomial")
    best_AIC <- Inf
    best_var <- NULL
    best_model <- NULL
    
    for(var in predictors){
      formula <- as.formula(paste(response, "~", paste(c(selected_vars, var), collapse = "+")))
      
      current_model <- glm(formula, data = data, family = "binomial")
      
      current_AIC <- AIC(current_model)
      if (current_AIC < best_AIC) {
        best_AIC <- current_AIC
        best_var <- var
        best_model <- current_model
      }
    }
    
    p_value <- anova(best_model, prev_model, test='Chisq')$"Pr(>Chi)"[2]
    predictors <- setdiff(predictors, best_var)
    
    if (!is.na(p_value) && p_value < significance_level) {
      selected_vars <- c(selected_vars, best_var)
    } else {
      break
    }
  }
  return(selected_vars)
}

univariate_models <- function(response, predictors, data){
  univariates <- data.frame(matrix(nrow=length(predictors), ncol= 4))
  colnames(univariates) <- c("var", "pvalue", "AIC", "length")
  for(i in 1:length(predictors)){
    formula <- as.formula(paste(response, "~", paste(predictors[i])))
    uni <- glm(formula, family=binomial(link="logit"), data=data)
    
    univariates$var[i] <- predictors[i]
    deviance <- uni$null.deviance - uni$deviance
    dfs <- uni$df.null - uni$df.residual
    
    univariates$pvalue[i] <- signif(1-pchisq(deviance, df=dfs), 5)
    univariates$AIC[i] <- signif(AIC(uni), 0)

    univariates$length[i] <- length(data[[predictors[i]]])
  }
  return(univariates)
}

gen_visuals <- function(responseName, varName, data){
  Data <- data.frame(Status=data[[responseName]], Predictor=data[[varName]])
  #attach(Data)
  Data$Predictor <- standardize_vars(Data$Predictor)
  
  if (is.factor(Data$Predictor)){
    Data$Status = as.factor(Data$Status)
    Data$Status <- ifelse(Data$Status=="1", "Mislabelled", "Labelled")
    
    Table <- table(Data$Predictor,Data$Status)
    
    prop_successes <- Table[, 2] / (Table[, 1] + Table[, 2])
    prop_samples <- (Table[, 2] + Table[, 1]) / sum(Table)
    prop_SE <- sqrt((prop_successes * (1-prop_successes)) / (Table[, 1] + Table[, 2]))
    mosaic_data <- data.frame(Labeled = Table[, 1], Mislabelled = Table[, 2])
    custom_colors <- c("Labelled" = "#0072B2", "Mislabelled" = "darkred")

    p <- ggplot(Data) +
      geom_mosaic(aes(x=product(Status, Predictor), fill = Status)) +
      #scale_fill_gradient(low = "#0072B2", high = "darkred") +
      scale_fill_manual(values = custom_colors)
    labs(title = "Mosaic Plot of Mislabelling Rate Between Groups",
         x = deparse(varName),
         y = "Mislabelling Status") +
      theme_minimal() 
      

    ggplotly(p + theme(axis.text.x = element_text(angle = 60, hjust=1)))
    
    
  } else if (is.numeric(Data$Predictor)){
    Data <- na.omit(Data)
    Data$Status = as.numeric(Data$Status)
    color_palette <- c("#0072B2", "darkred")
    custom_theme <- theme_minimal() +
      theme(
        plot.title = element_text(size = 25, hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        legend.position = "none",  # Remove legend
        panel.grid.major = element_line(color = "lightgray", linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.margin = unit(c(1, 1, 1, 1), "cm")
      )
    
    Data$Status = Data$Status - 1
    
    p <- ggplot(Data, aes(x = Predictor, y = Status)) +
      geom_point(color = color_palette[1], size = 3) +
      stat_smooth(method = "glm", color = color_palette[2], se = TRUE, method.args = list(family = binomial)) +
      labs(
        y = "Mislabeling Rate",
        x = deparse(varName)
      ) +
      custom_theme
    ggplotly(p)
    
    
  }
}

univariate_analysis <- function(responseName, varName, significance, data){
  response <- data[[responseName]]
  var <- data[[varName]]
  print("VAR CLASS")
  print(class(var))
  
  if (is.factor(var)){
    Table <- table(var, response)
    prop_successes <- Table[, 2] / (Table[, 1] + Table[, 2])
    prop_samples <- (Table[, 2] + Table[, 1]) / sum(Table)
    prop_SE <- sqrt((prop_successes * (1-prop_successes)) / (Table[, 1] + Table[, 2]))
    mosaic <- data.frame(Labeled = Table[, 1], Mislabelled = Table[, 2]) #for mosaic plots
    odds <- prop_successes / (1-prop_successes)
    results <- data.frame(Labeled = Table[, 1],
                          Mislabelled = Table[, 2], 
                          Odds = signif(odds, 3),
                          logOdds = signif(log(odds), 3),
                          MislabelledProp = signif(prop_successes, 3), 
                          SE = signif(prop_SE, 3),
                          CIlb = signif((prop_successes-(qnorm(significance / 2, lower.tail=FALSE)*prop_SE)), 3), 
                          CIub = signif((prop_successes+(qnorm(significance / 2, lower.tail=FALSE)*prop_SE)), 3), 
                          SampleProp = signif(prop_samples,3))
    
    return(results)
  }
}

stan_models <- function(Data, response, predictors, prior){
  if(prior == "Informed"){
    stan <- "data {
      int<lower=0> N;
      int<lower=0> P;
      matrix[N, P] X;
      vector[P] betaMu;
      vector[P] betaSig;
  
      int<lower=0, upper=1> isMislabelled[N];
    }
    parameters{
      vector[P] beta;
    }
    model{
      for(i in 1:P){
        beta[i] ~ normal(betaMu[i], betaSig[i]);
      }
      isMislabelled ~ bernoulli_logit(X * beta);
    }"
  } else if (prior == "Uninformed"){
    stan <- "data {
      int<lower=0> N;
      int<lower=0> P;
      matrix[N, P] X;
      vector[P] betaMu;
      vector[P] betaSig;
      int<lower=0, upper=1> isMislabelled[N];
    }
    parameters{
      vector[P] beta;
    }
    model{
      beta ~ normal(0, 1);
      isMislabelled ~ bernoulli_logit(X * beta);
    }"
  } else {
    stan <- "data {
      int<lower=0> N;
      int<lower=0> P;
      matrix[N, P] X;
      vector[P] betaMu;
      vector[P] betaSig;
      int<lower=0, upper=1> isMislabelled[N];
    }
    parameters{
      vector[P] beta;
    }
    model{
      beta ~ uniform(-5, 5);
      isMislabelled ~ bernoulli_logit(X * beta);
    }"
  }
  
  Xform <- paste("~", paste(predictors, collapse = "+"))
  
  theFormula <- glm(as.formula(paste(response, Xform)), data=Data, family=binomial(link="logit"))
  
  X <- model.matrix(as.formula(Xform), data=Data)
  betaMu <- coef(theFormula)
  betaSig <- summary(theFormula)$coef[, "Std. Error"] * sqrt(nrow(Data))
  
  model <- stanc(model_code = stan, model_name = 'model')
  Model <- stan_model(stanc_ret = model)
  
  Fit <- sampling(
    Model,
    data = list(
      N=nrow(X),
      P=ncol(X),
      X=X,
      isMislabelled = Data$isMislabelled,
      betaMu = betaMu,
      betaSig = betaSig),
    chains=4,
    iter=2000,
    warmup=1000,
    cores=4
  )
  return(Fit)
}

ui <- fluidPage(
  dashboardBody(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
    ),
    tabsetPanel(
      tabPanel(
        "Raw Data",
        fluidRow(
          div(style = "display: flex; flex-direction: row;",
            column(4,
                   radioButtons("response", "Select response variable", choices = "No repsonse selected")
                 ),
            column(4,
                   textOutput("ResponseMessage")
                   ),
            column(4,
                   radioButtons("ALLNA", "Remove all missing values of the dataset?", c("No", "Yes")),
                   radioButtons("set_seed", "Would you like to set a seed?", c("No", "Yes")),
                   conditionalPanel(
                     condition = "input.set_seed == 'Yes'",
                     numericInput("seed", label="Enter any Integer", value=0)
                   )
                  ),
                                                                                    
          )
        ),
        fluidRow(
          column(4,
                 "Enter a Dataset",
                 fileInput("file", "Choose a CSV file"),
                 tags$hr()
          )
        ),
        fluidRow(
          column(12, DTOutput("head"))
        )
      ),
      tabPanel(
        "Model Building",
        fluidRow(
          column(4,
                 checkboxGroupInput("predictors", "Select predictors", choices = NULL),
                 #radioButtons("response", "Select response variable", choices = "No repsonse selected"),
                 textOutput("columns"),
                 textOutput("variables")
          ),
          column(4, 
                 fluidRow(
                    h3("Options"),
                    sliderInput("Sig", "Significance level", value=0.05, min=0.001, max=1),
                    actionButton("askSig", "More Info"),
                        
                     radioButtons("removeNA", "Remove NA values from selected predictors?", c("No", "Yes")),
                     
                    actionButton("askNA", "More Info"),
                    radioButtons("Reduced", "Remove Categories with less than 'n' observations?", c("No", "Yes")),
                    actionButton("askSmall", "More Info"),
                    conditionalPanel(
                       condition = "input.Reduced == 'Yes'",
                       numericInput("N", "n: ", value=5, min=1, max=20)
                   )
                 ),
                 fluidRow(
                   h3("Univariate Summaries"),
                   DTOutput("Summary"),
                   actionButton("askUnivariate", "More Info")
                 ),
                 fluidRow(
                   h3("Model Formula"),
                   textOutput("ModelFormula"),
                   actionButton("askFormula", "More Info")
                 ),
                 fluidRow(
                   h3("Multivariate Summary"),
                   dataTableOutput("Model"),
                   actionButton("askModel", "More Info")
                 ),
                 fluidRow(
                   h3("Model AIC"),
                   textOutput("ModelAIC"),
                   actionButton("askAIC", "More Info")
                 )
          ),
          column(4,
                 fluidRow(
                   h3("Warnings"),
                   htmlOutput("Warnings")
                 ))
        )
      ),
      tabPanel(
        "Visualization",
        column(4,
               fluidRow(
                 actionButton("click2", "Generate Visualization")
               ),
               radioButtons("plotPredictor", "Select predictors", choices = "No repsonse selected"),
        ),
        column(4,
               fluidRow(
                 h3("Plot"),
                 plotlyOutput("Plot", width="800px", height="600px")
               ),
               fluidRow(
                 h3("Summary Statistics"),
                 dataTableOutput("SumStats")
               )
        )
      ),
      tabPanel(
        "Bayesian Inference",
        column(4, 
               fluidRow(
                 actionButton("click", "Generate Bayesian Fit")
               ),
               fluidRow(
                 radioButtons("Priors", "Please Select a Prior", choices = c("Informed", "Uninformed", "Uniform"))
               ),
               fluidRow(
                 radioButtons("AnyOrForward", "Any Variables or Forward Selected ones?", choices = c("Any", "Selected"))
               )
        ),
        column(4, 
               fluidRow(
                 textOutput("Message")
               ),
               fluidRow(
                 dataTableOutput("BayesianFit")
               ),
               fluidRow(
                 plotOutput("TracePlot", width="800px", height="600px")
               )
        )
      ),
      tabPanel(
        "Predictions",
        column(4, 
               checkboxGroupInput("pred_predictors", "Select predictors", choices = NULL)
               ),
        column(4,
                radioButtons("fit_forward", "Forward Selected variables or any variable?", choices=c("No", "Yes")),
                sliderInput("cent", "How much data to train on? (%)", value=0.7, min=0.001, max=1),
                fluidRow(
                  h3("Accuracy"),
                  textOutput("accuracy")
                ),
                fluidRow(
                  h3("Sensitivity"),
                  textOutput("sensitivity")
                ),
                fluidRow(
                  h3("Specificity"),
                  textOutput("specificity")
                )
              )
            )
          )
        )
      )
  


server <- function(input, output, session){
  observeEvent(input$askSig, {
    shinyalert("Significance Level", 
              "The P-value is a measurement of significance of a model, defined as the probability that the found relationship is found by chance. 
               In our model-building process, we use a cutoff to only include variables with a low enough p-value. 
               The standard, and reccomended p-value cutoff is 0.05 (5%).
               ", type = "info")
  })
  
  observeEvent(input$askNA, {
    shinyalert("Remove Missing Values", 
              "There cannot be missing values in the variables chosen, a model cannot be made with missing values.
               It would be ideal to remove all missing values, but can reduce size of dataset. 
              It is good practice to only remove missing values in variables under consideration.
              ", type="info")
  })
  
  observeEvent(input$askSmall, {
    shinyalert("Remove Small Categories",
               "Including small categories can create problems, like separation, where entire categories have one Mislabelling status.
               This is a violation of Linear Regression Assumptions and thus analysis is incorrect. 
               Another potential issue is multicollinearity: When an entire category has one predictor value, which is also a Linear Regression Assumption.
               Small categories have higher chances for these issues to occur.
               A reccomended starting value when dealing with these issues is around 3% of total observations.
               ", type="info")
  })
  
  observeEvent(input$askUnivariate, {
    shinyalert("Univariate Summary",
               "The univariate summary reports each variable predicting mislabelling individually.
               It provides the p-value, a measure of significance. 
               The P-value is a measurement of significance of a model, defined as the probability that the found relationship is found by chance. 
               The standard, and reccomended p-value cutoff is 0.05 (5%).
               
               It also reports on the AIC (Akaike Information Criterion), which consideres goodness of fit AND complexity of model (number of parameters).
               A simple model is desired. The scale of the AIC is not important, only the difference between models. A minimal AIC is desired
               
               The length is the number of observations used in the models.
               ", type="info")
  })
  
  observeEvent(input$askFormula, {
    shinyalert("Model Formula", 
               "The Model Formula shows the response variable (in this case, mislabelled status), and the variables used to predict it. In this case, predictor variables are selected using forward variable selection.
               Forward variable selection goes like this. One variable with the lowest AIC is chosen. This single-variable model is tested for significance. if its p-value is less than the cutoff, it is chosen.
               A second variable is chosen with the lowest AIC when predicting mislabelling state together. If this two-predictor model predicts mislabelling rate significantly better than the first, the second variable is chosen.
               This process of adding variables and checking if its addition is justified continues until no more variables are deemed significant.
               ", type="info")
  })
  
  #observeEvent(input$ask)
  
  #observeEvent(input$)
  
  
  
  observe({
    predictorsList <- colnames(data())
    selectedPredictors <- isolate(input$predictors)
    selectedPred_predictors <- isolate(input$pred_predictors)
    #responseList <- setdiff(predictorsList, selectedPredictors)
    responseList <- is.binomial(data())
    updateCheckboxGroupInput(session, "predictors", choices = predictorsList, selected = selectedPredictors)
    updateCheckboxGroupInput(session, "pred_predictors", choices = predictorsList, selected = selectedPred_predictors)
    
    updateRadioButtons(session, "response", choices = responseList)
    
    updateRadioButtons(session, "plotPredictor", choices = predictorsList)
    
    if(input$set_seed == "Yes") set.seed(input$seed)
  })
  
  results <- reactive({
    if(length(input$pred_predictors) > 0){
      predict_obs(input$response, input$pred_predictors, data(), input$cent)
    }
  })
  
  Warnings <- reactive({
   # print("IN THE WARNINGS LOOP")
    warnings <- character(0)
    
    Ma <- "No Variables have been selected at the given significance level. Either select more variables or raise the significance level\n\n"
    Mb <- "Fitted probabilities of 0 or 1 occured. This means a certain category in the given variables give a 0% or 100% mislabelling rate, and it unwanted in probit regression.\n\n"
    Mc <- "There is major separation. This means categories in one variable perfectly match with categories in a second\n\n"
    Md <- "the VIF (or Variance Inflation Factor) is too high. This is due to collinearity, or multiple variables are aligned with each other. This is not wanted in linear regression\n\n"
    Me <- "No Response is Selected\n\n"
    Mf <- "Probabilities in the model are very high. Although a categorical variable might have a low P-value does not mean each category will. You might want to check the significance of each category in the multivariate model summary\n\n"
    Mg <- "There are missing values in one or more of your variables. Select the omit missing values option or choose different variables\n\n"
    
    if(length(vars()) == 0) warnings <- c(warnings, Ma)
    if(is.null(input$response)) warnings <- c(warnings, Me)
    
    if(length(input$predictors) > 0){
      df <- data()[, c(input$response, input$predictors)]
      if(nrow(df) != nrow(na.omit(df))){
        warnings <- c(warnings, Mg) 
      }
    }
    
    tryCatch({
      glm <- forward_variable_selection(input$response, vars(), data(), input$Sig)
      vif(glm)
    }, warning = function(warning_condition){
      #print("AYOOOO")
      #Mb <- "Fitted probabilities of 0 or 1 occured. This means a certain category in the given variables give a 0% or 100% mislabelling rate, and it unwanted in probit regression."
      #Mc <- "There is major separation. This means categories in one variable perfectly match with categories in a second"
      
      if(grepl(warning_condition$message, "glm.fit: fitted probabilities numerically 0 or 1 occurred")){
        #print("ABC")
        warnings <<- c(warnings, Mb)
      } 
      if(grepl(warning_condition$message, "glm.fit: algorithm did not converge")){
        #print("THERE")
        #print(Mc)
        warnings <<- c(warnings, Mc)
        #print(warnings)
      }
    }, error=function(error_condition){
      #Mc <- "There is major separation. This means categories in one variable perfectly match with categories in a second"
      
      if(grepl(error_condition$message, "there are aliased coefficients in the model")){
        print("HERE")
        warnings <<- c(warnings, Mc)
      }
    })
    #print(warnings)
    return(warnings)
  })
  
  dataFull <- reactive({
    req(input$file)
    df <- read.csv(input$file$datapath)
    ext <- tools::file_ext(input$file$datapath)
    switch(ext,
           csv=vroom::vroom(input$file$datapath, delim=","),
           validate("please upload a .csv file")
    )
    for(var in colnames(df)){
      print(var)
      print(class(df[[var]]))
      df[[var]] <- standardize_vars(df[[var]])
      print(class(df[[var]]))
    }
    
    return(df)
  })
  
  data <- reactive({
    # print(c(input$predictors, input$response))
    # a <- c("isMislabelled", "theYear")
    # print(a)
    if(input$removeNA == "Yes"){
      df <- dataFull()
      sub <- df[, c(input$response, input$predictors)]
      df <- df[complete.cases(sub), ]
    } else{
      df <- dataFull()
    }
    
    if(input$ALLNA == "Yes"){ 
      dff <- na.omit(dataFull())
    } else if(input$Reduced == "Yes"){
      for(pred in input$predictors){
        print(pred)
        print(class(df[[pred]]))
        if(is.character(df[[pred]])){
          df <- df %>%
            group_by(across(all_of(pred))) %>%
            filter(n() >= input$N) %>%
            ungroup()
        }
      }
      dff <- df
    } else {
      dff <- df
    }
    return(dff)
  })

  vars <- reactive({
    if (!is.null(input$response) && length(input$predictors) > 0) {
      df <- data()[, c(input$response, input$predictors)]
      if(nrow(df) == nrow(na.omit(df))){
        forward_variable_selection(input$response, input$predictors, data(), input$Sig)
      } else{
        NULL
      }
    } else {
      NULL
    }
  })
  
  pred_vars <- reactive({
    if(input$fit_forward == "Yes"){
      forward_variable_selection(input$response, input$pred_predictors, data(), input$Sig)
    } else{
      input$pred_predictors
    }
  })
  
  fitForward <- reactive({
    if (length(vars()) > 0){
      stan_models(data(), input$response, vars(), input$Priors)
    }
  })
  
  fitAny <- reactive({
    if (!is.null(input$response) && length(input$predictors) > 0) {
      stan_models(data(), input$response, input$predictors, input$Priors)
    } else {
      character(0)
    }
  })
  
  output$ResponseMessage <- renderText(
    "Please select the labelling status variable. Recall this variable must be binomial (eg. contains only 0's and 1's)"
  )
  
  output$head <- renderDT(data(), options = list(pageLength = 50))
  output$columns <- renderText({
    paste(input$predictors, collapse = ", ")
  })
  output$variables <- renderText({
    paste(vars(), collapse = ", ")
  })

  output$Warnings <- renderUI({
    HTML(paste(Warnings(), collapse = "<br/><br/>"))
  })
  
  output$Summary <- renderDT({
    if (!is.null(input$response) && length(input$predictors) > 0) {
      univariate_models(input$response, input$predictors, data())
    }
  })
  
  output$ModelFormula <- renderText({
    if(length(vars()) > 0){
      paste("Formula: ", input$response, "~", paste(vars(), collapse="+"))
    }
  })
  
  output$Model <- renderDataTable({
    
    if (length(vars()) > 0) {
      formula <- as.formula(paste(input$response, "~", paste(vars(), collapse="+")))
      model <- glm(formula, data=data(), family=binomial(link="logit"))
      
      # Extracting the coefficient summary and converting it to a data frame
      coef_summary <- as.data.frame(summary(model)$coefficients)
      
      # Adding the variable names as a new column
      coef_summary$Variable <- rownames(coef_summary)
      coef_summary[, 1:4] <- signif(coef_summary[, 1:4], 3)
      rownames(coef_summary) <- NULL
      
      # Reordering the columns for better display
      coef_summary <- coef_summary[, c("Variable", "Estimate", "Std. Error", "z value", "Pr(>|z|)")]
      
      # Creating a DataTable
      datatable(coef_summary)
    }
  })
  
  output$ModelAIC <- renderText({
    if(!is.null(input$response) && length(vars()) > 0){
      formula <- as.formula(paste(input$response, "~", paste(vars(), collapse="+")))
      model <- glm(formula, data=data(), family=binomial(link="logit"))


      paste("AIC: ", signif(AIC(model),3))
    }
    
    
  })
  
  output$myTable <- renderDT({
    datatable(data(), options = list(pageLength = 10)) %>%
      formatStyle(columns = c("var"), width = '100px') %>%
      formatStyle(columns = c("pvalue", "AIC"), textAlign = "right") %>%
      styleInterval(data()$pvalue, 
                    c('background-color: #FFCCCC', 'background-color: #99FF99'), 
                    values = c(0.05, 0.01)
      )
  })
  
  observeEvent(input$click2, {
    
    output$SumStats <- renderDataTable({
      if (!is.null(input$response) && length(input$response) > 0) {
        univariate_analysis(input$response, input$plotPredictor, 0.05, data())
      }
    })
    
    output$Plot <- renderPlotly({
      gen_visuals(input$response, input$plotPredictor, data())
    })
  })
  
  output$Message <- renderText({
    "Please be patient with the Bayesian Fit, proper convergence can take minutes to complete!"
  })
  
  observeEvent(input$click, {
    output$BayesianFit <- renderDT({
      if(input$AnyOrForward == "Selected"){
        summary(fitForward())$summary
      }else{
        summary(fitAny())$summary
      }
    })
    
    output$TracePlot <- renderPlot({
      if(input$AnyOrForward == "Selected"){
        traceplot(fitForward())
      } else{
        traceplot(fitAny())
      }
    })
  })
  
  output$accuracy <- renderText({
    results()[1]
  })
  output$sensitivity <- renderText({
    results()[2]
  })
  output$specificity <- renderText({
    results()[3]
  })
}

shinyApp(ui, server)
