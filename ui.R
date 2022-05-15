# Define UI for data upload app ----

shinyUI(navbarPage(theme = shinytheme("sandstone"),"Methylation Data",
                   tabPanel("Upload sample",
                            # Data upload Page
                            # App title ----
                            titlePanel("Upload your sample"),
                            
                            # Sidebar layout with input and output definitions ----
                            sidebarLayout( 
                                
                                # Sidebar panel for inputs ----
                                sidebarPanel(
                                    
                                    # Horizontal line ----
                                    tags$hr(),
                                    # Option for test data
                                    
                                   # radioButtons("test_up", "Choose Data", choices = c("Upload","Test")),
                                    #conditionalPanel(
                                    #    condition = "input.test_up == 'Upload'",
                                        
                                    # Input: Upload single sample
                                    
                                    fileInput("sample_up", "Upload Sample File",
                                              multiple = TRUE,
                                              accept = c(".csv", ".CGmap")),
                                   
                                    # Horizontal line ----
                                    tags$hr()
                                    
                                    
                                    
                                ),
                                # Main panel for displaying outputs ----
                                mainPanel(
                                    
                                    # Output: Data file ----
                                    verbatimTextOutput("sample_disp"),
                                    verbatimTextOutput("sample_summary")
                                )
                            )
                   ),
                 tabPanel("Age",
                          # Data upload Page
                          # App title ----
                              mainPanel(
                                  
                                  # Output: Data file ----
                                  plotOutput("age_sample")                              
                                  )
                          
                 ),
                 tabPanel("Sex",
                          mainPanel(
                              plotOutput("sex")                              
                          )
                          
                 
                ),
                tabPanel("Weight",
                        mainPanel(
                            plotOutput("wght")                              
                        )
                ),
                tabPanel("Sterilization Status",
                         mainPanel(
                             plotOutput("strst")
                         )
                ),
                tabPanel("Phylogenetic Tree",
                         titlePanel("Uploading genotype sample"),
                         
                         # Sidebar layout with input and output definitions ----
                         sidebarLayout( 
                             
                             # Sidebar panel for inputs ----
                             sidebarPanel(
                                 
                                 # Horizontal line ----
                                 tags$hr(),
                                 # Option for test data
                                 
                                 # radioButtons("test_up", "Choose Data", choices = c("Upload","Test")),
                                 #conditionalPanel(
                                 #    condition = "input.test_up == 'Upload'",
                                 
                                 # Input: Upload single sample
                                 
                                 fileInput("geno_up", "Upload Sample File",
                                           multiple = TRUE,
                                           accept = c(".xlsx")),
                                 numericInput("levels", "Enter number of levels",value = 4),
                                 numericInput("tree_width", "Enter plot width",             
                                              value = 1.4,
                                              min = 1,
                                              step = 0.1),
                                 
                                 # Horizontal line ----
                                 tags$hr(),
                                 
                             ),
                            mainPanel(
    
                             plotOutput("phylot", height = 1200)
                            )
                         
                         )
                    
                ),
                tabPanel("Behaviour Traits",
                         titlePanel("Uploading Behaviour sample"),
                         
                         mainPanel(
                           
                           plotOutput("behav")
                         )
                             
                            
                    )
            )
    )
