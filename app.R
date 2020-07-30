
# Load packages
library(shiny) # Build shiny app
library(shinyWidgets)
library(shinydashboard) # For shiny dashboard layout
library(shinydashboardPlus) # For shiny dashboard layout
library(DT) # to generate interactive datatables
library(dplyr) # for data wrangling
library(plotly) # to generate interactive plots
library('survival') # for Cox PH model
library('Hmisc') # for Cox PH model
library(shinycssloaders) # for loading spinner while app is loading
library(randomForest) # to generate RSF classifier models
library(gprofiler2) # for Functional Enrichment analysis of drive genes
library(shinyalert) # for error alert
library(shinyjs) # for error message
library(pROC) # for ROC curve

## import the data

clinical<-read.csv('/home/ruth/METABRIC/data_clinical_patient.txt', header=T, skip = 4, sep='\t', na.strings=c(""," ","NA")) # clinical patient data 
clinical <- na.omit(clinical) # remove missing values 



combined_data <- read.csv("/home/ruth/METABRIC/combined_data.csv") # combined clinical and genomic data for RF model (oversampled)
combined_data <- combined_data[,-1] # remove first column

gene_expression <-read.table('/home/ruth/METABRIC/data_mRNA_median_Zscores.txt', sep = '\t', check.names = FALSE, header = TRUE, row.names = NULL) # gene expression data 

# driver genes identified by OncodriveCLUSTL and OncodriveFML
driver=c("MAP2K4","ARID1A", "PIK3CA", "TBX3", "MAP3K1", "CBFB", "TP53", "KMT2C",
         "AKT1", "GATA3", "RUNX1", "PTEN", "CDH1", "NF1", "PIK3R1", "RB1", "CDKN1B",
         "NCOR1", "CDKN2A", "ERBB2", "ERBB3", "FOXO3", "SMAD4", "KRAS", "BRCA2",
         "BAP1", "GPS2", "AGTR2", "ZFP36L1", "MEN1", "CHEK2", "SF3B1", "AHNAK2",
         "SYNE1", "MUC16")

gene_expression=gene_expression[,-2]#remove Entrez_Gene_Id column
gene_expression = gene_expression %>% filter (gene_expression$Hugo_Symbol %in% driver) # filter for driver genes
rownames(gene_expression) = gene_expression$Hugo_Symbol # set Hugo Symbol as row names
gene_expression=gene_expression[,-1] #remove Hugo_Symbol column
gene_expression = as.matrix(gene_expression) # intersect of clinical & genomic

clin_exprs_intersect=intersect(colnames(gene_expression), clinical$PATIENT_ID) #intersect of samples between gene expression and clinical data 
clin_exprs_match<-match(clin_exprs_intersect,colnames(gene_expression)) # match intersect with gene expression column names 
gene_expression = gene_expression[,clin_exprs_match] # filter gene expression to only contain samples corresponding to clinical data 
clin_exprs_match2<-match(clin_exprs_intersect,clinical$PATIENT_ID) # match intersect with clinical data column names 
clinical_exp <- clinical[clin_exprs_match2,] # filter clinical data to only contain samples that have corresponding gene expression data 




shinyApp(
  # Define UI for application
    ui = dashboardPagePlus(md = TRUE,
                           enable_preloader = TRUE,
                           dashboardHeaderPlus(
                              title = h3(id ='head', "EHR", style = "font-family: Goblin One; font-weight: bold") 
                           ),
                           # Create side-bar menu with all tab options: 
                           dashboardSidebar(
                             useShinyalert(),
                               sidebarSearchForm(textId = "searchpatient", buttonId = "searchButton",
                                                 label = "Enter Patient ID..."), # Search Bar for user to enter Patient ID for patient-level data 
                             br(),
                               sidebarMenu(
                                 menuItem(
                                   text = "Home",
                                   tabName = "home",
                                   icon = icon("info-circle")
                                 ),
                                   menuItem(
                                       text = "Clinical Data", 
                                       tabName = "clinicaldata",
                                       icon = icon("notes-medical"),
                                       startExpanded = TRUE,
                                       menuSubItem("Patient-Level", tabName = "Patientlevel"),
                                       menuSubItem("Cohort-Level", tabName = "Cohortlevel")
                                   ),
                                   menuItem(
                                       text = "Genomic Data", 
                                       tabName = "genomicdata",
                                       icon = icon("dna"),
                                       startExpanded = TRUE,
                                       menuSubItem("Patient-Level", tabName = "Patientlevelg"),
                                       menuSubItem("Cohort-Level", tabName = "Cohortlevelg")
                                   ),
                                   menuItem(
                                       text = "10-year survival prediction", 
                                       tabName = "boxelements",
                                       icon = icon("th"),
                                       startExpanded = TRUE,
                                       menuSubItem("Clinical Prediction", tabName = "clinical_ten"),
                                       menuSubItem("Genomic Prediction", tabName = "genomic_ten"),
                                       menuSubItem("Combined Prediction", tabName = "combined_ten")
                                   )
                               )
                           ),
                           dashboardBody(
                               setShadow(class = "dropdown-menu"),
                               setShadow(class = "box"),
                               
                               tabItems(
                                 # Home Tab - provides user with information on the app, contents of each tab etc. 
                                 tabItem(tabName = "home",
                                         boxPlus(
                                           title = h3("Welcome to EHR"), 
                                           closable = FALSE, 
                                           width = 12,
                                           status = "info", 
                                           solidHeader = TRUE, 
                                           collapsible = TRUE,
                                           h4("An integrated interactive shiny dashboard for visualisation and analysis of the Clinical and Genomic METABRIC Data and generation of ten-year breast cancer survival predictive models"),
                                           accordion(
                                             accordionItem(
                                               id = 1,
                                               title = "Clinical Data",
                                               color = "info",
                                               collapsed = TRUE,
                                               h4("Clinical Data can be viewed at a patient- and cohort-level in their corresponding tabs"),
                                               br(),
                                               h4("To view the data for a specific patient, enter their unique Patient ID into the searchbar on the left"),
                                               br(),
                                               h4("The following Clinical Variables are displayed:"),
                                               hr(),
                                               h4("Patient Information:"),
                                               h4(tags$li("Age of the patient at diagnosis time")),
                                               h4(tags$li("Inferred menopausal state at diagnosis time", em("(post/pre)"))),
                                               h4(tags$li("Overall Survival Time since diagnosis", em("(in months)"))),
                                               h4(tags$li("Overall Status", em("(Whether the patient is alive or deceased)"))),
                                               h4(tags$li("Vital Status", em("(Whether the patient is living, died-of-disease or died of other causes)"))),
                                               hr(),
                                               h4("Breast Cancer Classification:"),
                                               h4(tags$li("Integrative Cluster subtype", em("(Molecular subtype of the cancer based on gene expression (Dawson et al., 2013))"))),
                                               h4(tags$li("Claudin (PAM50) subtypes", em("(luminal A, luminal B, HER2-enriched, basal-like, Claudin-low and Normal)"))),
                                               h4(tags$li("Three Gene classifier subtype", em("(ER-/HER2-, ER+/HER2- High Prolif, ER+/HER2- Low Prolif,HER2+)"))),
                                               hr(),
                                               h4("Treatment Information:"),
                                               h4(tags$li("Chemotherapy", em("(Whether or not the patient had chemotherapy as a treatment (yes/no))"))),
                                               h4(tags$li("Breast cancer surgery type", em("(Mastectomy or Breast Conserving)"))),
                                               h4(tags$li("Radio-Therapy", em("(Whether or not the patient had radio-therapy (yes/no))"))),
                                               h4(tags$li("Hormone Therapy", em("(Whether or not the patient had hormone therapy (yes/no))"))),
                                               hr(),
                                               h4("Tumour Information:"),
                                               h4(tags$li("Nottingham prognostic index (NPI)", em("(used to determine prognosis following surgery for breast cancer. Its value is calculated using three pathological criteria: the size of the tumour; the number of involved lymph nodes; and the grade of the tumour)"))),
                                               h4(tags$li("Primary Tumour Laterality", em("(Whether it is involving the right breast or the left breast)"))),
                                               h4(tags$li("Histological Subtype", em("(Type of the cancer based on microscopic examination of the cancer tissue (It takes a value of 'Ductal/NST', 'Mixed', 'Lobular', 'Tubular/ cribriform', 'Mucinous', 'Medullary', 'Other', 'Metaplastic'))"))),
                                               h4(tags$li("Cellularity", em("(Cancer cellularity post chemotherapy, which refers to the amount of tumor cells in the specimen and their arrangement into clusters)"))),
                                               h4(tags$li("ER measured by IHC", em("(if estrogen receptors are expressed on cancer cells measured using immune-histochemistry (positive/negative))"))),
                                               h4(tags$li("HER2 measured by SNP6", em("(if the cancer is positive for HER2 or not by using advance molecular techniques (Type of next generation sequencing))"))),
                                               h4(tags$li("Lymph Nodes Examined Positive for cancer")),
                                               
                                             ),
                                             accordionItem(
                                               id = 2,
                                               title = "Genomic Data",
                                               color = "info",
                                               collapsed = TRUE,
                                               h4("Pre-processing of the mRNA gene expression data:"),
                                               h4(tags$li("The METABRIC gene expression mRNA z-score data is filtered to contain only the 35 driver genes identified by OncodriveCLUST and OncodriveFML.")),
                                               hr(),
                                               h4("35 driver genes are classified as either up-regulated or down-regulated for each individual patient:"),
                                               h4(tags$li("Genes with z-score mRNA expression values > median value for the given gene are classified as up-regulated")),
                                               h4(tags$li("Genes with z-score mRNA expression values < median value for the given gene are classified as down-regulated")),
                                               br(),
                                               h4("Genomic Patient-level tab displays the 35 driver genes and their assigned regulation status corresponding to the Patient ID entered into the Search Bar"),
                                               hr(),
                                               h4("Enrichment Analysis:"),
                                               h4(tags$li("Functional Enrichment Analysis of the 35 driver genes is performed, The driver genes are mapped to known functional information sources and significantly enriched biological terms and pathways can be identified.")),
                                               hr(),
                                               h4("Cox Proportional Hazard Model:"),
                                               h4(tags$li("a Cox proportional hazard model is used to identify which of the 35 driver genes were significantly associated with survival rate in the breast cancer patients and the results are displayed in a table"))
                                               ),
                                             accordionItem(
                                               id = 4,
                                               title = "Survival Prediction Models",
                                               color = "info",
                                               collapsed = TRUE,
                                               h4("In the 10-year survival prediction tabs, Random Survival Forest Classifier Models can be generated using the selected inputs to predict 10-year breast cancer survival. Models can be generated using only clinical predictors, only genomic predictors and combined predictors")
                                             ),
                                             accordionItem(
                                               id = 5,
                                               title = "Access the Code",
                                               color = "info",
                                               collapsed = TRUE,
                                               socialButton(
                                                 url = "https://github.com/RuthColleran/EHR-METABRIC-Dashboard",
                                                 type = "github"
                                               )
                                             )
                                           )
                                         ),
                                         ),
                                   tabItem(
                                     # Clinical Data - patient-level: displays corresponding clinical data to Patient ID entered in the search bar 
                                     tabName = "Patientlevel",
                                           widgetUserBox(
                                               title = "Patient Information",
                                               subtitle = "Disease - Breast Cancer",
                                               type = NULL,
                                               width = 6,
                                               src = "https://image.flaticon.com/icons/svg/149/149070.svg",
                                               color = "aqua-active",
                                               br(),
                                               boxProfileItem("Patient ID:", textOutput("patient_profile_id")),
                                               br(),
                                               boxProfileItem("Age At Diagnosis:", textOutput("patient_age")),
                                               br(),
                                               boxProfileItem("Inferred Menopausal State At Diagnosis:", textOutput("patient_menopause")),
                                               br(),
                                               boxProfileItem("Overall Survival Time since Diagnosis (Months):", textOutput("patient_os_months")),
                                               boxProfileItem("Vital Status:", textOutput("patient_vital_status")),
                                           hr(),
                                           # Patient survival timeline
                                           boxProfileItemList(title = NULL,
                                               timelineBlock(
                                                 reversed = TRUE,
                                                 timelineItem(
                                                   title = "Breast Cancer Diagnosis",
                                                   icon = "comment-medical",
                                                   color = "olive",
                                                   em("Patient Diagnosed with Breast Cancer")
                                                 ),
                                                 timelineItem(
                                                   title = "Overall Survival Time Since Diagnosis (Months)",
                                                   icon = "id-card-alt",
                                                   color = "teal",
                                                   textOutput("patient_os_time")),
                                                 timelineItem(
                                                   title = "Current Vital Status:",
                                                   icon = "star-of-life",
                                                   color = "purple",
                                                   textOutput("patient_vital_timeline")
                                                 ),
                                                 timelineStart(color = "gray")
                                               ))),
                                           boxPlus(title = "Breast Cancer Classification",
                                                   closable = FALSE,
                                                   collapsible = TRUE,
                                                   status = "info",
                                                   width = 6,
                                                   solidHeader = TRUE,
                                             boxProfileItem("Integrative Cluster Subtype", textOutput("patient_intclust")),
                                             boxProfileItem("Claudin Subtype", textOutput("patient_claudin")),
                                             boxProfileItem("Three-Gene Subtype", textOutput("patient_three_gene"))
                                           ),
                                    
                                           boxPlus(title = "Treatment Information",
                                                   closable = FALSE,
                                                   collapsible = TRUE, 
                                                   status = "info", 
                                                   width = 6,
                                                   solidHeader = TRUE,
                                                   boxProfileItem("Chemotherapy", textOutput("patient_chemo")),
                                                   boxProfileItem("Radio-therapy", textOutput("patient_radiotherapy")),
                                                   boxProfileItem("Hormone Therapy", textOutput("patient_hormonetherapy")),
                                                   boxProfileItem("Breast Surgery Type", textOutput("patient_breastsurgery"))),
                                           boxPlus(title = "Tumour Information",
                                                   closable = FALSE, 
                                                   collapsible = TRUE,
                                                   status = "info", 
                                                   solidHeader = TRUE,
                                                   width = 6,
                                                   boxProfileItem("Nottingham prognostic index (NPI)", textOutput("patient_npi")),
                                                   boxProfileItem("Primary Tumour Laterality", textOutput("patient_laterality")),
                                                   boxProfileItem("Histological Subtype", textOutput("patient_histological")),
                                                   boxProfileItem("Cellularity", textOutput("patient_cellularity")),
                                                   boxProfileItem("ER measured by IHC", textOutput("patient_ER_IHC")),
                                                   boxProfileItem("HER2 measured by SNP6", textOutput("patient_HER2_SNP6")),
                                                   boxProfileItem("Lymph Nodes Examined Positive", textOutput("patient_lymph_nodes")))
                                    
                                    ),
                               tabItem(
                                 # Clinical Data - cohort-level: displays visualisations of METABRIC cohort clinical data 
                                 tabName = "Cohortlevel",
                                       tabBox(width = 12,
                                       tabPanel("Cohort Overview",
                                                fluidRow(valueBoxOutput("total_patients"),
                                                         valueBoxOutput("os_alive_box"),
                                                         valueBoxOutput("os_deceased_box")),
                                       fluidRow(boxPlus(title = "Cohort Overall Survival Time",
                                                        closable = FALSE,
                                                        collapsible = TRUE,
                                                        status = "info", 
                                                        solidHeader = TRUE,
                                                        boxProfileItem("Minimum Overall Survival time (Months)", textOutput("os_min_cohort")),
                                                        br(),
                                                        boxProfileItem("Median Overall Survival time (Months)", textOutput("os_median_cohort")),
                                                        br(),
                                                        boxProfileItem("Maximum Overall Survival time (Months)", textOutput("os_max_cohort"))
                                                ),
                                                boxPlus(title = "Cohort Vital Status",
                                                        closable = FALSE,
                                                        collapsible = TRUE,
                                                        status = "info", 
                                                        solidHeader = TRUE,
                                                        boxProfileItem("Living", textOutput("os_vital_living")),
                                                        boxProfileItem("Died of Disease", textOutput("os_vital_died_disease")),
                                                        boxProfileItem("Died of Other Causes", textOutput("os_vital_died_other")))
                                                ),
                                       fluidRow(
                                         boxPlus(
                                           width = 4,
                                           title = "Overall Survival Time (Months)",
                                           collapsible = TRUE,
                                           solidHeader = TRUE,
                                           closable = FALSE,
                                           status = "info",
                                           plotlyOutput("os_time_plot")  %>% withSpinner(color="#0dc5c1")
                                         ),
                                         boxPlus(
                                           width = 4,
                                           title = "Age at Diagnosis",
                                           collapsible = TRUE,
                                           solidHeader = TRUE,
                                           closable = FALSE,
                                           status = "info",
                                           plotlyOutput("age_plot")  %>% withSpinner(color="#0dc5c1")
                                         ),
                                         boxPlus(
                                           width = 4,
                                           title = "Inferred Menopausal State",
                                           collapsible = TRUE,
                                           solidHeader = TRUE,
                                           closable = FALSE,
                                           status = "info",
                                           plotlyOutput("menopause_plot")  %>% withSpinner(color="#0dc5c1")
                                         )),
                                                ), 
                                       tabPanel("Tumour Information",
                                       fluidRow(box(
                                           width = 4,
                                           title = "Lymph Nodes Examined Positive",
                                           collapsible = TRUE,
                                           DT::dataTableOutput("lymph_nodes_table")  %>% withSpinner(color="#0dc5c1")
                                       ),
                                       box(
                                           width = 4,
                                           title = "Cellularity",
                                           collapsible = TRUE,
                                           plotlyOutput("cellularity_plot")  %>% withSpinner(color="#0dc5c1")
                                       ),
                                       box(
                                           width = 4,
                                           title = "NPI",
                                           collapsible = TRUE,
                                           plotlyOutput("NPI_plot")  %>% withSpinner(color="#0dc5c1")
                                       )),
                                       fluidRow(
                                       box(
                                         width = 4,
                                         title = "Tumor's Histological Subtype",
                                         collapsible = TRUE,
                                         DT::dataTableOutput("histological_table")  %>% withSpinner(color="#0dc5c1")
                                         ),
                                       box(
                                         width = 4,
                                         title = "ER status measured by IHC",
                                         collapsible = TRUE,
                                         plotlyOutput("ER_IHC_plot")  %>% withSpinner(color="#0dc5c1")
                                       ),
                                       box(
                                         width = 4,
                                         title = "HER2 status measured by SNP6",
                                         collapsible = TRUE,
                                         plotlyOutput("HER2_SNP6_plot")  %>% withSpinner(color="#0dc5c1")
                                       ))),
                                       tabPanel("Treatment Information",
                                                fluidRow(
                                       box(
                                           width = 6,
                                           title = "Chemotherapy",
                                           collapsible = TRUE,
                                           plotlyOutput("chemo_plot")  %>% withSpinner(color="#0dc5c1")
                                       ),
                                       box(
                                           width = 6,
                                           title = "Hormone Therapy",
                                           collapsible = TRUE,
                                           plotlyOutput("hormone_therapy_plot")  %>% withSpinner(color="#0dc5c1")
                                       )),
                                       fluidRow(
                                       box(
                                           width = 6,
                                           title = "Radio Therapy",
                                           collapsible = TRUE,
                                           plotlyOutput("radio_therapy_plot")  %>% withSpinner(color="#0dc5c1")
                                       ),
                                       box(
                                           width = 6,
                                           title = "Breast Surgery",
                                           collapsible = TRUE,
                                           plotlyOutput("breast_surgery_plot"))  %>% withSpinner(color="#0dc5c1")
                                       )),
                                       
                                       tabPanel("Breast Cancer Classification",
                                                fluidRow(
                                       box(
                                           width = 8,
                                           title = "Claudin Subtype",
                                           collapsible = TRUE,
                                           plotlyOutput("claudin_subtype_plot")  %>% withSpinner(color="#0dc5c1")
                                       ),
                                       box(
                                         width = 4,
                                         title = "3-gene classifier subtype",
                                         collapsible = TRUE,
                                         plotlyOutput("three_gene_plot")  %>% withSpinner(color="#0dc5c1")
                                         )),
                                       fluidRow(
                                       box(
                                           width = 8,
                                           title = "Integrative cluster",
                                           collapsible = TRUE,
                                           plotlyOutput("int_cluster_plot")  %>% withSpinner(color="#0dc5c1")
                                           ),
                                       box(
                                           width = 4,
                                           title = "Primary Tumor Laterality",
                                           collapsible = TRUE,
                                           plotlyOutput("laterality_plot")  %>% withSpinner(color="#0dc5c1")
                                           )))
                                     )),
                               tabItem(
                                 # Genomic Data - patient-level: displays 35 driver genes and regulation status corresponding to Patient ID entered in searchbar
                                 tabName = "Patientlevelg",
                                       boxPlus(title = "Upregulated Driver Genes",
                                               collapsible = TRUE,
                                               closable = FALSE,
                                               status = "success",
                                               solidHeader = TRUE,
                                               DT::dataTableOutput("upregulated_genes_table")
                                                    ),
                                                    boxPlus(
                                                      title = "Downregulated Driver Genes",
                                                      collapsible = TRUE,
                                                      closable = FALSE,
                                                      status = "primary", 
                                                      solidHeader = TRUE,
                                                      DT::dataTableOutput("downregulated_genes_table")
                                                    )),
                               # Genomic Data - Cohort-level: Functional Enrichment Analysis of 35 driver genes and Cox PH model to identfy driver genes significantly associated with survival
                               tabItem(tabName = "Cohortlevelg",
                                       boxPlus(
                                         title= "Enrichment Analysis of 35 driver genes",
                                         width = 12,
                                         collapsible = TRUE,
                                         solidHeader = TRUE,
                                         closable = FALSE,
                                         status = "info",
                                         plotlyOutput("g_profiler_plot")  %>% withSpinner(color="#0dc5c1") # interactive plot of functional enrichment analysis
                                       ),
                                       boxPlus(
                                         title = "Cox Proportional Hazard Model",
                                         closable = FALSE,
                                         collapsible = TRUE,
                                         status = "info", 
                                         solidHeader = TRUE,
                                         width = 6,
                                         DT::dataTableOutput("surv_gene_results_table")  %>% withSpinner(color="#0dc5c1")
                                       ),
                                       boxPlus(title= "Enrichment Analysis of Driver Genes Significantly Associated with Survival",
                                               closable = FALSE,
                                               collapsible = TRUE,
                                               status = "info", 
                                               solidHeader = TRUE,
                                               width = 6,
                                        plotlyOutput("g_profiler_surv_plot")  %>% withSpinner(color="#0dc5c1")
                                       )),
                               tabItem(
                                 # RSF classifier model for 10-year breast cancer survival using only clinical variables
                                 tabName = "clinical_ten",
                                       boxPlus(
                                         title = "RSF Classifier: Clinical Predictors",
                                         width = 5,
                                         collapsible = TRUE,
                                         closable = FALSE,
                                         solidHeader = TRUE,
                                         status = "info",
                                         checkboxGroupInput(
                                           "clinical_predictors_ten", # Input: user specifies clinical predictors used to generate model
                                           h3("Select clinical predictors to include:"),
                                           c("Age at Diagnosis" = "AGE_AT_DIAGNOSIS",
                                             "NPI" = "NPI",
                                             "Lymph Nodes Examined Positive" = "LYMPH_NODES_EXAMINED_POSITIVE",
                                             "Cellularity" = "CELLULARITY",
                                             "Chemotherapy" = "CHEMOTHERAPY",
                                             "ER IHC" = "ER_IHC",
                                             "HER2 SNP6" = "HER2_SNP6",
                                             "Hormone Therapy" = "HORMONE_THERAPY",
                                             "Inferred Menopausal State" = "INFERRED_MENOPAUSAL_STATE",
                                             "Integrative Cluster" = "INTCLUST",
                                             "Claudin Subtype" = "CLAUDIN_SUBTYPE",
                                             "Three Gene Status" = "THREEGENE",
                                             "Laterality" = "LATERALITY",
                                             "Radiotherapy" = "RADIO_THERAPY",
                                             "Histological Subtype" = "HISTOLOGICAL_SUBTYPE",
                                             "Breast Surgery" = "BREAST_SURGERY"
                                           )
                                         ),
                                         hr(),
                                         h4("You have selected:"),
                                         (verbatimTextOutput("selected_clin")), # display selected predictors
                                         useShinyjs(),
                                         hidden(h4("Please Select Predictors above to generate model", # display message if user does not select predictors
                                                   id = "message1",
                                                   style = "bold;color:blue;")),
                                         actionButton("clinical_ten_button", "Generate Model", icon("paper-plane"), 
                                                      style="color: #fff; background-color: #70DBDB; border-color: #2e6da4")), # Action button to generate RSF model
                                       tabBox(width = 7,
                                              tabPanel("Model Information",
                                       fluidRow(
                                         boxPlus(
                                         title = "Model Information",
                                         collapsible = TRUE,
                                         solidHeader = TRUE,
                                         closable = FALSE,
                                         status = "info",
                                         width = 12,
                                         (verbatimTextOutput("clinical_model_summary")) # summary of RSF model
                                       ))),
                                       tabPanel("Variable Importance",
                                         fluidRow(boxPlus(
                                           title = "Variable Importance",
                                           width = 12,
                                           collapsible = TRUE,
                                           solidHeader = TRUE,
                                           closable = FALSE,
                                           status = "info",
                                           DT::dataTableOutput("var_imp_clinical_table") # table of variable importance for RSF model
                                         )))),
                                       boxPlus(
                                         title = "Receiver Operating Characteristic (ROC) curve",
                                         collapsible = TRUE,
                                         solidHeader = TRUE,
                                         closable = FALSE,
                                         width = 7,
                                         status = "info",
                                         plotOutput("roc_clinical_plot"), # plot ROC curve from test set prediction
                                         hr(),
                                         h3(strong("AUC:"), textOutput("auc_clin")) # AUC from test set prediction
                                       ),
                                       
                                       ),
                               tabItem(
                                 # RSF classifier model for 10-year breast cancer survival using only genomic variables
                                 tabName = "genomic_ten",
                                       box(
                                         title = "RSF Classifier: Genomic Predictors",
                                         width = 5,
                                         collapsible = TRUE,
                                         solidHeader = TRUE,
                                         status = "info",
                                         closable = FALSE,
                                         checkboxGroupInput(
                                           "genomic_predictors_ten", # Input: user specifies genomic predictors used to generate model
                                           h3("Select driver genes to include in prediction:"),
                                           c("MAP2K4"="MAP2K4",
                                             "TBX3"="TBX3", 
                                             "MAP3K1"="MAP3K1", 
                                             "KMT2C"="KMT2C", 
                                             "AKT1"="AKT1", 
                                             "PTEN"="PTEN", 
                                             "NF1"="NF1", 
                                             "CDKN1B"="CDKN1B", 
                                             "SMAD4"="SMAD4",
                                             "KRAS"="KRAS", 
                                             "BAP1"="BAP1", 
                                             "MEN1"="MEN1", 
                                             "SYNE1"="SYNE1", 
                                             "MUC16"="MUC16"
                                           )),
                                         hr(),
                                         h4("You have selected:"),
                                         (verbatimTextOutput("selected_genes")), # display selected predictors
                                         useShinyjs(),
                                         hidden(h4("Please Select Predictors above to generate model", # display message if user does not select predictors
                                                   id = "message2",
                                                   style = "bold;color:blue;")),
                                         actionButton("genomic_ten_button", "Generate Model", icon("paper-plane"), 
                                                      style="color: #fff; background-color: #70DBDB; border-color: #2e6da4")), # Action button to generate RSF model
                                       tabBox(width = 7,
                                              tabPanel("Model Information",
                                                       fluidRow(
                                                         boxPlus(
                                                           title = "Model Information",
                                                           collapsible = TRUE,
                                                           solidHeader = TRUE,
                                                           closable = FALSE,
                                                           status = "info",
                                                           width = 12,
                                                           (verbatimTextOutput("genomic_model_summary")) # summary of RSF model
                                                         ))),
                                              tabPanel("Variable Importance",
                                                       fluidRow(boxPlus(
                                                         title = "Variable Importance",
                                                         width = 12,
                                                         collapsible = TRUE,
                                                         solidHeader = TRUE,
                                                         closable = FALSE,
                                                         status = "info",
                                                         DT::dataTableOutput("var_imp_genomic_table") # table of variable importance for RSF model
                                                       )))),
                                       boxPlus(
                                         title = "Receiver Operating Characteristic (ROC) curve", 
                                         collapsible = TRUE,
                                         solidHeader = TRUE,
                                         closable = FALSE,
                                         width = 7,
                                         status = "info",
                                         plotOutput("roc_genomic_plot"), # plot ROC curve from test set prediction
                                         hr(),
                                         h3(strong("AUC:"), textOutput("auc_genomic")) # AUC from test set prediction
                                       ),      
                               ),
                               tabItem(
                                 # RSF classifier model for 10-year breast cancer survival using combined clinical and genomic variables
                                 tabName = "combined_ten",
                                       boxPlus(
                                         title = "RSF Classifier: Combined Predictors",
                                         width = 5,
                                         collapsible = TRUE,
                                         closable = FALSE,
                                         status = "info",
                                         solidHeader = TRUE,
                                         checkboxGroupInput(
                                           "combined_predictors_ten", # Input: user specifies combined predictors used to generate model
                                           h3("Select combined predictors to include:"),
                                           c("Age at Diagnosis" = "AGE_AT_DIAGNOSIS",
                                             "NPI" = "NPI",
                                             "Lymph Nodes Examined Positive" = "LYMPH_NODES_EXAMINED_POSITIVE",
                                             "Cellularity" = "CELLULARITY",
                                             "Chemotherapy" = "CHEMOTHERAPY",
                                             "ER IHC" = "ER_IHC",
                                             "HER2 SNP6" = "HER2_SNP6",
                                             "Hormone Therapy" = "HORMONE_THERAPY",
                                             "Inferred Menopausal State" = "INFERRED_MENOPAUSAL_STATE",
                                             "Integrative Cluster" = "INTCLUST",
                                             "Claudin Subtype" = "CLAUDIN_SUBTYPE",
                                             "Three Gene Status" = "THREEGENE",
                                             "Laterality" = "LATERALITY",
                                             "Radiotherapy" = "RADIO_THERAPY",
                                             "Histological Subtype" = "HISTOLOGICAL_SUBTYPE",
                                             "Breast Surgery" = "BREAST_SURGERY",
                                             "MAP2K4"="MAP2K4",
                                             "TBX3"="TBX3", 
                                             "MAP3K1"="MAP3K1", 
                                             "KMT2C"="KMT2C", 
                                             "AKT1"="AKT1", 
                                             "PTEN"="PTEN", 
                                             "NF1"="NF1", 
                                             "CDKN1B"="CDKN1B", 
                                             "SMAD4"="SMAD4",
                                             "KRAS"="KRAS", 
                                             "BAP1"="BAP1", 
                                             "MEN1"="MEN1", 
                                             "SYNE1"="SYNE1", 
                                             "MUC16"="MUC16"
                                           )
                                         ),
                                         hr(),
                                         h4("You have selected:"),
                                         (verbatimTextOutput("selected_combined")),  # display selected predictors
                                         useShinyjs(),
                                         hidden(h4("Please Select Predictors above to generate model", # display message if user does not select predictors
                                                   id = "message3",
                                                   style = "bold;color:blue;")),
                                         actionButton("combined_ten_button", "Generate Model", icon("paper-plane"), 
                                                      style="color: #fff; background-color: #70DBDB; border-color: #2e6da4")), # Action button to generate RSF model
                                       tabBox(width = 7,
                                              tabPanel("Model Information",
                                                       fluidRow(
                                                         boxPlus(
                                                           title = "Model Information",
                                                           collapsible = TRUE,
                                                           solidHeader = TRUE,
                                                           closable = FALSE,
                                                           width = 12,
                                                           status = "info",
                                                           (verbatimTextOutput("combined_model_summary")) # summary of RSF model
                                                           ))),
                                              tabPanel("Variable Importance",
                                                       fluidRow(
                                                         boxPlus(
                                                           title = "Variable Information", 
                                                           collapsible = TRUE,
                                                           solidHeader = TRUE,
                                                           closable = FALSE,
                                                           width = 12,
                                                           status = "info",
                                                           DT::dataTableOutput("var_imp_combined_table") # table of variable importance for RSF model
                                                           )))),
                                 boxPlus(
                                   title = "Receiver Operating Characteristic (ROC) curve",
                                   collapsible = TRUE,
                                   solidHeader = TRUE,
                                   closable = FALSE,
                                   width = 7,
                                   status = "info",
                                   plotOutput("roc_combined_plot"), # plot ROC curve from test set prediction
                                   hr(),
                                   h3(strong("AUC:"), textOutput("auc_combined")) # AUC from test set prediction
                                   ),
                                 )
                               ),
    )),
    
    server = function(input, output) {
      
      
      observeEvent(input$searchButton, {
        # if empty/invalid patient ID entered into searchbar and action button is pressed display error message
        if (is.null(input$searchpatient)) {
          shinyalert("Oops!", "Please Enter a Valid Patient ID", type = "error") 
        } else if(!(input$searchpatient %in% clinical$PATIENT_ID)){
          shinyalert("Invalid Patient ID Entered", "Please Enter a Valid Patient ID", type = "error")
          }
    
      })            
        ### Clinical Data: Patient-level           
        patient_for_analysis <- reactive({
          # filter clinical data based on patient ID entered in searchbar
            clinical %>% filter(PATIENT_ID == input$searchpatient) 
        })
        
        output$patient_profile_id <- renderText({
          # display Patient ID corresponding to patient ID query in search bar
          patient_data <- patient_for_analysis()
          patient_data$PATIENT_ID <- as.character(patient_data$PATIENT_ID)
          
        })
        
        output$patient_age <- renderText({
          # display patient age at diagnosis corresponding to patient ID query in search bar
            patient_data <- patient_for_analysis()
            patient_data$AGE_AT_DIAGNOSIS
        })
        
        output$patient_menopause <- renderText({
          # display patient inferred menopausal state at diagnosis corresponding to patient ID query in search bar
          patient_data <- patient_for_analysis()
          patient_data$INFERRED_MENOPAUSAL_STATE <- as.character(patient_data$INFERRED_MENOPAUSAL_STATE)
          patient_data$INFERRED_MENOPAUSAL_STATE
        })
        
        output$patient_os_months <- renderText({
          # display patient OS months corresponding to patient ID query in search bar
          patient_data <- patient_for_analysis()
          patient_data$OS_MONTHS
        })
        output$patient_vital_status <- renderText({
          # display patient vital status corresponding to patient ID query in search bar
          patient_data <- patient_for_analysis()
          patient_data$VITAL_STATUS <- as.character(patient_data$VITAL_STATUS)
          patient_data$VITAL_STATUS
        })
        
        output$patient_os_time <- renderText({
          # display patient OS time (months) in timeline (corresponding to patient ID query in search bar)
          patient_data <- patient_for_analysis()
          patient_data$OS_MONTHS
        })
        output$patient_vital_timeline <- renderText({
          # display patient vital status in timeline (corresponding to patient ID query in search bar)
          patient_data <- patient_for_analysis()
          patient_data$VITAL_STATUS <- as.character(patient_data$VITAL_STATUS)
          patient_data$VITAL_STATUS
        })
        
        # Patient-level: Breast Cancer Classifcation
        output$patient_intclust <- renderText({
          # display patient integrative cluster corresponding to patient ID query in search bar
          patient_data <- patient_for_analysis()
          patient_data$INTCLUST <- as.character(patient_data$INTCLUST)
          patient_data$INTCLUST
        })
        
        output$patient_claudin <- renderText({
          # display patient claudin subtype corresponding to patient ID query in search bar
          patient_data <- patient_for_analysis()
          patient_data$CLAUDIN_SUBTYPE <- as.character(patient_data$CLAUDIN_SUBTYPE)
          patient_data$CLAUDIN_SUBTYPE
        })
        output$patient_three_gene <- renderText({
          # display patient three gene status corresponding to patient ID query in search bar
          patient_data <- patient_for_analysis()
          patient_data$THREEGENE <- as.character(patient_data$THREEGENE)
          patient_data$THREEGENE
        })
        
        # Treatment: patient-level
        output$patient_chemo <- renderText({
          # display patient chemotherapy info corresponding to patient ID query in search bar
            patient_data <- patient_for_analysis()
            patient_data$CHEMOTHERAPY <- as.character(patient_data$CHEMOTHERAPY)
            patient_data$CHEMOTHERAPY
        })
        
        output$patient_radiotherapy <- renderText({
          # display patient radiotherapy info corresponding to patient ID query in search bar
            patient_data <- patient_for_analysis()
            patient_data$RADIO_THERAPY <- as.character(patient_data$RADIO_THERAPY)
            patient_data$RADIO_THERAPY
        })
        
        output$patient_breastsurgery <- renderText({
          # display patient breast surgery info corresponding to patient ID query in search bar
            patient_data <- patient_for_analysis()
            patient_data$BREAST_SURGERY <- as.character(patient_data$BREAST_SURGERY)
            patient_data$BREAST_SURGERY
        })
        
        output$patient_hormonetherapy <- renderText({
          # display patient hormone therapy info corresponding to patient ID query in search bar
          patient_data <- patient_for_analysis()
          patient_data$HORMONE_THERAPY <- as.character(patient_data$HORMONE_THERAPY)
          patient_data$HORMONE_THERAPY
        })
        
        # Tumour Information: patient-level
        output$patient_ER_IHC <- renderText({
          # display patient ER status measured by IHC, corresponding to patient ID query in search bar
          patient_data <- patient_for_analysis()
          patient_data$ER_IHC <- as.character(patient_data$ER_IHC)
          patient_data$ER_IHC
        })
        
        output$patient_HER2_SNP6 <- renderText({
          # display patient HER2 status measured by SNP6, corresponding to patient ID query in search bar
          patient_data <- patient_for_analysis()
          patient_data$HER2_SNP6 <- as.character(patient_data$HER2_SNP6)
          patient_data$HER2_SNP6
        })
        
        output$patient_histological <- renderText({
          # display patient tumour histological subtype, corresponding to patient ID query in search bar
            patient_data <- patient_for_analysis()
            patient_data$HISTOLOGICAL_SUBTYPE <- as.character(patient_data$HISTOLOGICAL_SUBTYPE)
            patient_data$HISTOLOGICAL_SUBTYPE
        })
        
        output$patient_laterality <- renderText({
          # display patient primary tumour laterality, corresponding to patient ID query in search bar
            patient_data <- patient_for_analysis()
            patient_data$LATERALITY <- as.character(patient_data$LATERALITY)
            patient_data$LATERALITY
        })
        
        output$patient_cellularity <- renderText({
          # display patient tumour cellularity, corresponding to patient ID query in search bar
            patient_data <- patient_for_analysis()
            patient_data$CELLULARITY <- as.character(patient_data$CELLULARITY)
            patient_data$CELLULARITY 
        })
        
        output$patient_npi <- renderText({
          # display patient NPI, corresponding to patient ID query in search bar
            patient_data <- patient_for_analysis()
            patient_data$NPI 
        })
        
        output$patient_lymph_nodes <- renderText({
          # display patient's number of lymph nodes examined positive for cancer, corresponding to patient ID query in search bar
          patient_data <- patient_for_analysis()
          patient_data$LYMPH_NODES_EXAMINED_POSITIVE
        })
        
        
        
        ### Clinical Data: cohort-level
        
        # Cohort Overview:
      
        output$total_patients <- renderValueBox({
          # display total number of patients in a value box
          valueBox(nrow(clinical), h4("Total Number of Patients"), icon = icon("user-friends"),
                   color = "aqua"
          )
        })
        
        
        output$os_time_plot <- renderPlotly({
          # Overall survival time plot (histogram)
          os_time_plot <- ggplot(clinical, aes(x=clinical$OS_MONTHS))+
            geom_histogram(color="darkblue", fill="aquamarine3", binwidth = 12)+
            scale_y_continuous(name = "Patients", breaks = seq(0,100, by=25), labels = seq(0,100, by=25))+
            xlab("Overall Survival Time (Months)")+
            ylab("Number of patients")+
            theme_classic() +
            theme(
              panel.grid.major.y = element_line(size = 0.25, 
                                                linetype = 'solid',
                                                colour = "gray"),
              axis.title = element_text(size =14),
              axis.ticks = element_blank(),
              axis.line = element_blank(), 
              axis.text.y = element_text(size = 12), 
              axis.text.x = element_text(size = 12), # adjust font size
              legend.title = element_blank(), 
              legend.position = "bottom")
          
          ggplotly(os_time_plot)
          
        })
        output$age_plot <- renderPlotly({
          # Age at diagnosis (cohort-level) plot (histogram)
          age_plot <- ggplot(clinical, aes(x=clinical$AGE_AT_DIAGNOSIS))+
            geom_histogram(color="darkblue", fill="lightblue", binwidth = 4)+
            xlab("Age at Diagnosis")+
            ylab("Number of patients")+
            theme_classic() +
            theme(
              panel.grid.major.y = element_line(size = 0.25, 
                                                linetype = 'solid',
                                                colour = "gray"),
              axis.title = element_text(size =14),
              axis.ticks = element_blank(),
              axis.line = element_blank(), 
              axis.text.y = element_text(size = 12), 
              axis.text.x = element_text(size = 12), # adjust font size
              legend.title = element_blank(), 
              legend.position = "bottom"
            )
          
          ggplotly(age_plot)
            
        })
        
        menopause <- clinical %>% # prepare data for menopause plot
          group_by(INFERRED_MENOPAUSAL_STATE) %>%
          dplyr::summarise(Patients = n())
        
        output$menopause_plot <- renderPlotly({
          # Bar Chart of inferred menopausal state at cohort-level
          menopause_plot<- ggplot(data = menopause, mapping = aes(x = reorder(INFERRED_MENOPAUSAL_STATE, Patients), y= Patients, fill = INFERRED_MENOPAUSAL_STATE)) +
            geom_col(col = "#FFFFFF", width = 1, size = 4)+
            scale_fill_manual(values = c("deepskyblue3","darkslategray3"),
                              name = NULL)+
            scale_y_continuous(name = "Patients", breaks = seq(0,1250, by=250), labels = seq(0,1250, by=250))+
            xlab("Inferred Menopausal State")+
            theme_classic() +
            theme(
              panel.grid.major.y = element_line(size = 0.25, 
                                                linetype = 'solid',
                                                colour = "gray"),
              axis.title = element_text(size =14),
              axis.ticks = element_blank(),
              axis.line = element_blank(), 
              axis.text.y = element_text(size = 12), 
              axis.text.x = element_text(size = 12), # adjust font size
              legend.title = element_blank(), 
              legend.position = "none"
            )
          
          ggplotly(menopause_plot)
        })

        output$os_min_cohort<- renderText({
          # minimum overall survival time
          min(clinical$OS_MONTHS)
        })
        
        output$os_median_cohort<- renderText({
          # median overall survival time
          median(clinical$OS_MONTHS)
        })
        
        output$os_max_cohort<- renderText({
          # maximum overall survival time
          max(clinical$OS_MONTHS)
        })
        
        os_status_alive <- filter(clinical, clinical$OS_STATUS == "0:LIVING") # filter to only contain patients that are living
        
        output$os_alive_cohort <- renderText({
          # number of patients living
          nrow(os_status_alive)
        })
        
        output$os_alive_box <- renderValueBox({
          # number of patients living displayed in value box
          valueBox(nrow(os_status_alive), h4("Patients Living"), icon = icon("procedures"),
            color = "aqua"
          )
        })
        
        os_status_deceased <- filter(clinical, clinical$OS_STATUS == "1:DECEASED") # filter to only contain patients that are deceased
        
        output$os_deceased_cohort <- renderText({
          # number of patients deceased
          nrow(os_status_deceased)
        })
        
        output$os_deceased_box <- renderValueBox({
          # number of patients deceased displayed in value box
          valueBox(nrow(os_status_deceased), h4("Patients Deceased"), icon = icon("hospital"),
                   color = "aqua"
          )
        })
        output$os_vital_living <- renderText({
          # cohort-level vital status: living
          os_living <- filter(clinical, clinical$VITAL_STATUS == "Living")
          nrow(os_living)
        })
        
        output$os_vital_died_disease <- renderText({
          # cohort-level vital status: died of disease
          os_died_diease <- filter(clinical, clinical$VITAL_STATUS == "Died of Disease")
          nrow(os_died_diease)
        })
        
        output$os_vital_died_other <- renderText({
          # cohort-level vital status: died of other causes
          os_died_other <- filter(clinical, clinical$VITAL_STATUS == "Died of Other Causes")
          nrow(os_died_other)
        })
        
        ## Tumour info: cohort-level

        lymph_nodes_summary <- clinical %>% 
            group_by(LYMPH_NODES_EXAMINED_POSITIVE) %>%
            summarise(Samples = n()) # data wrangling for number of lymph nodes examined positive for entire cohort
        
        output$lymph_nodes_table<- DT::renderDataTable({
          # display in table
            lymph_nodes_summary %>% datatable(colnames = c("Lymphs Nodes", "Patients"), options = list(lengthMenu = c(6, 12,24,28), pageLength = 6))
            
        })
        
        cellularity_data <- clinical %>%
            group_by(CELLULARITY) %>%
            summarise(Patients = n()) # Cellularity info (cohort-level)
        
        output$cellularity_plot <- renderPlotly({
          # bar chart of cellularity info (cohort-level)
          cellularity_plot <- ggplot(data = cellularity_data , mapping = aes(x = reorder(CELLULARITY, Patients), y= Patients, fill = CELLULARITY)) +
            geom_col(col = "#FFFFFF", width = 1, size = 4)+
            scale_fill_manual(values = c("dodgerblue4", "steelblue1", "steelblue"),
                              name = NULL)+
            scale_y_continuous(name = "Patients", breaks = seq(0,800, by=200), labels = seq(0,800, by=200))+
            theme_classic() +
            theme(
              panel.grid.major.y = element_line(size = 0.25, 
                                                linetype = 'solid',
                                                colour = "gray"),
              axis.title.y = element_text(size =14),
              axis.title.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank(), 
              axis.text.y = element_text(size = 12), 
              axis.text.x = element_text(size = 14), # adjust font size
              legend.title = element_blank(), 
              legend.position = "none"
            )
          
          ggplotly(cellularity_plot)
        })
       
        output$NPI_plot <- renderPlotly({
          # histogram of NPI values (cohort-level)
          npi_plot <- ggplot(clinical, aes(x=clinical$NPI))+
            geom_histogram(color="darkblue", fill="palegreen2", binwidth = 1)+
            scale_y_continuous(name = "Patients", breaks = seq(0,550, by=50), labels = seq(0,550, by=50))+
            xlab("NPI")+
            ylab("Number of patients")+
            theme_classic() +
            theme(
              panel.grid.major.y = element_line(size = 0.25, 
                                                linetype = 'solid',
                                                colour = "gray"),
              axis.title = element_text(size =14),
              axis.ticks = element_blank(),
              axis.line = element_blank(), 
              axis.text.y = element_text(size = 12), 
              axis.text.x = element_text(size = 12), # adjust font size
              legend.title = element_blank(), 
              legend.position = "bottom")
          
          ggplotly(npi_plot)
        })
        
        ER_status_IHC <- clinical %>%
          group_by(ER_IHC) %>%
          summarise(Patients =n()) # ER status measured by IHC (cohort-level)
        
        output$ER_IHC_plot <- renderPlotly({
          # bar chart of ER status measured by IHC (cohort-level)
          ER_status_plot <- ggplot(data = ER_status_IHC , mapping = aes(x = ER_IHC, y= Patients, fill = ER_IHC)) +
            geom_col(col = "#FFFFFF", width = 1, size = 4)+
            scale_fill_manual(values = c("violetred", "slateblue3"),
                              name = NULL)+
            scale_y_continuous(name = "Patients", breaks = seq(0,1200, by=200), labels = seq(0,1200, by=200))+
            theme_classic() +
            theme(
              panel.grid.major.y = element_line(size = 0.25, 
                                                linetype = 'solid',
                                                colour = "gray"),
              axis.title.y = element_text(size =14),
              axis.title.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank(), 
              axis.text.y = element_text(size = 12), 
              axis.text.x = element_text(size = 14), # adjust font size
              legend.title = element_blank(), 
              legend.position = "none"
            )
          
          ggplotly(ER_status_plot)
          
        })
        
        HER2_snp6 <- clinical %>%
          group_by(HER2_SNP6) %>%
          summarise(Patients =n()) # HER2 status measured by SNP6
        
        output$HER2_SNP6_plot <- renderPlotly({
          # bar chart of HER2 status measured by SNP6 (cohort-level)
          her2_snp6_plot<-  ggplot(data =  HER2_snp6 , mapping = aes(x = HER2_SNP6, y= Patients, fill = HER2_SNP6)) +
            geom_col(col = "#FFFFFF", width = 1, size = 4)+
            scale_fill_manual(values = c("limegreen", "palevioletred", "lightskyblue", "steelblue"),
                              name = NULL)+
            scale_y_continuous(name = "Patients", breaks = seq(0,1200, by=200), labels = seq(0,1200, by=200))+
            theme_classic() +
            theme(
              panel.grid.major.y = element_line(size = 0.25, 
                                                linetype = 'solid',
                                                colour = "gray"),
              axis.title.y = element_text(size =14),
              axis.title.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank(), 
              axis.text.y = element_text(size = 12), 
              axis.text.x = element_text(size = 14), # adjust font size
              legend.title = element_blank(), 
              legend.position = "none"
            )
          
          ggplotly(her2_snp6_plot)
        })
        
        
        histological <-  clinical %>%
          group_by(HISTOLOGICAL_SUBTYPE) %>%
          summarise(Patients =n()) # histological subtype (cohort-level)
        
        output$histological_table <- DT::renderDataTable({
          # display histological subtype info from entire cohort in datatable
          histological %>% datatable(colnames = c("Tumor's histologic subtype", "Number of Patients"), options = list(lengthMenu = c(5, 7), pageLength = 7))
        })
        # Treatment info: cohort-level
        chemo <- clinical %>%
            group_by(CHEMOTHERAPY) %>%
            summarise(Patients =n())
        
        output$chemo_plot <- renderPlotly({
          # Bar chart of cohort level chemotherapy info
          chemo_plot <- ggplot(data = chemo, mapping = aes(x = CHEMOTHERAPY, y= Patients, fill = CHEMOTHERAPY )) +
            geom_col(col = "#FFFFFF", width = 1, size = 4)+
            scale_fill_manual(values = c("YES" = "steelblue4","NO"="skyblue"),
                              name = NULL)+
            scale_y_continuous(name = "Patients", breaks = seq(0,1250, by=250), labels = seq(0,1250, by=250))+
            xlab("Chemotherapy")+
            theme_classic() +
            theme(
              panel.grid.major.y = element_line(size = 0.25, 
                                                linetype = 'solid',
                                                colour = "gray"),
              axis.title = element_text(size =14),
              axis.ticks = element_blank(),
              axis.line = element_blank(), 
              axis.text.y = element_text(size = 12), 
              axis.text.x = element_text(size = 12), # adjust font size
              legend.title = element_blank(), 
              legend.position = "none"
            )
          
          ggplotly(chemo_plot)
        })
        
        hormone_therapy <- clinical %>%
            group_by(HORMONE_THERAPY) %>%
            summarise(Patients =n())
        
        output$hormone_therapy_plot <- renderPlotly({
          # Bar chart of cohort level hormone therapy info
          hormone_therapy_plot <- ggplot(data = hormone_therapy, mapping = aes(x = HORMONE_THERAPY, y= Patients, fill = HORMONE_THERAPY )) +
            geom_col(col = "#FFFFFF", width = 1, size = 4)+
            scale_fill_manual(values = c("YES" = "steelblue4","NO"="skyblue"),
                              name = NULL)+
            scale_y_continuous(name = "Patients", breaks = seq(0,1250, by=250), labels = seq(0,1250, by=250))+
            xlab("Hormone Therapy")+
            theme_classic() +
            theme(
              panel.grid.major.y = element_line(size = 0.25, 
                                                linetype = 'solid',
                                                colour = "gray"),
              axis.title = element_text(size =14),
              axis.ticks = element_blank(),
              axis.line = element_blank(), 
              axis.text.y = element_text(size = 12), 
              axis.text.x = element_text(size = 12), # adjust font size
              legend.title = element_blank(), 
              legend.position = "none"
            )
          
          ggplotly(hormone_therapy_plot)
        })
        
        radio_therapy <- clinical %>%
            group_by(RADIO_THERAPY) %>%
            summarise(Patients =n())
        
        output$radio_therapy_plot <- renderPlotly({
          # Bar chart of cohort level radiotherapy info
          radio_therapy_plot <- ggplot(data = radio_therapy, mapping = aes(x = RADIO_THERAPY, y= Patients, fill = RADIO_THERAPY )) +
            geom_col(col = "#FFFFFF", width = 1, size = 4)+
            scale_fill_manual(values = c("YES" = "steelblue4","NO"="skyblue"),
                              name = NULL)+
            scale_y_continuous(name = "Patients", breaks = seq(0,1250, by=250), labels = seq(0,1250, by=250))+
            xlab("Radio Therapy")+
            theme_classic() +
            theme(
              panel.grid.major.y = element_line(size = 0.25, 
                                                linetype = 'solid',
                                                colour = "gray"),
              axis.title = element_text(size =14),
              axis.ticks = element_blank(),
              axis.line = element_blank(), 
              axis.text.y = element_text(size = 12), 
              axis.text.x = element_text(size = 12), # adjust font size
              legend.title = element_blank(), 
              legend.position = "none"
            )
          
          ggplotly(radio_therapy_plot)
        })
        
        breast_surgery <- clinical %>%
            group_by(BREAST_SURGERY) %>%
            summarise(Patients =n())
        
        output$breast_surgery_plot <- renderPlotly({
          # Bar chart of cohort level breast surgery info
          breast_surgery_plot <- ggplot(data = breast_surgery, mapping = aes(x = BREAST_SURGERY, y= Patients, fill = BREAST_SURGERY )) +
            geom_col(col = "#FFFFFF", width = 1, size = 4)+
            scale_fill_manual(values = c("plum2","plum4"),
                              name = NULL)+
            scale_y_continuous(name = "Patients", breaks = seq(0,1000, by=200), labels = seq(0,1000, by=200))+
            xlab("Breast Surgery Type")+
            theme_classic() +
            theme(
              panel.grid.major.y = element_line(size = 0.25, 
                                                linetype = 'solid',
                                                colour = "gray"),
              axis.title = element_text(size =14),
              axis.ticks = element_blank(),
              axis.line = element_blank(), 
              axis.text.y = element_text(size = 12), 
              axis.text.x = element_text(size = 12), # adjust font size
              legend.title = element_blank(), 
              legend.position = "none"
            )
          
          ggplotly(breast_surgery_plot)
        })
        
        ## Breast cancer classification info: cohort-level
        claudin_subtype <- clinical %>%
            group_by(CLAUDIN_SUBTYPE) %>%
            summarise(Patients =n())
        
        output$claudin_subtype_plot <- renderPlotly({
          # Bar chart of cohort level claudin subtype
          claudin_plot <-ggplot(data = claudin_subtype, mapping = aes(x = reorder(CLAUDIN_SUBTYPE, Patients), y= Patients, fill = CLAUDIN_SUBTYPE )) +
            geom_col(col = "#FFFFFF", width = 1, size = 4)+
            scale_y_continuous(name = "Patients", breaks = seq(0,600, by=100), labels = seq(0,600, by=100))+
            xlab("Claudin Subtype")+
            theme_classic() +
            theme(
              panel.grid.major.y = element_line(size = 0.25, 
                                                linetype = 'solid',
                                                colour = "gray"),
              axis.title = element_text(size =14),
              axis.ticks = element_blank(),
              axis.line = element_blank(), 
              axis.text.y = element_text(size = 12), 
              axis.text.x = element_text(size = 12), # adjust font size
              legend.title = element_blank(), 
              legend.position = "none")
          
          ggplotly(claudin_plot)
        })
        
        int_cluster <- clinical %>%
            group_by(INTCLUST) %>%
            summarise(Patients =n())
        output$int_cluster_plot<- renderPlotly({
          # Bar chart of cohort leve integrative cluster
          intclust_plot <-ggplot(data = int_cluster, mapping = aes(x = reorder(INTCLUST, Patients), y= Patients)) +
            geom_col(color="#FFFFFF", fill="aquamarine3", width = 1, size = 1)+
            scale_y_continuous(name = "Patients", breaks = seq(0,250, by=50), labels = seq(0,250, by=50))+
            xlab("Integrative Cluster")+
            theme_classic() +
            theme(
              panel.grid.major.y = element_line(size = 0.25, 
                                                linetype = 'solid',
                                                colour = "gray"),
              axis.title = element_text(size =14),
              axis.ticks = element_blank(),
              axis.line = element_blank(), 
              axis.text.y = element_text(size = 12), 
              axis.text.x = element_text(size = 12), # adjust font size
              legend.title = element_blank(), 
              legend.position = "none")
          
          ggplotly(intclust_plot)
            
        })
        three_gene <- clinical %>%
            group_by(THREEGENE) %>%
            summarise(Patients =n())
        
        output$three_gene_plot <- renderPlotly({
          # Bar chart of cohort level three-gene classifier subtype
          three_gene_plot <- ggplot(data = three_gene, mapping = aes(x = reorder(THREEGENE, Patients), y= Patients, fill = THREEGENE)) +
            geom_col(color="#FFFFFF", width = 1, size = 4)+
            scale_fill_manual(values = c("cornflowerblue","blueviolet", "chartreuse2", "hotpink3"),
                              name = NULL)+
            scale_y_continuous(name = "Patients", breaks = seq(0,600, by=100), labels = seq(0,600, by=100))+
            theme_classic() +
            coord_flip()+
            theme(
              panel.grid.major.x = element_line(size = 0.25, 
                                                linetype = 'solid',
                                                colour = "gray"),
              axis.title.x = element_text(size =14),
              axis.title.y = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank(), 
              axis.text.y = element_text(size = 12), 
              axis.text.x = element_text(size= 11),
              legend.title = element_blank(), 
              legend.position = "none")
          
          ggplotly(three_gene_plot)
        })
        
        laterality <- clinical %>%
          group_by(LATERALITY) %>%
          summarise(Patients =n())
        
        output$laterality_plot <- renderPlotly({
          # Bar chart of cohort level primary tumour laterality
          laterality_plot <- ggplot(data = laterality, mapping = aes(x = LATERALITY, y= Patients, fill = LATERALITY)) +
            geom_col(col = "#FFFFFF", width = 1, size = 4)+
            scale_fill_manual(values = c("springgreen3","skyblue2"),
                              name = NULL)+
            scale_y_continuous(name = "Patients", breaks = seq(0,800, by=200), labels = seq(0,800, by=200))+
            xlab("Primary Tumour Laterality")+
            theme_classic() +
            theme(
              panel.grid.major.y = element_line(size = 0.25, 
                                                linetype = 'solid',
                                                colour = "gray"),
              axis.title = element_text(size =14),
              axis.ticks = element_blank(),
              axis.line = element_blank(), 
              axis.text.y = element_text(size = 12), 
              axis.text.x = element_text(size = 14), # adjust font size
              legend.title = element_blank(), 
              legend.position = "none"
            )
          
          ggplotly(laterality_plot)
        })
        
        ### Genomic Data
        # Enrichment Analysis - gprofiler2
        gostres <- gost(query = driver, 
                          organism = "hsapiens", ordered_query = FALSE, 
                          multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                          measure_underrepresentation = FALSE, evcodes = FALSE, 
                          user_threshold = 0.05, correction_method = "g_SCS", 
                          domain_scope = "annotated", custom_bg = NULL, 
                          numeric_ns = "", sources = NULL, as_short_link = FALSE)
          
        output$g_profiler_plot <- renderPlotly({
          # plot functional enrichment analysis result
          gostplot(gostres, capped = TRUE, interactive = TRUE)})

        
        # Classify gene regulation status: >median is up-regulated genes and <median is down regulated genes
        event_rna <- apply(gene_expression,1, function(x) ifelse(x > median(x),"up","down"))
        event_rna <- as.data.frame(event_rna) #should be as data frame
        
        genomic_cohort_data <- event_rna
        
        genes_cohort <- colnames(genomic_cohort_data) 
        
        ### Patient-level genomic
        patient_genomic_data <- t(genomic_cohort_data) # transpose dataframe
        patient_genomic_data  <- as.data.frame(patient_genomic_data) # make as dataframe
        
        
        ## Up-regulated genes
        patient_upregulated_genes <- reactive({
          upregulated <- select(patient_genomic_data, one_of(input$searchpatient))
          upregulated %>%
            filter(upregulated[,1] == "up")
        })
        
        output$upregulated_genes_table<- DT::renderDataTable({
          # Data table of upregulated genes corresponding to patient ID entered in search bar
          req(input$searchpatient)
          if(!(input$searchpatient %in% colnames(patient_genomic_data))){
            return(NULL)
          }
          up_genes <- patient_upregulated_genes()
          up_genes %>% datatable()
          
        })
        
        # down-regulated genes
        patient_downregulated_genes <- reactive({
          # Data table of downregulated genes corresponding to patient ID entered in search ba
          downregulated <- select(patient_genomic_data, one_of(input$searchpatient))
          downregulated %>%
            filter(downregulated[,1] == "down")
        })
        
        output$downregulated_genes_table<- DT::renderDataTable({
          req(input$searchpatient)
          if(!(input$searchpatient %in% colnames(patient_genomic_data))){
            return(NULL)
          }
          down_genes <- patient_downregulated_genes()
          down_genes %>% datatable()
          
        })
        
        
        
        ## Identify driver genes associated with survival
        # create new column ‘status’ with survival event as binary
        clinical_exp$status <- clinical_exp$OS_STATUS
        clinical_exp <- mutate(clinical_exp, status = recode(status, '0:LIVING'= 0, "1:DECEASED" = 1)) # Recode status
        
        #add time and event columns of clinical_exp to event_rna
        event_rna <- cbind(event_rna,clinical_exp$OS_MONTHS) #time
        event_rna <- cbind(event_rna,clinical_exp$status) #event
        colnames(event_rna)[36:37] <- c("time", "event") # rename columnes
        
        # Cox PH model
        set.seed(420)
        surv_formulas <- sapply(driver,
                                function(x) as.formula(paste('Surv(time, event)~', as.factor(x))))
        surv_models <- lapply(surv_formulas, function(x){coxph(x, data = event_rna)})
        # results
        surv_results <- lapply(surv_models,
                               function(x){ 
                                 x <- summary(x)
                                 p.value<-signif(x$wald["pvalue"], digits=2)
                                 wald.test<-signif(x$wald["test"], digits=2)
                                 beta<-signif(x$coef[1], digits=2);#coeficient beta
                                 HR <-signif(x$coef[2], digits=2);#exp(beta)
                                 HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                                 HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                                 HR <- paste0(HR, " (", 
                                              HR.confint.lower, "-", HR.confint.upper, ")")
                                 res<-c(beta, HR, wald.test, p.value)
                                 names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                               "p.value")
                                 return(res)
                               })
        res <- t(as.data.frame(surv_results, check.names = FALSE)) # results from Cox PH model
        res <- as.data.frame(res) # Make dataframe of results
        
        output$surv_gene_results_table<- DT::renderDataTable({
          # display results of Cox PH in datatable
          res %>% datatable(options = list(lengthMenu = c(10, 20, 30, 35), pageLength = 10)) %>%
            formatStyle(
              'p.value',
              target = 'row', # highlight rows with p value <0.05 in green
              backgroundColor = styleInterval(0.05, c('limegreen', 'white'))
            )
            
          
        })
        
        ## driver genes significantly associated with survival
        res$p.value <- as.numeric(as.character(res$p.value)) # change from factor so filter can be applied
        surv_sig_genes <- filter(res, res$p.value <= 0.05) # filter dataframe to only contain genes with p value <0.05
        
        sig_genes <- rownames(surv_sig_genes) # create vector of names of genes identified by Cox PH to be associated with survival
        # Enrichment Analysis of the genes identified to be significantly associated with survival
        gostres_surv <- gost(query = sig_genes, 
                        organism = "hsapiens", ordered_query = FALSE, 
                        multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                        measure_underrepresentation = FALSE, evcodes = FALSE, 
                        user_threshold = 0.05, correction_method = "g_SCS", 
                        domain_scope = "annotated", custom_bg = NULL, 
                        numeric_ns = "", sources = NULL, as_short_link = FALSE)
        
        output$g_profiler_surv_plot <- renderPlotly({
          # display results of functional enrichment analysis of genes significantly associated with survival in interactive plot
          gostplot(gostres_surv, capped = TRUE, interactive = TRUE)})
        
        ### 10-year Survival Prediction
        
        ## Clinical-only Prediction

        output$selected_clin <- renderPrint({ input$clinical_predictors_ten }) # Print the predictors selected by user 
        
        observeEvent(input$clinical_ten_button, {
          # Display message if user has not specified input
          if (is.null(input$clinical_predictors_ten)) {
            shinyjs::show("message1")
          } else {
            shinyjs::hide("message1")}
          
          req(input$clinical_predictors_ten)
          set.seed(103) # set seed 
          # Shuffle row indices: rows
          rows_ten_combined <- sample(nrow(combined_data))
          # Randomly order data
          combined_data <- combined_data[rows_ten_combined, ]
          
          ten_year_surv_clinical <- select(combined_data, one_of(input$clinical_predictors_ten, "ten_year")) # filter data to contain predictors selected by user
          
          set.seed(111)
          train <- sample(nrow(ten_year_surv_clinical), 0.8*nrow(ten_year_surv_clinical), replace = FALSE) # split data into 80/20 (training/test)
          TrainSet <- ten_year_surv_clinical[train,] # training data
          ValidSet <- ten_year_surv_clinical[-train,] # test data
          
          
          set.seed(121)
          clinical_model <- randomForest(ten_year ~ ., data = TrainSet, importance = TRUE) # RSF model
          
          rf_prediction_clinical <- predict(clinical_model, ValidSet, type = "prob") # predict with RSF model using testing data 
          ROC_rf_clin <- roc(ValidSet$ten_year, rf_prediction_clinical[,2]) # ROC curve for prediction on test set 
          ROC_rf_auc_clin <- auc(ROC_rf_clin) # AUC for prediction on test set
          
          output$clinical_model_summary <- renderPrint({ clinical_model }) # print summary of the generated model
          
          
          output$auc_clin <- renderText({
            # print AUC 
            ROC_rf_auc_clin
          })
          
          output$roc_clinical_plot <- renderPlot({
            # plot ROC
            plot(ROC_rf_clin, col = "aquamarine3", main = "ROC For Random Forest (Clinical Predictors only)")
          })
          
          # variable importance
          varimp_clinical <- importance(clinical_model)
          varimp_clinical <- as.data.frame(varimp_clinical, rownames = TRUE) # make as dataframe

          output$var_imp_clinical_table <- DT::renderDataTable({
            # display variable importance in datatable
            varimp_clinical %>% round(digits = 3) %>%
              datatable(rownames = TRUE, options = list(lengthMenu = c(5, 10, 16), pageLength = 16))
            
          })
          
        }) 
        
        ## Genomic-only prediction
        output$selected_genes <- renderPrint({ input$genomic_predictors_ten }) # Print the predictors selected by user 
        
        observeEvent(input$genomic_ten_button, {
          # Display message if user has not specified input
          if (is.null(input$genomic_predictors_ten)) {
            shinyjs::show("message2")
          } else {
            shinyjs::hide("message2")}
          
          req(input$genomic_predictors_ten)
          # need to shuffle dataframe before we can perform prediction
          # Set seed
          set.seed(42)
          # Shuffle row indices: rows
          rows_ten_combined <- sample(nrow(combined_data))
          # Randomly order data
          combined_data <- combined_data[rows_ten_combined, ]
          
          ten_year_genomic <- select(combined_data, one_of(input$genomic_predictors_ten, "ten_year"))  # filter data to contain predictors selected by user
          
          set.seed(100)
          train <- sample(nrow(ten_year_genomic), 0.8*nrow(ten_year_genomic), replace = FALSE) # split data into 80/20 (training/test)
          TrainSet <- ten_year_genomic[train,] # training set
          ValidSet <- ten_year_genomic[-train,] # testing set
          
          
          set.seed(101)
          model_genomic <- randomForest(ten_year ~ ., data = TrainSet, importance = TRUE) # RSF model
          
          rf_prediction_genomic <- predict(model_genomic, ValidSet, type = "prob") # predict with RSF model using testing data 
          ROC_rf_g <- roc(ValidSet$ten_year, rf_prediction_genomic[,2]) # ROC curve for prediction on test set 
          ROC_rf_auc_g <- auc(ROC_rf_g) # AUC for prediction on test set
          
          output$genomic_model_summary <- renderPrint({ model_genomic }) # print summary of the generated model
          
          output$auc_genomic <- renderText({
            # print AUC 
            ROC_rf_auc_g
          })
          
          output$roc_genomic_plot <- renderPlot({
            # plot ROC
            plot(ROC_rf_g, col = "cornflowerblue", main = "ROC For Random Forest (Genomic Predictors only)")
          })
          
          # variable importance
          varimp_genomic <- importance(model_genomic)
          varimp_genomic <- as.data.frame(varimp_genomic) # make as dataframe
          
          #varimp_genomic <- varimp_genomic %>% mutate_if(is.numeric, round, digits=3) # round to 3 decimal places
          
          output$var_imp_genomic_table <- DT::renderDataTable({
            # display variable importance in datatable
            varimp_genomic %>% round(digits = 3) %>%
              datatable(options = list(lengthMenu = c(5, 10, 16), pageLength = 16))
            
          })
        }) 
        
        ## Combined prediction
        output$selected_combined <- renderPrint({ input$combined_predictors_ten }) # Print the predictors selected by user 
        
        observeEvent(input$combined_ten_button, {
          # Display message if user has not specified input
          if (is.null(input$combined_predictors_ten)) {
            shinyjs::show("message3")
          } else {
            shinyjs::hide("message3")}
          
          req(input$combined_predictors_ten)
          # need to shuffle dataframe before we can perform prediction
          # Set seed
          set.seed(42)
          # Shuffle row indices: rows
          rows_ten_combined <- sample(nrow(combined_data))
          # Randomly order data
          combined_data <- combined_data[rows_ten_combined, ]
          
          ten_year_combined <- select(combined_data, one_of(input$combined_predictors_ten, "ten_year"))  # filter data to contain predictors selected by user
          
          set.seed(100)
          train <- sample(nrow(ten_year_combined), 0.8*nrow(ten_year_combined), replace = FALSE) # split data into 80/20 (training/test)
          TrainSet <- ten_year_combined[train,] # training set
          ValidSet <- ten_year_combined[-train,] # testing set
          
          
          set.seed(201)
          model_combined <- randomForest(ten_year ~ ., data = TrainSet, importance = TRUE) # RSF model
          
          rf_prediction_combined <- predict(model_combined, ValidSet, type = "prob") # predict with RSF model using testing data 
          ROC_rf_combined <- roc(ValidSet$ten_year, rf_prediction_combined[,2]) # ROC curve for prediction on test set 
          ROC_rf_auc_combined <- auc(ROC_rf_combined) # AUC for prediction on test set
          
          output$combined_model_summary <- renderPrint({ model_combined }) # print summary of the generated model
          
          output$auc_combined <- renderText({
            # print AUC 
            ROC_rf_auc_combined
          })
          
          output$roc_combined_plot <- renderPlot({
            # plot ROC
            plot(ROC_rf_combined, col = "cadetblue4", main = "ROC For Random Forest (Combined Predictors)")
          })
          
          # variable importance
          varimp_combined <- importance(model_combined)
          varimp_combined <- as.data.frame(varimp_combined) # make as dataframe
          
          output$var_imp_combined_table <- DT::renderDataTable({
            # display variable importance in databtable
            varimp_combined %>% round(digits = 3) %>%
              datatable(options = list(lengthMenu = c(5, 10, 15, 20, 25, 31), pageLength = 15))
        }) 
        })
    }
)
