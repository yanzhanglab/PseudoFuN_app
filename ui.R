library('shiny')
library('networkD3')
source('helpers.R')
library('shinyWidgets')
library(markdown)
library(colourpicker)

# Define UI for app that draws a histogram ----


navbarPage(  
  # App title ----
  # titlePanel("PseudoFuN DB Search"),
  title=div(a(img(src="http://icons.iconarchive.com/icons/graphicloads/100-flat/256/zoom-search-2-icon.png",
                  height = 30,
                  style = "margin:0px 0px; padding-bottom: 5px"),
              "PseudoFuN DB Search", href="")
           ),
  
  
  tabPanel("Search Engine",
    # Sidebar layout with input and output definitions ----
    sidebarLayout(
      
      # Sidebar panel for inputs ----
      sidebarPanel(width = 3,
        
        # Input: Slider for the number of bins ----
        h4("Choose a database", style="color: STEELBLUE"),
        selectInput("db",
                    label = "",
                    choices = c("BlastDB",
                                "CUDAlign18",
                                "CUDAlign54",
                                "CUDAlign135",
                                "CUDAlign198"),
                    selected = "BlastDB"),
        
        textInput("gene",
                  h4("Enter a gene", style="color: STEELBLUE"),
                  value = "PTEN"),
        
        h4("GO Analysis", style="color: STEELBLUE"),
        awesomeRadio("go",
                    label = "",
                    choices = list("Do Not Run GO Analysis"=0,
                                   "Run GO Analysis: Biological Process"=1,
                                   "Run GO Analysis: Molecular Function"=2,
                                   "Run GO Analysis: Cellular Component"=3),
                    status = "primary"),
        conditionalPanel(condition="input.go > 0",
                         helpText("Fisher classic will be run by default. KS Classic and KS elim are optional and may cost longer time (few minutes)."),
                         prettyCheckbox(inputId = "run.ks", label = "Run KS Classic", value = F, icon = icon("check")),
                         prettyCheckbox(inputId = "run.ks.elim", label = "Run KS elim", value = F, icon = icon("check")),
                         prettyCheckbox(inputId = "inc0", label = "Include GO terms without any assigned genes", value = FALSE, icon = icon("check"))
                         ),
        # checkboxInput("inc0", "Include GO terms without any assigned genes", value = FALSE)
        
        actionButton("action1", "Confirm and Run", style="color: WHITE; background-color: DODGERBLUE")
      ),
      
      # Main panel for displaying outputs ----
      mainPanel(
        tabsetPanel(id = "tabs",
          tabPanel("GO Analysis",
                   h4("GO Analysis (all networks)", style="color: STEELBLUE"),
                   h4("To view DB search results please select Network tab", style="color: STEELBLUE; font-size: 14px"),
                   DT::dataTableOutput("GOtable"),
                   tags$head(tags$script('// disable download at startup.
                                         $(document).ready(function() {
                                         $("#download_go").attr("disabled", "true").attr("onclick", "return false;");
                                         Shiny.addCustomMessageHandler("download_go", function(message) {
                                         $("#download_go").removeAttr("disabled").removeAttr("onclick");
                                         });
                                         })')),
                   downloadButton('download_go', 'Download GO Analysis Results (CSV)'))
                   )
        )
      )
  ), # end of tabPanel "Search Engine"
  
  tabPanel("Circos Plot",
           sidebarLayout(
             position = "left",
             sidebarPanel(
               width = 3,
               h4("Change Plot Setting", style="color: STEELBLUE"),
               helpText("Rendering Circos Plot could cost up to 1 minutes. Please be patient."),
               prettyCheckbox(inputId = "circos_param_genelink", label = "Show Pair-wise Links", value = T, icon = icon("check")),
               prettyCheckbox(inputId = "circos_param_genesymbol", label = "Show Gene Symbols", value = T, icon = icon("check")),
               fluidRow(
                 column(6, numericInput(inputId="font.scale", label="Font Scale:", value = 1, min = 0.2, max=10, step = 0.1)),
                 column(6, numericInput(inputId="link.width", label="Link Width:", value = 2, min = 1, max=100, step = 0.1))
               ),
               colourInput(inputId="color.picker",label="Choose Link Color:",value="pink",
                           showColour = "both",palette = "square"),
               actionButton("action2", "Refresh", style="color: WHITE; background-color: DODGERBLUE")
               
               
             ),
             mainPanel(
               h2("Please be patient plots may take a few seconds to render.", style="color: STEELBLUE; font-size: 14px"),
               uiOutput("Circos_plot_gene_db_name"),
               h2("hg38:", style="color: STEELBLUE; font-size: 22px"),
               plotOutput("circos_plot_component_hg38", width = 800, height = 600),
               h2("hg19:", style="color: STEELBLUE; font-size: 22px"),
               plotOutput("circos_plot_component_hg19", width = 800, height = 600)
             ) # end of mainPanel
           ) # end of sidebarLayout
  ),
  tabPanel("TCGA Expression",
           sidebarLayout(
             position = "left",
             sidebarPanel(
               width = 3,
               h4("TCGA Cancer Selection", style="color: STEELBLUE"),
               helpText("Please select the cancer to conduct gene expression analysis on the selected network.
                        If no search has been conducted, please conduct a search using the Search Engine tab.
                        To conduct the TCGA expression analysis on the search results, the TCGA Expression
                        button must be selected under the  Search Engine tab."),
               selectInput(inputId = "TCGA_cancer", label = "TCGA Cancers",
                           choices = list("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM"),
                           selected = "BRCA", multiple = FALSE),
               selectInput(inputId = "miRNA_pred_targ", label = "Number algorithms predicting miRNA targetting of gene",
                           choices = c(0,1,2,3),
                           selected = 0, multiple = FALSE),
               h4("Differential Pseudogene Expression (DPgE)", style="color: STEELBLUE"),
               helpText("The DPgE table is pregenerated using a linear model and the FDR is calculated for each
                        pseudogene in the TCGA cancer chosen above. In this case (since there are only tumor and normal groups), the F value is
                        equivalent to the T statistic in the form (F = T^2). Below a cutoff FDR value can be selected for 
                        the displayed DPgE table."),
               sliderInput(inputId = "DGE_cutoff_value", label = "FDR cutoff: 1e-value",
                           min = 10,max = 200, value = 15,step = 1),
               actionButton("action3", "Confirm", style="color: WHITE; background-color: DODGERBLUE")
               
             ),
             mainPanel(
               h2("TCGA Expression Panel", style="color: STEELBLUE; font-size: 22px"),
               uiOutput("TCGA_Expression_gene_db_cancer_name"),
               h2("Please be patient plots may take a few seconds to render.", style="color: STEELBLUE; font-size: 14px"),
               fluidRow(
                 column(6, h4("Normal Sample Coexpression", style="color: STEELBLUE")),
                 column(6, h4("Tumor Sample Coexpression", style="color: STEELBLUE"))
               ),
               fluidRow(
                 column(6, plotOutput("normal_heatmap", width = 500, height = 400)),
                 column(6, plotOutput("tumor_heatmap", width = 500, height = 400))
               ),
               h4("Differential Expression Tumor vs. Normal", style="color: STEELBLUE"),
               plotOutput("pseudo_boxplot", width = 1200, height = 500),
               h4("Gene and Pseudogene miRNA Associations", style="color: STEELBLUE"),
               plotOutput("correlation_plot", width = 600, height = 600),
               h4("DPgE Table", style="color: STEELBLUE"),
               DT::dataTableOutput("DGEtable"),
               fluidRow(
                column(5, downloadButton('download_TCGA_exp', 'Download Data'))
               )
             ) # end of mainPanel
           ) # end of sidebarLayout
  ),
  tabPanel("Read Me",
           # fluidRow(
           #   column(width = 8, offset = 1,
           includeMarkdown("./README.md")
           #   )
           # )
  ),
  tabPanel("Tutorial",
           h3("Tutorial", style="color: STEELBLUE; padding-bottom: 20px"),
           h4("Google Slides", style="text-align: center; color: STEELBLUE; padding-bottom: 20px"),
           tags$div(
             HTML("<iframe src=\"https://docs.google.com/presentation/d/e/2PACX-1vRqdNS0ZHg8z7VTZqBwleJ2MVuzTKXDWqnKNPrRPl3wiXSNY3YxKyvRNdcmo8Pwk6l2q5uzw2VdoHgv/embed?start=false&loop=false&delayms=3000\" frameborder=\"0\" width=\"960\" height=\"569\" allowfullscreen=\"true\" mozallowfullscreen=\"true\" webkitallowfullscreen=\"true\"></iframe>"),
             style="text-align: center; padding: 20px"
           )
  ),
  tabPanel("About",
           h3("About Us", style="color: STEELBLUE; padding-bottom: 20px"),
           "The Yan Zhang Lab at OSUMC studies statistical and computational methods and their applications to genomic and proteomic research, such as (1) functional and evolutionary impact of structural variations (such as insertions, deletions and retroduplications) in both normal and abnormal populations; (2) association and eQTL analysis of structural variations; (3) integrative analysis of genomic and proteomic data; (4) statistical modeling of biological networks, integrating data of multiple levels.",
           tags$div(
             tags$img(src='images/osumc_logo.png',
                      height="125",
                      alt="OSUMC", class="center", style="padding: 30px"),
             tags$img(src='images/IUSM2.png',
                      height="100",
                      alt="IUSM", class="center", style="padding: 30px"),
             style="text-align: center; padding: 20px"
           ),
           h4("Development Team", style="color: STEELBLUE"),
           tags$ul(
             tags$li("Travis Johnson"),
             tags$li("Eric Franz"),
             tags$li("Zhi Huang")
           ),
           h4("Prof. Yan Zhang's Laboratory", style="color: STEELBLUE"),
           tags$ul(
             tags$li("Yan Zhang"),
             tags$li("Travis Johnson"),
             tags$li("Sihong Li")
           )
  ),
  
  tags$head(
    tags$script(HTML("(function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){(i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
                     m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
                     })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');
                     ga('create', 'UA-113406500-2', 'auto');
                     ga('send', 'pageview');")),
    tags$script('Shiny.addCustomMessageHandler("myCallbackHandler",
                function(typeMessage) {console.log(typeMessage)
                if(typeMessage == "tab_circos"){
                  console.log("Go to Circos Plot panel.");
                $("a:contains(Circos Plot)").click();
                }
                if(typeMessage == "tab_TCGA_Expression"){
                  console.log("Go to TCGA Expression panel.");
                $("a:contains(TCGA Expression)").click();
                }
                });'),
    tags$script(HTML("document.title = 'PseudoFuN DB Search';"))# rename the title by JS
  ),
  
  tags$div(
    p(a("PseudoFuN", href=""), "Version Beta | ", a("OSUMC",href="https://wexnermedical.osu.edu", target="_blank"), " | ", a("IUSM",href="https://medicine.iu.edu/", target="_blank"), style="color: grey; font-size: 12px"), 
    p("Questions and feedback: travis.johnson@osumc.edu | ", a("Report Issue", href="https://github.com/yanzhanglab/PseudoFuN_app/issues", target="_blank"), " | ", a("Github", href="https://github.com/yanzhanglab/PseudoFuN_app", target="_blank"), style="color: grey; font-size: 12px"),
    style="text-align: center; padding-top: 40px"
  ),
  
  tags$div(
    
    # Shiny shows the outer conditionalPanel as long as the document hasn't
    # loaded; the inner rmd_loader is shown by rmd_loader.js as soon as
    # we've been waiting a certain number of ms
    shiny::conditionalPanel(
      "!output.__reactivedoc__",
      tags$div(
        id = "rmd_loader_wrapper",
        tags$div(id = "rmd_loader", style = "display: none",
                 tags$p("Please wait... PseudoFuN is initializing ..."),
                 tags$p("PseudoFuN is a novel database and query tool for homologous PseudoGene and coding Gene (PGG) families.
It supports dynamic search, graphical visualization and functional analysis of pseudogenes and coding genes based on the PGG families. 
                        This work sets a start point for functional analysis of potentially regulatory pseudogenes.
                        
                        This Shiny app supports queries of pseudogene/coding gene names. It is implemented by Travis S Johnson and Zhi Huang. 
                        The underelying pseudogene alignment database is generated by Travis S Johnson and Sihong Li using the Ohio Supercomputer Center (OSC). 
                        A second query tool which allows for direct sequence queries is also available through the OSC developed by Eric Franz and Travis S Johnson. 
                        
                        The online Shiny app is still under improvement. The app has been tested on OSX El Captitan 10.11.6 and R version 3.4.2. 
                        The website version has been tested in Chrome, Firefox, and Safari.
                        
                        Travis Johnson, Sihong Li, Eric Franz, Zhi Huang, Shuyu Dan Li, Moray J Campbell, Kun Huang, Yan Zhang. PseudoFuN: a resource to derive functional potentials in pseudogenes. Submitted.
                        
                        ## What does PseudoFuN do?
                        Pseudogene-gene families are sets of homologous genes and pseudogenes
                        based on sequence similarity. This app displays a set of the more than 26000
                        pseudogene-gene families that were generated using our pipeline. This database allows a
                        user to query gene symbols, Ensembl IDs and Entrez IDs in our database. The
                        pseudogene-gene family networks containing these queries are then displayed interactively
                        using networkD3. The genes in these networks can also be used to run Gene Ontology analysis
                        with topGO to assess the functional potentials of these pseudogenes. 
                        A more comprehensive query tool is also developed in association with the Ohio Supercomputer Center which
                        allows users to not only query gene names but also query sequences and then include the new query into
                        the assigned PGG network. The query will be performed on GPU-based clusters. This application also allows users to interactively view the GO
                        information for genes in the network individually, run multiple queries simultaneously and
                        download the results to a csv file.  For access to the more expansive tool please contact
                        Travis S Johnson (travis.johnson@osumc.edu) or Dr. Yan Zhang (yan.zhang@osumc.edu).")
                 ))),
    shiny::uiOutput("__reactivedoc__")
  )
)# end of navbar page


