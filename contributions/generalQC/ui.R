menuList <- list(
  menuItem("General QC",
    tabName = "generalQC", startExpanded = FALSE,
    menuSubItem("UMI histogram", tabName = "umiHist"),
    menuSubItem("Sample histogram", tabName = "sampleHist"),
    menuSubItem("PC variance", tabName = "variancePC"),
    menuSubItem("TSNE plot", tabName = "tsnePlot"),
    menuSubItem("Umap", tabName = "umapPlot")
  )
)

tabList <- list(
  tabItem(
    "umiHist",
    tags$h3("Histogram of UMI counts"),
    fluidRow(column(
      10,
      offset = 1,
      plotOutput("plotUmiHist") %>% withSpinner()
    ))
  ),

  tabItem(
    "sampleHist",
    tags$h3("Histogram of cells per sample"),
    fluidRow(column(
      10,
      offset = 1,
      plotOutput("plotSampleHist") %>% withSpinner()
    ))
  ),

  tabItem(
    "variancePC",
    tags$h3("Variance of PCs"),
    fluidRow(column(
      10,
      offset = 1,
      plotOutput("variancePCA") %>% withSpinner()
    ))
  ),
  tsnePlotTab = tabItem(
    tabName = "tsnePlot",
    fluidRow(div(h3("TSNE Plot"), align = "center")),
    br(),
    fluidRow(
      column(
        3,
        numericInput("tsneDim", "Tsne dimensions", 3, min = 3, max = 3)
      ),
      column(
        3,
        numericInput("tsnePerplexity", "Tsne tsnePerplexity", 30, min = 1, max = 100)
      ),
      column(
        3,
        numericInput("tsneTheta", "Tsne tsneTheta", 0.5, min = 0.1, max = 10)
      ),
      column(
        3,
        numericInput("tsneSeed", "Tsne tsneSeed", 1, min = 1, max = 10000)
      )
    ),
    # fluidRow(column(
    #   12,
    #  )),
    # fluidRow(column(
    #   12,
    #
    # )),
    # fluidRow(column(
    #   12,
    #
    # )),
    fluidRow(
      column(3, selectInput(
        "dim3D_x",
        label = "X",
        choices = c("tsne1", "tsne2", "tsne3"),
        selected = "tsne1"
      )),
      column(
        3,
        selectInput(
          "dim3D_y",
          label = "Y",
          choices = c("tsne1", "tsne2", "tsne3"),
          selected = "tsne2"
        )
      ),
      column(
        3,
        selectInput(
          "dim3D_z",
          label = "Z",
          choices = c("tsne1", "tsne2", "tsne3"),
          selected = "tsne3"
        )
      ),
      column(
        3,
        selectInput(
          "col3D",
          label = "colored by",
          choices = c("tsne1"),
          selected = "tsne1"
        )
      )
    ),
    fluidRow(column(
      12,
      jqui_resizable(plotlyOutput("tsne_main"))
    )),

    fluidRow(column(
      10,
      offset = 1,
      tableSelectionUi("cellSelectionTSNEMod")
    ))
  ),
  umapTab <- tabItem(

    # unused parameters:
    # learning_rate = 1,
    # scale = FALSE,
    # init_sdev = NULL,
    # repulsion_strength = 1, a = NULL,
    # b = NULL, nn_method = NULL, n_trees = 50,
    # search_k = ifelse(n_refine_iters > 0, n_neighbors * n_trees, 2 *
    #                     n_neighbors * n_trees),
    # n_refine_iters = 0, approx_pow = FALSE,
    # y = NULL, target_n_neighbors = n_neighbors,
    # target_metric = "euclidean", target_weight = 0.5, pca = NULL,
    # pca_center = TRUE, ret_model = FALSE, ret_nn = FALSE,
    # n_sgd_threads = 0, grain_size = 1, verbose = getOption("verbose",
    # TRUE))
    tabName = "umapPlot",
    fluidRow(checkboxInput("activateUMAP", "activate Umap projection", FALSE)),
    fluidRow(
      column(
        3,
        selectInput(
          "um_randSeed",
          label = "random seed",
          choices = c(1:100), selected = "1"
        )
      ),
      column(
        3,
        selectInput(
          "um_n_neighbors",
          label = "N Neighbors",
          choices = c(2:100), selected = "15"
        )
      ),
      column(
        3,
        selectInput(
          "um_n_components",
          label = "N components",
          choices = c(2:20), selected = "2"
        )
      ),
      column(
        3,
        selectInput(
          "um_negative_sample_rate",
          label = "negative_sample_rate",
          choices = c(1:50), selected = "5"
        )
      )
    ),
    fluidRow(
      column(
        3,
        selectInput(
          "um_metric",
          label = "metric",
          choices = c("euclidean", "manhattan", "cosine", "hamming"),
          selected = "euclidean"
        )
      ),
      column(
        3,
        selectInput(
          "um_n_epochs",
          label = "n_epochs",
          choices = c(1:1000), selected = "200"
        )
      ),
      # selectInput(
      #   "um_alpha", label = "alpha",
      #   choices = seq(0.1,10,0.1), selected = "1.0"
      # ),
      column(
        3,
        selectInput(
          "um_init",
          label = "init",
          choices = c("spectral", "random"), selected = "spectral"
        )
      ),
      column(
        3,
        selectInput(
          "um_spread",
          label = "spread",
          choices = c(1:10), selected = "1"
        )
      )
    ),
    fluidRow(
      column(
        3,
        selectInput(
          "um_min_dist",
          label = "min_dist",
          choices = seq(0.05, 0.5, 0.01), selected = "0.01"
        )
      ),
      column(
        3,
        selectInput(
          "um_set_op_mix_ratio",
          label = "set_op_mix_ratio",
          choices = seq(0, 1, 0.1), selected = "1"
        )
      ),
      column(
        3,
        selectInput(
          "um_local_connectivity",
          label = "local_connectivity",
          choices = 1:20, selected = "1"
        )
      ),
      column(
        3,
        selectInput(
          "um_bandwidth",
          label = "bandwidth",
          choices = c(1:20), selected = "1"
        )
      )
    ),
    # fluidRow(
    #   column(
    #     3,
    #     # selectInput(
    #     #   "um_gamma", label = "gamma",
    #     #   choices = seq(0,10,0.2), selected = "1"
    #     # ),
    # 
    #     # selectInput(
    #     #   "um_umap1", label = "dim 1 to plot",
    #     #   choices = paste0("UMAP", 1:20), selected = "UMAP1"
    #     # ),
    #     # selectInput(
    #     #   "um_umap2", label = "dim 2 to plot",
    #     #   choices = paste0("UMAP", 1:20), selected = "UMAP2"
    #     # ),
    #   )
    # ),
    fluidRow(column(
      12,
      clusterUI("umap_main")
    ))
  )
)
