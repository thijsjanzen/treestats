#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#
library(ape)
library(ggplot2)
library(shiny)
library(nlme)

tree_data <- read.table("https://raw.githubusercontent.com/thijsjanzen/treestats/shiny2/data/emp_stats.txt",
                        header = TRUE)
available_stats <- colnames(tree_data)[2:55]

phy_tree <- ape::read.tree("https://raw.githubusercontent.com/thijsjanzen/treestats/shiny2/data/phy_tree.txt")

sim_data <- read.table("https://raw.githubusercontent.com/thijsjanzen/treestats/shiny2/data/sim_data.txt",
                       header = TRUE)



sidebarPanel2 <- function(..., out = NULL, width = 4)
{
  div(class = paste0("col-sm-", width),
      tags$form(class = "well", ...),
      out
  )
}



# Define UI for application that draws a histogram
ui <- fluidPage(

  # Application title
  titlePanel("Comparison of tree statistics"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel2(
                                selectInput(inputId = "x_axis",
                                           label = "Statistic1",
                                           selected = "colless",
                                           choices = sort(available_stats)),
                               selectInput(inputId = "y_axis",
                                           label = "Statistic2",
                                           selected = "sackin",
                                           choices = sort(available_stats)),
                               selectInput(inputId = "coloring",
                                           label = "Coloring",
                                           choices = c("None", "Size", "Taxonomic group")),
      out = HTML('This Shiny app provides a means to explore the found correlations in <a href =https://doi.org/10.1101/2024.01.24.576848>Janzen 2024</a>. Empirical correlations are based on 215 Empirical trees from <a href = https://doi-org.proxy-ub.rug.nl/10.1111/ele.13382>Condamine et al. 2019</a>. Simulated data shown is a random subset of 500 trees per diversification model, see the original paper for results using a much larger dataset. <br><br>
              There are three types of plots available: <br>
                1) Raw correlations, without correction. <br>
                2) Residual correlations, corrected for tree size and phylogenetic relatedness. <br>
                3) Correlations on simulated data, trees were simulated with fixed size using four different diversification models.')
    ),

    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(type = "tabs", id = "tabs1",
        tabPanel("Empirical data: Raw values", value = 1, plotOutput("distPlot")),
        tabPanel("Empirical data: Residuals",  value = 2, plotOutput("corPlot")),
        tabPanel("Simulated data", value = 3, plotOutput("simPlot"))
      ),

    )

  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  output$distPlot <- renderPlot({
    x <- tree_data[which(colnames(tree_data) == input$x_axis)]
    y <- tree_data[which(colnames(tree_data) == input$y_axis)]
    s <- tree_data$number_of_lineages
    tt <- tree_data$Taxa

    to_plot <- cbind(x, y, s, tt)

    colnames(to_plot) <- c("x", "y", "Tree Size", "Taxa")
    local_cor <- paste0("Pearson correlation = ", round(cor(x, y), 2))

    if (input$coloring == "Size") {
      p1 <- ggplot(to_plot, aes(x = x, y = y, col = (`Tree Size`))) +
        geom_point(size = 2) +
        scale_color_viridis_c(option = "C") +
        stat_smooth(method = "lm", col = "#416894", fill = "#416894") +
        xlab(input$x_axis) +
        ylab(input$y_axis) +
        theme_classic() +
        labs(col = "Tree Size") +
        ggtitle(local_cor) +
        theme(legend.text = element_text(size=14),
              legend.title = element_text(size=16),
              axis.title.x = element_text(size=16),
              axis.title.y = element_text(size=16),
              plot.title = element_text(size = 16, face = "bold"))
    }
    if (input$coloring == "None") {
      p1 <- ggplot(to_plot, aes(x = x, y = y)) +
        geom_point(col = "#416894", size = 2) +
        stat_smooth(method = "lm", col = "#416894", fill = "#416894") +
        xlab(input$x_axis) +
        ylab(input$y_axis) +
        theme_classic() +
        ggtitle(local_cor) +
        theme(legend.position = "none") +
        theme(legend.text = element_text(size=14),
              legend.title = element_text(size=16),
              axis.title.x = element_text(size=16),
              axis.title.y = element_text(size=16),
              plot.title = element_text(size = 16, face = "bold"))
    }
    if (input$coloring == "Taxonomic group") {
      p1 <- ggplot(to_plot, aes(x = x, y = y, col = Taxa)) +
        geom_point(size = 2) +
        stat_smooth(method = "lm", col = "#416894", fill = "#416894") +
        scale_color_brewer(type = "div", palette = 3) +
        xlab(input$x_axis) +
        ylab(input$y_axis) +
        theme_classic() +
        ggtitle(local_cor) +
        labs(col = "Taxonomic\nGroup") +
        theme(legend.text = element_text(size=14),
              legend.title = element_text(size=16),
              axis.title.x = element_text(size=16),
              axis.title.y = element_text(size=16),
              plot.title = element_text(size = 16, face = "bold")) +
        guides(colour = guide_legend(override.aes = list(size=3)))
    }

    p1
  })

  output$corPlot <- renderPlot({
    x <- unlist(tree_data[which(colnames(tree_data) == input$x_axis)])
    y <- unlist(tree_data[which(colnames(tree_data) == input$y_axis)])
    z <- tree_data$number_of_lineages
    tt <- tree_data$Taxa

    sp <- tree_data$Family
    bm <- ape::corBrownian(value = 1, phy = phy_tree, form =~ sp)

    a1 <- nlme::gls(x~z, correlation = bm)
    a2 <- nlme::gls(y~z, correlation = bm)

    xvals <- a1$residuals
    yvals <- a2$residuals
    to_plot <- cbind(xvals, yvals, z, tt)
    colnames(to_plot) <- c("x", "y", "Size", "Taxa")
    to_plot <- as.data.frame(to_plot)
    to_plot$x <- as.numeric(to_plot$x)
    to_plot$y <- as.numeric(to_plot$y)
    to_plot$Size <- as.numeric(to_plot$Size)

    local_cor <- paste0("Pearson correlation = ", round(cor(xvals, yvals), 2))

    if (input$coloring == "Taxonomic group") {
      p1 <- ggplot(to_plot, aes(x = x, y = y, col = as.factor(Taxa))) +
        geom_point(size = 2) +
        scale_color_brewer(type = "div", palette = 3) +
        stat_smooth(method = "lm", col = "#416894", fill = "#416894") +
        xlab(paste("residual:", input$x_axis)) +
        ylab(paste("residual:", input$y_axis)) +
        theme_classic() +
        ggtitle(local_cor) +
        labs(col = "Taxonomic\nGroup") +
        theme(legend.text = element_text(size=14),
              legend.title = element_text(size=16),
              axis.title.x = element_text(size=16),
              axis.title.y = element_text(size=16),
              plot.title = element_text(size = 16, face = "bold")) +
        guides(colour = guide_legend(override.aes = list(size=3)))
    }
    if (input$coloring == "None") {
      p1 <- ggplot(to_plot, aes(x = x, y = y)) +
        geom_point(col = "#416894", size = 2) +
        stat_smooth(method = "lm", col = "#416894", fill = "#416894") +
        xlab(paste("residual:", input$x_axis)) +
        ylab(paste("residual:", input$y_axis)) +
        theme_classic() +
        ggtitle(local_cor) +
        theme(legend.text = element_text(size=14),
              legend.title = element_text(size=16),
              axis.title.x = element_text(size=16),
              axis.title.y = element_text(size=16),
              plot.title = element_text(size = 16, face = "bold"))
    }
    if (input$coloring == "Size") {
      p1 <- ggplot(to_plot, aes(x = x, y = y, col = Size)) +
        geom_point(size = 2) +
        scale_color_viridis_c(option = "C") +
        stat_smooth(method = "lm", col = "#416894", fill = "#416894") +
        xlab(paste("residual:", input$x_axis)) +
        ylab(paste("residual:", input$y_axis)) +
        theme_classic() +
        ggtitle(local_cor) +
        theme(legend.text = element_text(size=14),
              legend.title = element_text(size=16),
              axis.title.x = element_text(size=16),
              axis.title.y = element_text(size=16),
              plot.title = element_text(size = 16, face = "bold"))
    }
    p1
  })

  output$simPlot <- renderPlot({

    x <- unlist(sim_data[which(colnames(sim_data) == input$x_axis)])
    y <- unlist(sim_data[which(colnames(sim_data) == input$y_axis)])
    used_model <- sim_data$model

    to_plot <- cbind(x, y, used_model)
    colnames(to_plot) <- c("x", "y", "model")
    to_plot <- as.data.frame(to_plot)
    to_plot$x <- as.numeric(to_plot$x)
    to_plot$y <- as.numeric(to_plot$y)
    local_cor <- paste0("Pearson correlation = ", round(cor(x, y), 2))
    ggplot(to_plot, aes(x = x, y = y, col = as.factor(model))) +
      geom_point(size = 2) +
      stat_smooth(method = "lm", col = "#416894", fill = "#416894") +
      scale_color_brewer(type = "div", palette = 3) +
      xlab(input$x_axis) +
      ylab(input$y_axis) +
      labs(col = "Diversification\nModel") +
      theme_classic() +
      ggtitle(local_cor) +
      theme(legend.text = element_text(size=14),
            legend.title = element_text(size=16),
            axis.title.x = element_text(size=16),
            axis.title.y = element_text(size=16),
            plot.title = element_text(size = 16, face = "bold")) +
      guides(colour = guide_legend(override.aes = list(size=3)))
  })
}

# Run the application
shinyApp(ui = ui, server = server)

