#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#
library(tidyverse)
tree_data <- readr::read_tsv("/Users/thijsjanzen/found_stats.txt")
available_stats <- colnames(tree_data)[2:55]

library(shiny)


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Comparison of tree statistics"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            selectInput(inputId = "x_axis",
                        label = "Statistic1",
                        selected = "colless",
                        choices = available_stats),
            selectInput(inputId = "y_axis",
                        label = "Statistic2",
                        selected = "sackin",
                        choices = available_stats),
            checkboxInput("use_size",
                          "Color according to tree size?",
                          value = FALSE)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
        x <- found_stats %>%
            select(input$x_axis)
        y <- found_stats %>%
            select(input$y_axis)
        s <- found_stats %>%
            select(number_of_lineages)

        to_plot <- cbind(x, y, s)

        colnames(to_plot) <- c("x", "y", "size")
        if (!input$use_size) to_plot$size <- 1
        p1 <- ggplot(to_plot, aes(x = x, y = y, col = size)) +
            geom_point() +
            stat_smooth() +
            xlab(input$x_axis) +
            ylab(input$y_axis) +
            theme_classic()
        if (!input$use_size) p1 <- p1 + theme(legend.position = "none")
        p1
    })
}

# Run the application
shinyApp(ui = ui, server = server)



x <- unlist(found_stats %>%
    select(input$x_axis))
y <- unlist(found_stats %>%
    select(input$y_axis))

x <- unlist(as.vector(local_stats[stat1]))
y <- unlist(as.vector(local_stats[stat2]))
z <- unlist(as.vector(local_stats["number_of_lineages"]))

a1 <- nlme::gls(y~z, correlation = ape::corBrownian(1, local_tree))
a2 <- nlme::gls(x~z, correlation = ape::corBrownian(1, local_tree))

found_cor <- cor(a1$residuals, a2$residuals)
plot(a1$residuals~a2$residuals)
cor(a1$residuals,a2$residuals)
plot(x~y)
