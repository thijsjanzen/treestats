library(InteractiveComplexHeatmap)

# htShinyExample(5.5)
require(readr)
require(magrittr)
require(tidyr)
require(dplyr)
library(shiny)
library(circlize)
require(ggpubr)

library(treestats)

data(emp_data)


res.cor <- emp_data$cor

plot_data <- emp_data$data

col_fun <- colorRamp2(seq(-1, 1, length.out = 100), ggpubr::get_palette(palette = "RdBu", k = 100))


ht = Heatmap(res.cor, name = "Correlation",
             col = col_fun,
             rect_gp = gpar(type = "none"),
             cell_fun = function(j, i, x, y, w, h, fill) {
               grid.rect(x, y, w, h, gp = gpar(fill = "transparent", col = "grey"))
               grid.circle(x = x, y = y, r = abs(res.cor[i, j])/2 * min(unit.c(w, h)),
                           gp = gpar(fill = col_fun(res.cor[i, j]), col = col_fun(res.cor[i, j])))
             },
             show_row_dend = FALSE, show_column_dend = FALSE)

ht <- draw(ht)

ui = shiny::fluidPage(
  InteractiveComplexHeatmapOutput(title1 = "Heatmap on 209 Empirical trees",
                                  response = "click",
                                  output_ui = shiny::plotOutput("scatterplot", width = 400, height = 400),
                                  width1 = 800,
                                  height1 = 800)
)

click_action = function(df, output) {
  output$scatterplot = shiny::renderPlot({
    if (is.null(df)) {
      grid.text("You should click on heatmap cells.")
    } else {
      nm = colnames(plot_data)
      i1 = 1 + df$column_index
      i2 = 1 + df$row_index

      x = as.vector(plot_data[, nm[i1]])[[1]]
      y = as.vector(plot_data[, nm[i2]])[[1]]

      c1 <- which(rownames(res.cor) == nm[i1])
      c2 <- which(rownames(res.cor) == nm[i2])

      found_cor <- res.cor[c1, c2]
      focal_color <- col_fun(found_cor)

      plot(x, y, xlab = nm[i1], ylab = nm[i2],
           main = paste0("Correlation = ", sprintf('%.3f', found_cor)),
           pch = 16, col = focal_color, cex = 0.5)
      abline(lm(y~x), col = focal_color, lwd = 2)
    }
  })
}


server = function(input, output, session) {
  makeInteractiveComplexHeatmap(input, output, session, ht,
                                click_action = click_action)
}

shiny::shinyApp(ui, server)
