library(shiny)
library(shinydashboard)


server <- function(input, output, session)
{
  output$instructions <-
    renderImage({

      my_file <- system.file("apps/diversity/app_diagram.png", package = "microbiome.onglab")

      width  <- session$clientData$output_instructions_width

      list(
        src = my_file,
        contentType = "image/png",
        width = width,
        alt = "diversity app instructions"
      )


    }, deleteFile = FALSE)



}
