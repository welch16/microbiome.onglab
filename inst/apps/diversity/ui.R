

## ui.R ##
library(shiny)
library(shinydashboard)

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Instructions", tabName = "instructions"),
    menuItem("Diversity", tabName = "diversity"),
    menuItem("PCoA", tabName = "pcoa"),
    menuItem("Relative Abundance", tabName = "relative"),
    fileInput(inputId = "stats_file", label = "A. Sample info."),
    fileInput(inputId = "pcoa_file", label = "B. PCoA obj."),
    fileInput(inputId = "taxa_file", label = "C. Variable info.")
    # menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
    # menuItem("Widgets", icon = icon("th"), tabName = "widgets",
    #          badgeLabel = "new", badgeColor = "green")
  )
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "instructions",
            h1("Instructions"),
            imageOutput("instructions"),
            h3("The idea of this app is to compare three data from three different sources:"),
            h3("  A. Sample information (i.e. a table with sample as first column)"),
            h3("  B. Principal Coordinate Analysis info. (the samples must be the same as in A.)"),
            h3("  C. Variable based information, e.g. the taxa labels")),
    tabItem(tabName = "diversity",
            h2("Diversity and Richness")),
    tabItem(tabName = "pcoa",
            h2("Principal Coordinates Analysis")),
    tabItem(tabName = "relative",
            h2("Relative abundance"))
    # tabItem(tabName = "dashboard",
    #         h2("Dashboard tab content")
    # ),
    #
    # tabItem(tabName = "widgets",
    #         h2("Widgets tab content")
    )
)

# Put them together into a dashboardPage
dashboardPage(
  dashboardHeader(title = "Diversity analysis"),
  sidebar,
  body
)
