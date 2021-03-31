#' Shiny app User interface
#'
#' @return UI
#' @import shinycssloaders
#' @import DT
#' @rawNamespace import(plotly, except = c(slice,config))
#' @rawNamespace import(shinydashboard, except = c(valueBox,valueBoxOutput,renderValueBox))
#' @rawNamespace import(shiny, except = c(dataTableOutput,renderDataTable))
#' @export
#'
#' @examples
#' ui<-ui()
ui <- function() {
  UI <- shinydashboard::dashboardPage(
    title = 'Proteogenomic Isoforms Group Comparison',
    skin = "yellow",
    shinydashboard::dashboardHeader(
      title = shiny::span(
        tags$img(
          src = "CurieIsoAndSpe.png",
          height = 48,
          width = 48
        ),
        "Proteogenomic Isoforms Group Comparison"
      ),
      titleWidth = "25%"
    ),
    shinydashboard::dashboardSidebar(
      shiny::uiOutput("SampleGroupOrder"),
      shiny::uiOutput("Group"),
      shiny::checkboxInput(
        inputId = "Proteotypic",
        label = "Proteotypic peptides only",
        value = FALSE,
        width = NULL
      ),
      shiny::checkboxInput(
        inputId = "mRNA",
        label = "From mRNA",
        value = TRUE,
        width = NULL
      ),
      shiny::checkboxInput(
        inputId = "UNIPROT",
        label = "From UNIPROT",
        value = TRUE,
        width = NULL
      ),
      shiny::conditionalPanel(
        condition = "input.UNIPROT == true && input.mRNA == true",
        shiny::selectInput(
          "Match",
          "Correspondance",
          c("Perfect match",
            "Blast",
            "Both"),
          selected = "Both"
        )
      ),
      shiny::conditionalPanel(
        condition = "input.UNIPROT == true",
        shiny::selectInput(
          "Bank",
          "Protein Bank",
          list("Canonical",
               "Isoform",
               "TrEMBL",
               "All"),
          selected = "All"
        )
      ),
      shiny::conditionalPanel(
        condition = "input.Bank == 'Isoform' | input.Match== 'Blast'",
        shiny::checkboxInput(
          inputId = "Counterpart",
          label = "With counterpart",
          value = TRUE,
          width = NULL
        )
      ),
      shiny::actionButton('clear', 'Clear rows selection')
    ),
    shinydashboard::dashboardBody(
      htmlwidgets::getDependency('sparkline'),
      shinycssloaders::withSpinner(DT::dataTableOutput("Peptides"), color = "#FFCA05"),
      shiny::fixedRow(
        shiny::column(
          6,
          shinydashboard::box(
            title = "Present peptides informations",
            status = "warning",
            solidHeader = TRUE,
            shiny::htmlOutput("Description")
          )
        ),
        shiny::column(
          6,
          shinycssloaders::withSpinner(plotly::plotlyOutput("scatterPlotly"), color = "#FFCA05")
        )
      )
    )
  )
  return(UI)
}