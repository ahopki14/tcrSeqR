library(shiny)

timepoints <- levels(plot_ds$type)
metrics <- names(plot_ds)[11:16]
x_vals <- names(plot_ds)[6:8]


ui <- fluidPage(
		titlePanel('immunoSeqR Shiny'),
		sidebarLayout(
		  sidebarPanel(
		  	selectInput(inputId='metric', label='Metric', choices=metrics),
		  	selectInput(inputId='type', label='Time Point', choices=timepoints),
		  	selectInput(inputId='x_val', label='Split By', choices=x_vals)
		  ),
		  mainPanel(
	       	    plotOutput(outputId="p")
		  )
		)
)


server <- function(input, output){
	output$p <- renderPlot({
		g <- iseqr_plot_factor(plot_ds, 
				       metric=input$metric, 
				       x_val=input$x_val, 
				       type=input$type, 
				       labels=F)
		g + xlab(input$x_val)
		#hist(rnorm(input$t))
	})

}
shinyApp(ui=ui, server=server)
