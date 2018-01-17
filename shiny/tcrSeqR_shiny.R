library(shiny)

#load a tcr object called 'ds'

timepoints <- levels(ds$type)
metrics <- c("Clonality", "Number.of.Expanded.Clones",
	     "Richness","Total.Sequences","Morisita",
	     "Log2.Fold.Change.in.Clonality","Log2.Fold.Change.in.Richness")
x_vals <- c('response', 'age','os.m','arm','nonskin.irae','sex')


ui <- fluidPage(
		titlePanel('tcrSeqR Shiny'),
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
		if(class(colData(ds)[,input$x_val])=='factor' | 
		   class(colData(ds)[,input$x_val])=='character' ){
			g <- iseqr_plot_factor(ds, 
					       metric=input$metric, 
					       by=input$x_val, 
					       type=input$type, 
					       labels=F)
			g + xlab(input$x_val)
		}else if(class(colData(ds)[,input$x_val])=='integer' | 
			 class(colData(ds)[,input$x_val])=='numeric'){
			g <- iseqr_plot_metrics(ds, 
					       metric=input$metric, 
					       by=input$x_val, 
					       type=input$type, 
					       labels=F)
			g + xlab(input$x_val)

		}
	})

}
shinyApp(ui=ui, server=server)
