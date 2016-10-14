#!/usr/bin/Rscript

.libPaths("/home/julius/R/x86_64-pc-linux-gnu-library/3.2")
library(ggplot2, quietly = T, warn.conflicts = F)
library(cowplot, quietly = T, warn.conflicts = F)
library(shiny, quietly = T, warn.conflicts = F)

df = read.table("/tmp/output_baf.csv", h=T)
calls = read.table("/tmp/output_baf.csv_calls", h=F)
genes = read.table("/tmp/output_baf.csv_genes", h=F)
genes$ln = genes$V3 - genes$V2
genes = genes[genes$ln>0,]
genes = genes[order(genes$ln, decreasing=T),]
genes = genes[!duplicated(genes$V4), ]
genes = genes[order(genes$V2),]
genes$ys = seq_along(genes$V1) %% 10

minlrr = min(df$lrr)
maxlrr = max(df$lrr)
nudgey = (maxlrr - minlrr)*0.05
defxlims = c(min(df$pos), max(df$pos))

ui = fluidRow(
        h4("Brush and double-click to zoom, double-click outside the brush to reset"),
        plotOutput("mainplot",
                   dblclick = "plot_dblclick",
                   brush = brushOpts(id = "plot_brush", direction="x", resetOnNew=T)
        ),
        actionButton(
                inputId = 'close', type = "button",
                onclick = "setTimeout(function(){window.close();},500);", "Stop app" ## dar neveikia
        )
)

server = function(input, output, session) {
        xlims = reactiveValues(x=defxlims)
        brush = NULL
        
        observe({ if(input$close>0) stopApp() })
        
        observeEvent(input$plot_dblclick, {
                brush = input$plot_brush

                if (!is.null(brush)) {
                        brushcorr = c(brush$xmin, brush$xmax)*(1.02)-0.02
                        xlims$x = brushcorr*(xlims$x[2]-xlims$x[1]) + xlims$x[1]
                } else {
                        xlims$x = defxlims
                }
        })
        
        output$mainplot = renderPlot({
                nudgex = (xlims$x[2] - xlims$x[1])*0.01
                pb = ggplot(df) + geom_point(aes(x=pos, y=baf), col="blue") +
                        geom_rect(data=calls, aes(xmin=V2, xmax=V3, ymin=0, ymax=1), alpha=0.1, col="gray60") + 
                        coord_cartesian(xlim=xlims$x, expand=FALSE) +
                        ylab("BAF") + xlab("position") + theme_gray()
                pb = ggdraw(switch_axis_position(pb, axis="x"))
                pl = ggplot(df) + geom_point(aes(x=pos, y=lrr), col="red") +
                        geom_rect(data=calls, aes(xmin=V2, xmax=V3, ymin=minlrr, ymax=maxlrr), alpha=0.1, col="gray60") + 
                        geom_text(data=calls, aes(x=V2, y=minlrr, label=V6), nudge_x=-nudgex, nudge_y=nudgey) +
                        coord_cartesian(xlim=xlims$x, expand=FALSE) +
                        ylab("LRR") + xlab(NULL) + theme_gray() +
                        theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
                pg = ggplot(genes) + geom_rect(aes(xmin=V2, xmax=V3, ymin=-1, ymax=max(ys)+1), fill="lightblue") +
                        geom_rect(data=calls, aes(xmin=V2, xmax=V3, ymin=-1, ymax=max(genes$ys)+1), alpha=0.1, col="gray60") + 
                        geom_text(aes(x=(V2+V3)/2, y=ys, label=V4)) +
                        geom_segment(aes(x=(V2+V3)/2, xend=(V2+V3)/2, y=ys-0.8, yend=ys-1.5)) +
                        coord_cartesian(xlim=xlims$x, expand=FALSE) +
                        xlab(NULL) + ylab("GENES") + theme_gray() +
                        theme(axis.text.y=element_text(color="white"), axis.ticks.y=element_line(color="white"),
                              panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank())
                
                plot_grid(pb, pl, pg, ncol = 1)        
        })
}

app = shinyApp(ui=ui, server=server)
runApp(app, host="192.168.0.133", port=3535)
