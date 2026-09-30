options(stringsAsFactors = FALSE)
require(survival)
require(survminer)
require(RColorBrewer)
require(gridExtra)
library(grid)


#' @description Plot Kaplan-Meier survival curves and perform a log-rank test.
#' @param clinical.data A data frame with classification labels. It must contain at least four columns:
#'   {sample_ID}, {event} (numeric or logical event indicator), {time} (numeric follow-up time), and {sample.label} (a factor with group labels).
#' @param upper.time Numeric. Upper limit of follow-up time, in the same unit as {clinical.data}.
#'   Default is NULL. If used, samples with time greater than this value are removed.
#' @param xscale Character. Allowed values are "d_m", "d_y", "m_d", "m_y", "y_d", and "y_m", where d = day,
#'   m = month, and y = year. For example, xscale = "d_m" converts the x-axis unit from days to months. If NULL, no conversion is performed.
#' @param unit.xlabel Character. Time unit shown on the x-axis. One of "year", "month", "week", or "day". It should match the time unit in the final plot.
#' @param surv.median.line Character. Whether to draw horizontal/vertical lines at median survival.
#'   One of "none", "hv", "h", or "v". v: vertical, h: horizontal.
#' @param risk.table Logical. Whether to draw a number-at-risk table. Default is TRUE.
#' @param pval Logical. Whether to show the log-rank p-value. Default is TRUE.
#' @param conf Logical. Whether to show confidence intervals. Default is FALSE.
#' @param main Character. Main title of the plot. Default is NULL.
#' @param colors Character vector. Colors used for groups. If NULL, the Set1 palette is used.
#' @param survival.event Character. Type of survival event, such as "Overall Survival"or "Progress Free Survival".
#'
#' @return A ggsurvplot object.
#' @export


plot.surv <- function(clinical.data,
                      upper.time=NULL,
                      xscale = c( "d_m", "d_y", "m_d", "m_y", "y_d", "y_m"),
                      unit.xlabel = c("year", "month", "week", "day"),
                      surv.median.line = c("none", "hv", "h", "v"),
                      risk.table = c(TRUE, FALSE),
                      pval = c(TRUE, FALSE),
                      conf = c(FALSE, TRUE),
                      main=NULL,
                      colors=NULL,
                      survival.event=c("Overall Survival", "Progress Free Survival")
)
{
    survival.event <- survival.event[1]
    unit.xlabel <- unit.xlabel[1]

    if(!is.null(upper.time))
    {
        clinical.data <- clinical.data[clinical.data$time <= upper.time,]
    }
    
    xSL <- data.frame(xScale=c(1,7,30,365.25),xLab=c("Days","Weeks","Months","Years"), stringsAsFactors=FALSE)
    switch(unit.xlabel, year={xScale <- 365.25}, month={xScale <- 30}, week={xScale <- 7}, day={xScale <- 1})
    xLab <- xSL[which(xSL[,1]==xScale),2]
    
    t.name <- levels(clinical.data$sample.label)
    colors <- RColorBrewer::brewer.pal(n = 9, name = "Set1")
    t.col <- colors
    
    km.curves <- survfit(Surv(time, event)~sample.label, data=clinical.data)
    legend.content <- substr(names(km.curves$strata),start = 14,stop = 1000)
    
    ggsurv <- ggsurvplot(
                         km.curves,               
                         data = clinical.data,             
                         palette = t.col,
                         
                         risk.table = risk.table[1],       
                         pval = pval[1],          
                         surv.median.line = surv.median.line[1],  
                         title = main,     
                         font.main = 15,                  
                         xlab = paste("Time","in",xLab,sep = " "),  
                         ylab = survival.event,  
                         conf.int = conf,
						 
                         legend.title = "", 
                         legend.labs = legend.content, 
                         legend = c(0.8,0.9), 
                         font.legend = 9,    
                         
                         tables.theme = theme_cleantable(),
                         risk.table.title = "No. at risk:",
                         risk.table.y.text.col = T, 
                         risk.table.y.text = FALSE, 
                         tables.height = 0.15,    
                         risk.table.fontsize = 3  
                        );
    ggsurv$plot <- ggsurv$plot + theme(plot.title = element_text(hjust = 0.5), plot.margin = unit(c(5.5, 5.5, 5.5, 50), "points"))
    ggsurv$table <- ggsurv$table + theme(plot.title = element_text(hjust = -0.04), plot.margin = unit(c(5.5, 5.5, 5.5, 50), "points"))

    print(ggsurv)
}
