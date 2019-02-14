##' manhattan plot
##'
##' manhattan plot. Input a dataframe including CHR, BP, P, SNP
##' @title manhattan plot
##' @param dataframe, same as qqman package, 4 columns including SNP, CHR, BP and P
##' @param title title
##' @param maxy maximum of y axis
##' @param suggestiveline 0
##' @param genomewideline default -log10(5e-8)
##' @param axisSize default 20
##' @param labelSize default 5
##' @param annotate F
##' @param SNPlogic which SNP need to be annotated in the plot
##' @author Caleb
##' @export
##' @examples
##' gwasResults1 <- data.frame(SNP=paste0("rs1",1:40), CHR=1, BP = 1:400, P=runif(400))
##' gwasResults2 <- data.frame(SNP=paste0("rs2",1:40), CHR=2, BP = 1:600, P=runif(600))
##' gwasResults <- rbind(gwasResults1, gwasResults2)
##' manhattan(gwasResults)
##' manhattan(gwasResults, annotate=T, SNPlogic=c(100,200,300,310))
##'
manhattan <- function(dataframe, maxy=NULL, oneColor = "red", suggestiveline=0, genomewideline=-log10(5e-8), axisSize = 20, labelSize = 5, annotate=F, SNPlogic=NULL){
    
	if (annotate & is.null(SNPlogic)) stop("You requested annotation but provided no SNPlist!")
	
	d=dataframe
	
	stopifnot("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d) )
	
    d=d[d$CHR %in% 1:23, ]
    d$logp = -log10(d$P)
	
    d$pos=NA
    ticks=NULL
    lastbase=0
	
	sortedCHR <- sort(unique(d$CHR))
	numchroms=length(sortedCHR)
	
    for (i in sortedCHR) {
      if (i==sortedCHR[1]) {
        d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP
      }	else {
        lastbase=lastbase+tail(subset(d,CHR==i-1)$BP, 1)
        d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP+lastbase
      }
      ticks=c(ticks, d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
    }
    ticklim=c(min(d$pos),max(d$pos))
 
 	ncolor <- 6
	mycols0 <- palette()[1:ncolor]
	if(numchroms == 1){
		mycols0[1] <- oneColor
	} else {
		mycols0[1] <- "grey"		
	}
    d$color <- with(d, mycols0[CHR%%ncolor + 1])
	
	
	if(is.null(maxy)){
		maxy <- max(d$logp)
	}
	
	
	## visualize the manhattan plot
	plot <- ggplot(data=d,aes(x=pos,y=logp,col=color)) + 
	geom_point() +
	scale_x_continuous(name="Chromosome", breaks=ticks, labels=(unique(d$CHR))) +
	scale_y_continuous(limits=c(0,maxy), breaks=0:maxy, labels=0:maxy) + 
	scale_colour_manual(values=mycols0) + 
	coord_cartesian(xlim = c(min(d$pos), max(d$pos)), ylim = c(0, maxy*1.05), expand = FALSE) + 
	labs(x = "Chromosome", y = expression(-log[10](italic(p)))) + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
	theme(text = element_text(size=axisSize), legend.position="none") 
    if (suggestiveline) plot=plot+geom_hline(yintercept=suggestiveline,colour="black", linetype="dashed")
    if (genomewideline) plot=plot+geom_hline(yintercept=genomewideline,colour="black")
		
	if (annotate){
		d.annotate=d[SNPlogic, ]
		plot <- plot + geom_text_repel(data=d.annotate, aes(x=pos,y=logp,label=SNP), size=labelSize,colour=I("black"))
	} 
	
	plot
	

}


