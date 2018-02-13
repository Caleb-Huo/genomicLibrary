##' manhattan plot
##'
##' manhattan plot
##' @title manhattan plot
##' @param dataframe, same as qqman package, 4 columns including SNP, CHR, BP and P
##' @param title title
##' @param max.y max
##' @param suggestiveline 0
##' @return genomewideline default -log10(5e-8)
##' @return size.x.labels default 9
##' @return size.y.labels default 10
##' @return annotate F
##' @return SNPlist NULL
##' @author Caleb
##' @export
##' @examples
##' gwasResults <- data.frame(SNP=paste0("rs",1:40), CHR=1, BP = 1:40, P=runif(40))
##' manhattan(gwasResults)
##'
manhattan = function(dataframe, title=NULL, max.y="max", suggestiveline=0, genomewideline=-log10(5e-8), size.x.labels=9, size.y.labels=10, annotate=F, SNPlist=NULL) {
  
  if (annotate & is.null(SNPlist)) stop("You requested annotation but provided no SNPlist!")
  
  d=dataframe
  
  #limit to only chrs 1-23?
  d=d[d$CHR %in% 1:23, ]
  
  if ("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d) ) {
    
    d=na.omit(d)
    d=d[d$P>0 & d$P<=1, ]
    d$logp = -log10(d$P)
    
    d$pos=NA
    ticks=NULL
    lastbase=0
    
    #new 2010-05-10
    numchroms=length(unique(d$CHR))
    if (numchroms==1) {
      d$pos=d$BP
    } else {
      
      for (i in unique(d$CHR)) {
        if (i==1) {
          d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP
        }	else {
          lastbase=lastbase+tail(subset(d,CHR==i-1)$BP, 1)
          d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP+lastbase
        }
        ticks=c(ticks, d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
      }
      ticklim=c(min(d$pos),max(d$pos))
      
    }
    
    mycols=rep(c("gray10","gray60"),max(d$CHR))
    mycols=rainbow(max(d$CHR))
	set.seed(32608)
	mycols <- sample(mycols)
	
    if (max.y=="max") maxy=ceiling(max(d$logp)) else maxy=max.y
    if (maxy<8) maxy=8
    
    if (annotate) d.annotate=d[d$SNP %in% SNPlist, ]
    
    if (numchroms==1) {
      plot=qplot(pos,logp,data=d,ylab=expression(-log[10](italic(p))), xlab=paste("Chromosome",unique(d$CHR),"position"))
    }	else {
      plot=qplot(pos,logp,data=d, ylab=expression(-log[10](italic(p))) , colour=factor(CHR))
      plot=plot+scale_x_continuous(name="Chromosome", breaks=ticks, labels=(unique(d$CHR)))
      plot=plot+scale_y_continuous(limits=c(0,maxy), breaks=1:maxy, labels=1:maxy)
      plot=plot+scale_colour_manual(values=mycols)
    }
    
    if (annotate) 	plot=plot + # geom_point(data=d.annotate, colour=I("grey50")) + 
      geom_text(data=d.annotate, label=d.annotate$SNP, nudge_x = 0.25, nudge_y = 0.25, check_overlap = T, colour=I("grey50"))
    
    #plot=plot + theme() 
    #plot=plot + theme(title=title)
    plot=plot+ theme_bw() +
    theme(
      legend.position = "none",
      #panel.background=theme_blank(), 
      #panel.grid.minor=theme_blank(),
      axis.text.x=element_text(size=size.x.labels, colour="grey50"), 
      axis.text.y=element_text(size=size.y.labels, colour="grey50"), 
      axis.title.x=element_text(size=size.x.labels, colour="grey50"), 
      axis.title.y=element_text(size=size.y.labels, colour="grey50"), 
      #axis.ticks=theme_segment(colour=NA)
    ) +
	labs(title = title)
    
    if (suggestiveline) plot=plot+geom_hline(yintercept=suggestiveline,colour="blue", alpha=I(1/3))
    if (genomewideline) plot=plot+geom_hline(yintercept=genomewideline,colour="grey50")
    
    plot
    
  }	else {
    stop("Make sure your data frame contains columns CHR, BP, and P")
  }
}

