cor.tests <- function(xdatas, ydata1s, ydata2s, ydata3s, trait, metric){
	
	out.cors <- replicate(4, matrix(NA, ncol=5, nrow=3), simplify=FALSE)
	for(i in 1:length(xdatas)){
		
		xdata1 <- xdatas[[i]]
		ydata1 <- ydata1s[[i]]
		ydata2 <- ydata2s[[i]]
		ydata3 <- ydata3s[[i]]
	
		for(j in 1:length(trait)){
			
			out.cors[[i]][j,1]  <- cor.test(xdata1[,metric][xdata1$link.rule == trait[j]],
				ydata1[,metric][ydata1$link.rule == trait[j]], method="spearman")$estimate
			
			out.cors[[i]][j,2]  <- cor.test(xdata1[,metric][xdata1$link.rule == trait[j]],
				ydata2[,metric][ydata2$link.rule == trait[j]], method="spearman")$estimate	
			
			out.cors[[i]][j,3]  <- cor.test(xdata1[,metric][xdata1$link.rule == trait[j]],
				ydata3[,metric][ydata3$link.rule == trait[j]], method="spearman")$estimate	
				
			out.cors[[i]][j,4]  <- cor.test(ydata1[,metric][ydata1$link.rule == trait[j]],
				ydata2[,metric][ydata2$link.rule == trait[j]], method="spearman")$estimate
				
			out.cors[[i]][j,5]  <- cor.test(ydata1[,metric][ydata1$link.rule == trait[j]],
				ydata3[,metric][ydata3$link.rule == trait[j]], method="spearman")$estimate		
				

			colnames(out.cors[[i]]) <- c("fun.eco.cor", 
										 "fun.samp1000.cor", 
										 "fun.samp100.cor", 
										 "eco.samp1000.cor", 
										 "eco.samp100.cor")
		}
	}
	
	names(out.cors) <- names(xdatas)
	return(out.cors)
}