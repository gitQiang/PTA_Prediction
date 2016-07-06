library(astsa)


plot_PTA_Price <- function(data){
        
        for(i in 1:4){
                
                tmpdata <- groupPredict(data,i)
                tmpdata1 <- sapply(1:ncol(tmpdata), function(i) na.approx(tof(tmpdata[,i]),maxgap=5,na.rm=FALSE))
                colnames(tmpdata1) <- colnames(tmpdata)
                rownames(tmpdata1) <- rownames(tmpdata)
                tmpdata1 <- fraction_NA(tmpdata1,pNA=0.5)
                
                x <- diff(tmpdata1[,1],1)
                plot_acfs(x)
                
                x <- log(tmpdata1[,1])
                plot_acfs(x)
                
                x <- diff(log(tmpdata1[,1]),1)
                plot_acfs(x)
        }
        
}

plot_acfs <- function(x){
        x <- as.numeric(x)
        x <- x[!is.na(x)]

        plot(x,type="l",main="Series:x",ylab="",xlab="")
        squared_x <- x^2
        abs_x <- abs(x)
        acf2(x)
        acf2(squared_x)
        acf2(abs_x)
}

