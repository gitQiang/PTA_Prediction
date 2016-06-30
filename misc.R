csv_clean <- function(filename,wfile=""){
        
        ## file format: date and other indexs
        oneTable <- read.csv(filename,strip.white=TRUE)
        
        ### repalce special chars
        oneTable[,1] <- gsub("年","-",oneTable[,1])
        oneTable[,1] <- gsub("月","-",oneTable[,1])
        oneTable[,1] <- gsub("日","",oneTable[,1])
        oneTable[,1] <- gsub("/","-",oneTable[,1])
        for(j in 2:ncol(oneTable)) oneTable[,j] <- gsub(",","",oneTable[,j])
        
        ### date order and format
        n.unit <- length(unlist(strsplit(oneTable[1,1],"-")))
        if(n.unit==1){ ## in year unit
                oneTable[,1] <- paste(oneTable[,1],"-1-1",sep="")   
        }else if(n.unit==2){
                oneTable[,1] <- paste(oneTable[,1],"1",sep="")          
        }
        date <- as.Date(oneTable[,1])
        oneTable[,1] <- paste(year(date),month(date),day(date),sep="-")
        oneTable <- oneTable[order(date), ]
        
        if(wfile!="") write.csv(oneTable,file=wfile,row.names=FALSE,quote = FALSE)
        
        oneTable
        
}

fill_NA <- function(tmpTable,len=4){
        
        tmpNew <- tmpTable[,1]
        subs <- rep(FALSE,ncol(tmpTable))
        subs[1] <- TRUE
        for(i in 2:ncol(tmpTable)){
                tmp1 <- tof(tmpTable[,i])
                tmp2 <- na.approx(tmp1,maxgap=len,na.rm=FALSE) 
                if(!all(is.na(tmp2))){
                        tmpNew <- cbind(tmpNew,tmp2)
                        subs[i] <- TRUE
                }
        }
        colnames(tmpNew) <- colnames(tmpTable)[subs]
        
        tmpNew
}

data_filling <- function(){
        
        path <- "D:/data/恒逸/恒逸共享/调研数据整理/data_v1"
        filenames <- list.files(path,".csv",full.names=TRUE)
        useYears <- 2002:2016
        useDates <- seq.Date(as.Date("2002-1-1"),as.Date("2016-12-31"),by="day")
        useMonths <- paste(year(useDates),month(useDates),sep="-")
        dataM <- c()
        factors <- c()
        
        for(i in 1:length(filenames)){
                tmpTable <- csv_clean(filenames[i])
                tmpDate <- as.Date(tmpTable[,1])
                tmpTable <- tmpTable[year(tmpDate) %in% useYears,  ]
                tmpTable <- fill_NA(tmpTable,len=4)
                
                if(ncol(tmpTable) >= 2){
                        tmpDate <- as.Date(tmpTable[,1])
                        tmpName <- colnames(tmpTable)
                        
                        for(j in 2:ncol(tmpTable)){
                                onex <- rep(NA,length(useDates))
                                onex[match(tmpDate,useDates)] <- tmpTable[,j]
                                onex <- fill_onex(useMonths,tof(onex))
                                onex <- na.approx(onex,maxgap=62,na.rm=FALSE)
                                dataM <- cbind(dataM,onex)
                                factors <- c(factors,colnames(tmpTable)[j])
                        }
                }
                
        }
        rownames(dataM) <- as.character(useDates)
        colnames(dataM) <- factors

        dataM
}

fill_onex <- function(useMonths,onex){
        uniMon <- unique(useMonths)
        for(i in 1:length(uniMon)){
                tmpsub <- which(useMonths==uniMon[i])
                tmpmon <- onex[tmpsub]
                if(!all(is.na(tmpmon))){
                        if(sum(is.na(tmpmon)) <= 15 ){ 
                                tmpmon  <- na.approx(tmpmon,maxgap=15,na.rm=FALSE)
                        }else{
                                tmpmon[is.na(tmpmon)]  <- mean(tmpmon[!is.na(tmpmon)])  
                        }
                        onex[tmpsub] <- tmpmon
                }
        }
        onex
}

groupPredict <- function(data,i){
        
        useDates <- as.Date(rownames(data))
        if(i==1){
                tmpdata <- data
        }else if(i==2){
                gs <- paste(year(useDates),week(useDates),sep="-")
                ws <- unique(gs)
        }else if(i==3){
                gs <- paste(year(useDates),month(useDates),sep="-")
        }else if(i==4){
                gs <- paste(year(useDates),quarter(useDates),sep="-")
        }
        
        if(i > 1){
                tmpdata <- sapply(1:ncol(data), function(i)  sapply(unique(gs), function(ix) mean(data[gs==ix,i])) )
                colnames(tmpdata) <- colnames(data)                
        }
        
        if(i==2) rownames(tmpdata) <- ws
        
        tmpdata
}

oneDimPredict <- function(data,targetIndex,fre,per,sflag){
        
        jobid <- 6
        jobid <- jobTrace(jobid,trace)
        
        
        sdv <- sapply(1:ncol(data), function(i) sd(data[!is.na(data[,i]),i])) ## delete constant columns
        data <- data[, sdv!=0]
        newdata <- timelag_data(data,targetIndex,fre=fre)$newdata ## time lags for all factors
        sub <- max(delete_NA(newdata[,targetIndex]))-round(0.1*nrow(newdata)) ## delete discontinuous factors
        newdata <- newdata[,!is.na(newdata[sub,])]
        ### add time variable
        tmp <- colnames(newdata)[-targetIndex]
        newdata <- cbind(newdata[,targetIndex],1:nrow(newdata),newdata[,-targetIndex])
        colnames(newdata) <- c("Target","Time",tmp)
        
        jobid <- jobTrace(jobid,trace) # 2: the regression model with arima errors

        newdata <- newdata[!is.na(newdata[,1]), ]
        vNA <- rowSums(is.na(newdata))
        vNA <- vNA[1:(which(vNA==min(vNA))[1]-1)]
        dNA <- sort(unique(vNA),decreasing = TRUE)
        
        i=1
        startT <- which(vNA==dNA[i])[1]
        tmpnewdata <- newdata[startT:nrow(newdata), !is.na(newdata[startT, ])]
        endT <- min(which(is.na(tmpnewdata),arr.ind=TRUE)[,1])
        endT <- ifelse(endT==Inf,nrow(tmpnewdata),endT)
        tmpnewdata <- tmpnewdata[1:(endT-1), !is.na(tmpnewdata[endT-1, ])]
                
        jobid <- jobTrace(jobid,trace)
        ### output
        perform <- list();k=1;
        perform[[k]] <- oneModel(tmpnewdata,targetIndex,per=per,fre=fre,sflag=sflag,plotP=TRUE)
        k <- k+1
        pseudoR <- pseudoPredict(tmpnewdata,per,targetIndex)
        perform[[k]] <- pseudoR
        
        perform
}

oneModel <- function(tmpnewdata,sub=1,per=per,fre=fre,sflag=sflag,plotP=TRUE){

        #### output 
        R2 <- 1:per ### R^2 values
        preds <- 1:per ### predict values
        residuals <- 1:per
        para <- c()
        
        n2 <- nrow(tmpnewdata) ### input tmpnewdata length
        labs0 <- rownames(tmpnewdata)      
        
        if(sflag==1) pdq <- c(1,1,2,0,0,0,1) #day
        if(sflag==2) pdq <- c(1,1,1,2,1,0,7) #week
        if(sflag==3) pdq <- c(0,1,1,0,1,1,11) #month
        if(sflag==4) pdq <- c(0,1,0,0,1,0,2) #year
        
        for(n1 in (n2-per):(n2-1)){
                if(sflag<4){
                        #pdq <- sarima_paraNew(tmpnewdata[1:n1,sub],fre=fre)
                        stsr <- arima(ts(tmpnewdata[1:n1,sub],frequency = pdq[7]),order=pdq[1:3],seasonal = list(order=pdq[4:6],period=pdq[7]))
                        tmpdata <- tmpnewdata[1:n1,]
                        tmpdata[,sub] <- stsr$residuals
                        xnew <-stepCV_hq(tmpdata,sub,cvf=1,dir="forward")
                        mlmr <- lm(paste(colnames(tmpdata)[sub]," ~ ",paste(xnew,sep="",collapse = "+"),sep=""), data=as.data.frame(tmpdata))
                        res1 <- mlmr$residuals
                }
                
                if(sflag==4){
                        tmpdata <- as.data.frame(tmpnewdata[1:n1,])
                        xnew <-stepCV_hq(tmpdata,sub,cvf=1,dir="forward")
                        mlmr <- lm( paste(colnames(tmpdata)[sub]," ~ ",paste(xnew,sep="",collapse = "+"),sep=""), data=tmpdata )
                        res <- ts(mlmr$residuals,frequency = fre)
                        pdq <- sarima_paraNew(res,fre=fre)
                        stsr <- arima(res,order=pdq[1:3],seasonal = list(order=pdq[4:6],period=pdq[7]))
                        res1 <- stsr$residuals
                }
                
                #### one step prediction
                preds[n1-(n2-per-1)] <- sum(c(1,tmpnewdata[n1+1,xnew]) * tof(as.vector(mlmr$coefficients),na.rm=FALSE,fill=TRUE,f="zero"))  +  predict(stsr, n.ahead=1)$pred[1]
                oneP <- pdq
                para <- cbind(para,oneP)
                R2[n1-(n2-per-1)] <- R_squared_hq(tmpnewdata[1:n1,sub],tmpnewdata[1:n1,sub]-res1)
                residuals[n1-(n2-per-1)] <- tmpnewdata[n1+1,sub] - preds[n1-(n2-per-1)]
                
        }
        
        if(plotP){
                main="";ylab="Price";xlab="";ord="topleft";
                #### plot predict result
                n1=n2-per
                ymax <- 1.1*max(c(tmpnewdata[(n1+1):n2,sub],as.vector(preds)))
                ymin <- 0.9*min(c(tmpnewdata[(n1+1):n2,sub],as.vector(preds)))
                labs=rownames(tmpnewdata)[(n1+1):n2]
                par(mai=c(1.2,1.2,1,1))
                plot(tmpnewdata[(n1+1):n2,sub],col=1,type="b",ylim=c(ymin,ymax),ylab=ylab,xlab=xlab,xaxt="n",lwd=3,main=main)
                lines(as.vector(preds),col=2,type="b",lwd=3)
                axis(1,at=1:(n2-n1),labels =FALSE)  
                
                pos <- 1:(n2-n1)-per/100
                if(length(pos) > 20){pos <- pos[seq(1,length(pos),length.out=20)];labs <- labs[seq(1,length(labs),length.out=20)];}
                text(pos, par("usr")[3]-0.11*(ymax-ymin), labels = labs, srt = 90, pos = 1, xpd = TRUE)
                legend(ord,legend=c("观察值","预测值"),col=1:2,lwd=2) 
        }
        
        jobid <- jobTrace(9,trace)
        
        list(R2=R2,preds=preds,residuals=residuals,para=para)
}

pseudoPredict <- function(tmpdata,per,targetIndex=1){
        R2 <- 1:per
        preds <- 1:per
        residuals <- 1:per
        n2 <- nrow(tmpdata)
        for(n1 in (n2-per):(n2-1)){
                R2[n1-(n2-per-1)] <- R_squared_hq(tmpdata[2:n1,targetIndex],data[1:(n1-1),targetIndex])
                preds[n1-(n2-per-1)] <- tmpdata[n1,targetIndex]
                residuals[n1-(n2-per-1)] <- tmpdata[n1+1,targetIndex] - preds[n1-(n2-per-1)]
        }
        
        list(R2=R2,preds=preds,residuals=residuals,para=-1)
}

sarima_paraNew <- function(x,fre=10,nlag=0){
        
        if(!is.ts(x)) x <- as.ts(x,frequency=fre)
        if(nlag==0) nlag <- min(length(x) - fre+1,200)
        ## reference: http://people.duke.edu/~rnau/arimrule.htm
        sd=1; ## one fixed value !!!!
        maxd <- 5
        
        tmp <- acf2_hq(x,max.lag=nlag)
        d=0;
        if(all(tmp[1:5,"ACF"]>0)){
                sdd <- sapply(1:maxd, function(i) sd(diff(x,i)))
                d <- which.min(sdd)
        }
        
        if(d>0){
                dx <- diff(x,d)
                tmp <- acf2_hq(dx,max.lag = min(nlag,length(dx)-2))
        }else{ dx <- x; }
        
        dashV <- -1/length(dx) + 1.96/sqrt(length(dx))
        p <- which(abs(tmp[,"PACF"])<dashV)[1] - 1
        p <- ifelse(p<0,0,p)
        
        q <- which(abs(tmp[,"ACF"])<dashV)[1] - 1
        q <- ifelse(q<0,0,q)
        
        peak0 <- which(tmp[,"ACF"] > dashV)
        dpeak0 <- diff(peak0,1)
        fres0 <- dpeak0[dpeak0 > fre/2 & dpeak0 < 1.5*fre]
        sma <- length(fres0)
        peak1 <- which(tmp[,"PACF"] > dashV)
        dpeak1 <- diff(peak1,1)
        fres1 <- dpeak1[dpeak1 > fre/2 & dpeak1 < 1.5*fre]
        sar <- length(fres1)
        fres <- unique(c(fres0,fres1))
        
        if(fre==1){sar=0;sma=0;fres=1;}
        
        minBIC <- 1e9
        pdq <- c(p,d,q,0,sd,0,fre)
        for(f in fres){
                tmp <- auto.arima(ts(x,frequency = f), d=d, D=sd, max.p=max(p,1), max.q=max(q,1), max.P=max(sar,1), max.Q=max(sma,1), max.d=d, max.D=sd, start.p=0, start.q=0, start.P=0, start.Q=0, seasonal=TRUE, ic="bic",stepwise = FALSE)
                onebic <- BIC(tmp)
                if(onebic < minBIC){ minBIC <- onebic; pdq <- tmp$arma[c(1,6,2,3,7,4,5)]; }
        }

        pdq
}

stepCV_hq <- function(data,sub=1,cvf=1,dir="backward"){
        ### which measure is used: 1:CV;2:AIC;3:AICc;4:BIC;5:AdjR2 
        
        Y <- colnames(data)[sub]
        X <- colnames(data)[-sub]
        CVs <- sapply(X, function(i) CV(lm( paste(Y,"~",i,sep=""),  data=as.data.frame(data) ))[cvf])
        if(cvf==5){CVs = 1-CVs;}
        
        
        if(dir=="backward"){
                fit0 <- lm(paste(Y,"~.",sep=""), data=as.data.frame(data))
                print(CV(fit0))
                CV0 <- ifelse(cvf<=4,CV(fit0)[cvf],1-CV(fit0)[cvf])
                tmp <- sort(CVs,decreasing = TRUE, index.return=TRUE)
                X <- X[tmp$ix]
                Xnew <- X;
                for(i in 1:length(X)){
                        tmpfit <-  lm(paste(Y," ~ ",paste(setdiff(Xnew,X[i]),sep="",collapse=" + "),sep=""), data=as.data.frame(data))
                        tmpCV <- ifelse(cvf<=4,CV(tmpfit)[cvf],1-CV(tmpfit)[cvf])
                        print(tmpCV)
                        print(CV0)
                        if(tmpCV < CV0){fit0 <- tmpfit; CV0 <- tmpCV; Xnew <- setdiff(Xnew,X[i]);}
                }
        }
        
        if(dir=="forward"){
                tmp <- sort(CVs, index.return=TRUE)
                X <- X[tmp$ix]
                Xnew <- X[1];
                fit0 <-  lm(paste(Y," ~ ",X[1],sep=""), data=as.data.frame(data))
                CV0 <- ifelse(cvf<=4,CV(fit0)[cvf],1 - CV(fit0)[cvf])
                for(i in 2:length(X)){
                        tmpfit <-  lm(paste(Y," ~ ",paste(c(Xnew,X[i]),sep="",collapse=" + "),sep=""), data=as.data.frame(data))
                        tmpCV <- ifelse(cvf<=4,CV(tmpfit)[cvf],1-CV(tmpfit)[cvf])
                        if(tmpCV < CV0){fit0 <- tmpfit; CV0 <- tmpCV; Xnew <- c(Xnew,X[i]);}
                }
        }
        
        Xnew
}

timelag_data <- function(data,targetIndex,per=20,fre=12,n.model=3){
        newdata <- data[,targetIndex]
        newNames <- c()
        cNames <- colnames(data)
        
        xcols <- (1:ncol(data))[-targetIndex]
        lags <- c()
        
        for(i in xcols){
                tmpsubs <- delete_NA(tof(data[,targetIndex]),tof(data[,i]))
                tmpdata <- data[tmpsubs,c(targetIndex,i)]
                
                tmpccf <- ccf(tmpdata[,2],tmpdata[,1],lag.max=min(30,nrow(data)),plot=FALSE,type="correlation")
                lag <- tmpccf$lag[which.max(abs(tmpccf$acf[tmpccf$lag <= fre]))]
                #lag <- max(-fre, tmpccf$lag[which.max(abs(tmpccf$acf[tmpccf$lag<0]))])
                if(length(lag)==0) lag=0
                
                bm=1
                Xone <- data[,i,drop=FALSE]
                None <- cNames[i]
                
                if(lag==0){
                        lag = -1
                }else if(lag > 0){
                        lag <- lag - fre
                        if(lag==0) lag <- -fre
                }
    
                newdata <- cbind(newdata,c(rep(NA,-lag), Xone[1:(nrow(data)+lag),1]));
                newNames <- c(newNames,None)
                lags <- c(lags,lag)
        }
        colnames(newdata) <- c("Target",newNames)
        rownames(newdata) <- rownames(data)
        lags <- cbind(xcols,lags)
        
        list(newdata=newdata,lags=lags)    
}

delete_NA <- function(oneV,twoV=NULL){
        if(length(twoV)==0){
                oneV <- tof(oneV)
                tmpsubs <- which(is.na(oneV))
                tmpsubs <- c(0,tmpsubs,length(oneV)+1)
                n.na <- length(tmpsubs)
                onesub <- which.max(tmpsubs[2:n.na]-tmpsubs[1:(n.na-1)])
                (tmpsubs[onesub]+1):(tmpsubs[onesub+1]-1)
        }else{
                tmpsubs <- which(is.na(oneV) | is.na(twoV))
                tmpsubs <- c(0,tmpsubs,length(oneV)+1)
                n.na <- length(tmpsubs)
                onesub <- which.max(tmpsubs[2:n.na]-tmpsubs[1:(n.na-1)])
                (tmpsubs[onesub]+1):(tmpsubs[onesub+1]-1)
        }
}

tof <- function(oneV,na.rm=FALSE,fill=FALSE,f="mean"){
        ## replace function of as.numeric
        options(warn=-1)
        oneV <- as.numeric(oneV)
        options(warn=0)
        
        if(na.rm) oneV <- oneV[!is.na(oneV)]
        if(fill){
                if(f=="mean") oneV[is.na(oneV)] <- mean(oneV[!is.na(oneV)])
                if(f=="min") oneV[is.na(oneV)] <- min(oneV[!is.na(oneV)])
                if(f=="median") oneV[is.na(oneV)] <- median(oneV[!is.na(oneV)])
                if(f=="zero") oneV[is.na(oneV)] <- 0
        }
        oneV
}

fraction_NA  <- function(tmp,pNA=0.5){
        mode(tmp) <- "numeric"
        ksubs <- rep(TRUE,ncol(tmp))
        ksubs <- sapply( 1:(ncol(tmp)), function(i) sum(is.na(tmp[,i])) <  pNA*nrow(tmp) )
        tmp <- tmp[,ksubs]
        tmp  
}

R_squared_hq <- function(y,y1){
        y_ <- mean(y)
        1-(sum((y-y1)^2)/sum((y-y_)^2))
}

acf2_hq <- function(series, max.lag = NULL){
        
        num = length(series)
        if (num > 49 & is.null(max.lag)) 
                max.lag = ceiling(10 + sqrt(num))
        if (num < 50 & is.null(max.lag)) 
                max.lag = floor(5 * log10(num))
        if (max.lag > (num - 1)) 
                stop("Number of lags exceeds number of observations")
        ACF = stats::acf(series, max.lag, plot = FALSE)$acf[-1]
        PACF = stats::pacf(series, max.lag, plot = FALSE)$acf
        
        ACF <- round(ACF, 2)
        PACF <- round(PACF, 2)
        return(cbind(ACF, PACF))
}

#=======job tracing and write log file
jobTrace <- function(job,trace=0){
        
        logFile = "PTArunning.log"
        jobNames <- c("Data Loading ... ...","Start Day Prediction ... ...","Start Week Prediction ... ...","Start Month Prediction ... ...", "Start Season Prediction ... ...", "Start Data integrating ... ...","Start Multiple factor selection ... ...","Start PTA-PPM module ... ...","End PTA-PPM module ... ...","","","")
        
        if(trace!=0) print(jobNames[job])
        printstr <- paste(Sys.time(), Sys.info()["user"],jobNames[job],sep="\t")
        cat(printstr, file=logFile, append=TRUE, sep = "\n")
        
        job <- job+1
        job
}
