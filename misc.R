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
        
        ## delete constant columns
        sdv <- sapply(1:ncol(data), function(i) sd(data[!is.na(data[,i]),i]))
        data <- data[, sdv!=0]
        
        ## time lags for all factors
        newdata <- timelag_data(data,targetIndex,fre=fre)$newdata
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
        perform[[k]] <- oneModel(tmpnewdata,per=per,fre=fre,sflag=sflag)
        k <- k+1
        pseudoR <- pseudoPredict(tmpnewdata,per,targetIndex)
        perform[[k]] <- pseudoR
        
        perform
}

oneModel <- function(tmpnewdata,per=per,fre=fre,sflag=sflag){
        
        #### output 
        R2 <- 1:per ### R^2 values
        preds <- 1:per ### predict values
        residuals <- 1:per
        para <- c()
        
        n2 <- nrow(data) ### input data length
        sub <- 1 ### default sub-index of y 
        labs0 <- rownames(data)      
        
        if(sflag==1) pdq <- c(1,1,2,0,0,0,1) #day
        if(sflag==2) pdq <- c(1,1,1,2,1,0,7) #week
        if(sflag==3) pdq <- c(0,1,1,0,1,1,11) #month
        if(sflag==4) pdq <- c(0,1,0,0,1,0,2) #year
        
        for(n1 in (n2-per):(n2-1)){
                if(sflag<4){
                        #pdq <- sarima_paraNew(data[1:n1,sub],fre=fre)
                        stsr <- arima(ts(data[1:n1,sub],frequency = pdq[7]),order=pdq[1:3],seasonal = list(order=pdq[4:6],period=pdq[7]))
                        tmpdata <- data[1:n1,]
                        tmpdata[,sub] <- stsr$residuals
                        xnew <-stepCV_hq(tmpdata,sub,cvf=1,dir="forward")
                        mlmr <- lm(paste(colnames(tmpdata)[sub]," ~ ",paste(xnew,sep="",collapse = "+"),sep=""),data=as.data.frame(tmpdata))
                        res1 <- mlmr$residuals
                        
                }
                
                if(sflag==4){
                        tmpdata <- as.data.frame(data[1:n1,])
                        xnew <-stepCV_hq(tmpdata,sub,cvf=1,dir="forward")
                        mlmr <- lm( paste(colnames(tmpdata)[sub]," ~ ",paste(xnew,sep="",collapse = "+"),sep=""),data=tmpdata )
                        res <- ts(mlmr$residuals,frequency = fre)
                        pdq <- sarima_paraNew(res,fre=fre)
                        stsr <- arima(res,order=pdq[1:3],seasonal = list(order=pdq[4:6],period=pdq[7]))
                        res1 <- stsr$residuals
                }
                
                #### one step prediction
                if(sflag<4){
                        oneP <- pdq
                        preds[n1-(n2-per-1)] <- sum(c(1,data[n1+1,xnew]) * tof(as.vector(mlmr$coefficients),na.rm=FALSE,fill=TRUE,f="zero"))  +  predict(stsr, n.ahead=1)$pred[1]
                }else if(sflag==4){
                        oneP <- pdq
                        preds[n1-(n2-per-1)] <- sum(c(1,data[n1+1,xnew]) * tof(as.vector(mlmr$coefficients),na.rm=FALSE,fill=TRUE,f="zero"))  +  predict(stsr, n.ahead=1)$pred[1]
                }
                
                R2[n1-(n2-per-1)] <- R_squared_hq(data[1:n1,sub],data[1:n1,sub]-res1)
                residuals[n1-(n2-per-1)] <- data[n1+1,sub] - preds[n1-(n2-per-1)]
                para <- cbind(para,oneP)
        }
        
        #### plot predict result
        n1=n2-per
        ymax <- 1.1*max(c(data[(n1+1):n2,sub],as.vector(preds)))
        ymin <- 0.9*min(c(data[(n1+1):n2,sub],as.vector(preds)))
        labs=rownames(data)[(n1+1):n2]
        par(mai=c(1.2,1.2,1,1))
        plot(data[(n1+1):n2,sub],col=1,type="b",ylim=c(ymin,ymax),ylab="Price",xlab="",xaxt="n",lwd=3,main=main)
        lines(as.vector(preds),col=2,type="b",lwd=3)
        axis(1,at=1:(n2-n1),labels =FALSE)  
        
        pos <- 1:(n2-n1)-per/100
        if(length(pos) > 20){pos <- pos[seq(1,length(pos),length.out=20)];labs <- labs[seq(1,length(labs),length.out=20)];}
        text(pos, par("usr")[3]-0.11*(ymax-ymin), labels = labs, srt = 90, pos = 1, xpd = TRUE)
        legend(ord,legend=c("观察值","预测值"),col=1:2,lwd=2) 
        
        jobid <- jobTrace(jobid=9,trace)
        
        list(R2=R2,preds=preds,residuals=residuals,para=para)
        
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

