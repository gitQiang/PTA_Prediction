#' PTA Price Prediction based on linear regression and ARIMA.
#' 
#' More details please contact: Qiang Huang. 
#' Copyright reserved by Qiang Huang.
#' 
#'====================================================
#' Qiang Huang, 2016/07/01, huangqiang@3golden.com.cn
#'====================================================
PTA_Prediction <- function(filenames=NULL,trace1=0,trans=0){
        # options
        #   filenames: filenames of add new indexes.
        #   trace1: print or not the running trace, 1 print, 0 not print
        #   trans: transformation for target index, 0 raw price, 1 diff(log,1)
        
        # output
        #   by now, we directly ouput to shiny web site.
        
        #==========================================================
        rm(list=ls())
        gc()
        filenames=NULL
        trace1=1
        trans <- 0
        setwd("D:/code/PTA_Prediction") ## only for testing ...
        
        # library and source files
        source("misc.R")
        jobid <- jobTrace(11,trace1);
        jobid <- jobTrace(12,trace1);
        library(lubridate) ## as.Date
        library(zoo) ## na.approx
        library(forecast) ##CV
        
        
        #==========================================================
        # default parameter settings
        startYear <- 2002
        endYear <- year(Sys.Date())
        useYears <- startYear:endYear
        
        
        #==========================================================
        # ? upload new indexes
        if(!is.null(filenames)){dataM <- addNewIndex(filenames,useYears);}else{dataM <- c();}
        
        jobid <- jobTrace(1,trace1)
        data <- data_filling()
        #load("data_2002_2015")
        data <- cbind(data,dataM)
        mode(data) <- "numeric"

        # specific setting to PTA Price
        target <- "市场价.PTA.对苯二甲酸."
        sub <- which(colnames(data)==target)
        tmp <- data[,-sub]
        data <- cbind(data[,sub],tmp)
        colnames(data)[1] <- "Target"
        
        # data index transform
        # tmp <- data_transform(data[,2:ncol(data)])
        # colnames(tmp) <- colnames(data)[2:ncol(data)]
        # data <- cbind(data[,1],tmp)
        
        #==========================================================
        ## one day, one week, one month and one quarter predictions
        fres <- c(1,7,12,4)
        pers <- c(50,30,20,10)
        precs <- list()
        
        results <- list()
        
        for(i in 1:4){
                jobtmp <- jobTrace(i+1,trace1)
                tmpdata <- groupPredict(data,i)
                tmpdata1 <- sapply(1:ncol(tmpdata), function(i) na.approx(tof(tmpdata[,i]),maxgap=5,na.rm=FALSE))
                colnames(tmpdata1) <- colnames(tmpdata)
                rownames(tmpdata1) <- rownames(tmpdata)
                tmpdata1 <- fraction_NA(tmpdata1,pNA=0.5)
                
                # target transform
                if(trans==1){
                        target0 <- tmpdata1[,1]
                        tmpsub <- which(!is.na(tmpdata1[,1]))
                        x1 <- Target_transform(tmpdata1[tmpsub,1])
                        tmpdata1[tmpsub[-1],1] <- x1
                        tmpdata1[tmpsub[1],1] <- NA
                }
                
                results[[i]] <- oneDimPredict(tmpdata1,targetIndex=1,fre=fres[i],per=pers[i],sflag=i,trace1=trace1,trans=trans)
                
                if(trans==1){
                        ind0 <- min(match(tmp[[1]]$labs,rownames(tmpdata1)))
                        obs0 <- Target_transform(tmp[[1]]$obs, target0[ind0-1])
                        pred0 <- Target_transform(tmp[[1]]$preds, target0[ind0-1])
                        plot_testing(obs0,pred0,tmp[[1]]$labs)
                }
                precs[[i]] <- precision_pred(results[[i]],p=0.05)
        }
        
        jobid <- jobTrace(10,trace1)
        
}

### 预留接口，用于获取实时数据
getUpdateData <- function(){}

### 预留接口，用于增加新的数据指标
addNewIndex <- function(filenames=NULL,useYears=2002:2016){
        
        startYear <- useYears[1]
        endYear <- useYears[length(useYears)]
        
        useDates <- seq.Date(as.Date(paste(startYear,"-1-1",sep="")),as.Date(paste(endYear,"-12-31",sep="")),by="day")
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
