csv_clean: read raw data table

fill_NA: fill NA for each column

data_filling: fill all data indexes into day data

fill_onex: sub-function of data_filling, fill one column of data

groupPredict: group data into day, week, month or season

oneDimPredict: pre-deal data, considering time-lag and index length and so on

oneModel: sub-function of oneDimPredict, to run our model, predicting or training 

pseudoPredict: predicted by previous value only

sarima_paraNew: auto-estimate parameters of SARIAM model

stepCV_hq: forward/backward stepwise to select indexes

timelag_data: add time-lag for each index

delete_NA: get continous non-NA values

tof: the replace function of as.numeric

fraction_NA: delete index with NA fraction >= pNA, default 0.5

R_squared_hq: get R^2 value

acf2_hq: replace acf2 function without plot

jobTrace: job tracing and write log file
