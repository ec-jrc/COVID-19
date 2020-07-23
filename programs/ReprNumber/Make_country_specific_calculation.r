##### Code for estimating country specific effective reproduction number with different methods  ####
#####  parts of the code reutilizes already partly written code by othe researcrhes, the sources for these are:
#####https://www.datacamp.com/community/tutorials/replicating-in-r-covid19
##### https://www.r-bloggers.com/effective-reproduction-number-estimation/
#####
######
######
##### 

### instructions to run the code  ####
### 1. set the homedirectory in your computer where to save the results
##   ADD HOMEDIRECTORY, for example "c:/DATA/COVID/RESULTS_CALCULATION"
Homedir<-"ADD DIRECTORY TO SAVE RESULTS HERE"

### 2. Uncomment the package installations if they are not installed

# install.packages("purrr")
# install.packages("knitr")
#install.packages("kableExtra")
# install.packages("viridis")
# install.packages("accelerometry")
# install.packages("tidyverse")
# install.packages("R0")
#install.packages("IRdisplay")
#install.packages("smoother")
#install.packages("HDInterval")

####### 3.   Specify the country the analysis should be performed for   #####
Country_Name<-"France"

##### 4.  Specify the date the analysis should start
starting_date <- "2020-03-22"

### 5. Execute all the code, outputs will be saved in the homedirectory as a csv file with name "Reproduction_number_estimates.csv"


#### Possible errors    #####
##### some error messages are displayed during the analysis, these are no real errors but mostly messages from different
##### estimation tools on assumptions taken or if some parts of the analysis takes longer than expected
##### only if the analysis stops and it not finished there can be an error. This error is due to the country specific data
##### this can either depend tht the starting date is taken too early, some methods do not work in this situation. Other possibility 
#### is that in some countries zero cases can be reported for a longer timeperiod.
############################






suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(kableExtra))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(accelerometry))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(IRdisplay))
suppressPackageStartupMessages(library(smoother))
suppressPackageStartupMessages(library(HDInterval))

starting_date <- "2020-03-22"
Country_Name<-"France"

### Set up distribution on generation time
GT_pmf <- structure( c(0, 0.1, 0.1, 0.2, 0.2, 0.2, 0.1, 0.1), names=0:7)
GT_obj <- R0::generation.time("empirical", val=GT_pmf)

GT_obj<-R0::generation.time("gamma", c(6.6, 1.5))

Homedir<-"H:/COVID/R0_technical_paper/Graphs_test"

## Read in ECDC data   ##
data_ECDC <- read.csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", na.strings = "", fileEncoding = "UTF-8-BOM")


#####  Define functions to calculate effective reproduction number #######

####### Functions to use for effective reproduction number calculations   #######

##### Kevin systrom method   ########

est_rt_Systrom <- function(covid_cases, state_selected) {
  
  
  R_T_MAX = 12
  r_t_range = seq(0, R_T_MAX, length = R_T_MAX*100 + 1)
  GAMMA = 1/7
  
  #' Compute new cases and smooth them
  smooth_new_cases_1 <- function(cases){
    cases %>%
      arrange(date) %>%
      mutate(new_cases = c(cases[1], diff(cases))) %>%
      mutate(new_cases_smooth = round(
        smoother::smth.gaussian(new_cases, window = 7, alpha=7/4, tails = FALSE)
      )) %>%
      select(state, date, new_cases, new_cases_smooth)
  }
  
  
  #' Start using data first when smoothed value is above 25
  smooth_new_cases_25 <- function(x){
    start_index <- which(x$new_cases_smooth > 25)[1]
    length_vector <-dim(x)[1]
    return(x[start_index:length_vector,])
  }
  
  
  compute_likelihood <- function(cases){
    likelihood <- cases %>%
      filter(new_cases_smooth > 0) %>%
      mutate(
        r_t = list(r_t_range),
        lambda = map(lag(new_cases_smooth, 1), ~ .x * exp(GAMMA * (r_t_range - 1))),
        likelihood_r_t = map2(new_cases_smooth, lambda, dpois, log = TRUE)
      ) %>%
      slice(-1) %>%
      select(-lambda) %>%
      unnest(c(likelihood_r_t, r_t))
  }
  
  
  compute_posterior <- function(likelihood){
    likelihood %>%
      arrange(date) %>%
      group_by(r_t) %>%
      mutate(posterior = exp(
        zoo::rollapplyr(likelihood_r_t, 7, sum, partial = TRUE)
      )) %>%
      group_by(date) %>%
      mutate(posterior = posterior / sum(posterior, na.rm = TRUE)) %>%
      # HACK: NaNs in the posterior create issues later on. So we remove them.
      mutate(posterior = ifelse(is.nan(posterior), 0, posterior)) %>%
      ungroup() %>%
      select(-likelihood_r_t)
  }
  
  
  
  
  
  # Estimate R_t and a 95% highest-density interval around it
  estimate_rt <- function(posteriors){
    posteriors %>%
      group_by(state, date) %>%
      summarize(
        r_t_simulated = list(sample(r_t_range, 10000, replace = TRUE, prob = posterior)),
        r_t_most_likely = r_t_range[which.max(posterior)]
      ) %>%
      mutate(
        r_t_lo = map_dbl(r_t_simulated, ~ hdi(.x)[1]),
        r_t_hi = map_dbl(r_t_simulated, ~ hdi(.x)[2])
      ) %>%
      select(-r_t_simulated)
  }
  
  estimates <- covid_cases %>%
    filter(state == state_selected) %>%
    smooth_new_cases_1() %>%
    smooth_new_cases_25()  %>%
    compute_likelihood() %>%
    compute_posterior() %>%
    estimate_rt()
  
  return(estimates)
  
}


##### cislaghi method



est_rt_Cislaghi <- function(ts, window_interval) {
  data1 <- accelerometry::movingaves(ts, window=5)
  names(data1) <- names(ts)[3:(length(ts)-2)]
  res <- sapply( (1+window_interval):length(data1), function(t) {
    data1[t]/data1[t-window_interval]
  })
  return(res)
  
}


###### JRC method ###########################


est_rt_JRC <- function(ts, window_interval) {
  res <- sapply( 1:(length(ts)-window_interval), function(t) {
    ((log(ts[t+window_interval])-log(ts[t]))/window_interval)*7 + 1
  })
  return(res)
  
}

#### RKI method   ##############################

est_rt_rki <- function(ts, GT=4L) {
  # Sanity check
  if (!is.integer(GT) | !(GT>0)) stop("GT has to be postive integer.")
  # Estimate, if s=1 is t-7
  res <- sapply( (2*GT):length(ts), function(t) {
    sum(ts[t-(0:(GT-1))]) / sum(ts[t-2*GT+1:GT])
  })
  names(res) <- names(ts)[(2*GT):length(ts)]
  return(res)
}


###### Exponential growth, Wallinga & Lipsitch   ########


R.from.r <- function (r, GT) {
  Tmax = length(GT$GT)
  R = 1/sum(GT$GT * (exp(-r * (0:(Tmax - 1)))))
}

est_rt_exp <- function(ts, GT_obj, half_window_width=3L) {
  # Loop over all time points where the sliding window fits in
  res <- sapply((half_window_width+1):(length(ts)-half_window_width), function(t) {
    # Define the sliding window
    idx <- (t-half_window_width):(t+half_window_width)
    data <- data.frame(Date=1:length(ts), y=as.numeric(ts))
    # Fit a Poisson GLM
    m <- glm( y ~ 1 + Date, family=poisson, data= data  %>% slice(idx))
    # Extract the growth rate 
    r <- as.numeric(coef(m)["Date"])
    
    # Equation 2.9 from Wallinga & Lipsitch
    R <- R.from.r(r, GT_obj)
    return(R)
  })
  names(res) <- names(ts)[(half_window_width+1):(length(ts)-half_window_width)]
  return(res)
}

######### Wallinga & Teunis   ##############



est_rt_wt <- function(ts1, GT_obj) {
  end1 <- length(ts1) 
  
  R0::est.R0.TD(ts1, GT=GT_obj, begin=1, end=as.numeric(end1),correct=TRUE, nsim=1000)
}





#### function to calculate country specific value


Est_R0 <- function(start_date, country_name) {



  data2<-data_ECDC[,c("dateRep", "cases", "countriesAndTerritories")]
  # Set cases with minus value to 0
  data2$cases[which(data2$cases <0)]<-0
  # change format of Dates 
  data2$dateRep <- as.Date(as.character(data2$dateRep), format="%d/%m/%Y")
  names(data2) <- c("Date", "y", "CountryName")
  
  data2<-data_ECDC[,c("dateRep", "cases", "countriesAndTerritories")]
  # Set cases with minus value to 0
  data2$cases[which(data2$cases <0)]<-0
  # change format of Dates 
  data2$dateRep <- as.Date(as.character(data2$dateRep), format="%d/%m/%Y")
  names(data2) <- c("Date", "y", "CountryName")
  eval(parse(text=paste("index1 <- which(data2$CountryName ==\"",country_name, "\")", sep="" )))
  out <- data2[index1,]
  # Get starting date for the country
  out<-out[dim(out)[1]:1,]

  out <- out[out$Date %in% seq(from=as.Date(start_date) , to=out$Date[length(out$Date)] , by=1),]
  
  
out3 <- out
    names(out3)<-c("date", "cases1", "state")
  
  out3$cases[1] <-out3$cases1[1]
  for (j in 2:dim(out3)[1]) {
    out3$cases[j] <-   out3$cases[j-1]+  out3$cases1[j]
  }
  
  
  #### cislaghi method  ######
  
  
  rt_Cislaghi_3 <- est_rt_Cislaghi(out$y  %>% setNames(out$Date), window_interval=3)
  rt_Cislaghi_4 <- est_rt_Cislaghi(out$y  %>% setNames(out$Date), window_interval=4)
  rt_Cislaghi_5 <- est_rt_Cislaghi(out$y  %>% setNames(out$Date), window_interval=5)
  rt_Cislaghi_6 <- est_rt_Cislaghi(out$y  %>% setNames(out$Date), window_interval=6)
  
  rt_Cislaghi_3_df <- data.frame(Date=as.Date(names(rt_Cislaghi_3)), R_hat=rt_Cislaghi_3, Method=" 3")
  rt_Cislaghi_4_df <- data.frame(Date=as.Date(names(rt_Cislaghi_4)), R_hat=rt_Cislaghi_4, Method=" 4")
  rt_Cislaghi_5_df <- data.frame(Date=as.Date(names(rt_Cislaghi_5)), R_hat=rt_Cislaghi_5, Method=" 5")
  rt_Cislaghi_6_df <- data.frame(Date=as.Date(names(rt_Cislaghi_6)), R_hat=rt_Cislaghi_6, Method=" 6")
  
  
  #### Kevins Systrom method    ######
  
  # r_t_range is a vector of possible values for R_t
  R_T_MAX = 12
  r_t_range = seq(0, R_T_MAX, length = R_T_MAX*100 + 1)
  
  # Gamma is 1/serial interval
  # https://wwwnc.cdc.gov/eid/article/26/6/20-0357_article
  GAMMA = 1/7
  
  
  eval(parse(text=paste(" rt_Systrom <- est_rt_Systrom(covid_cases=out3, state_selected=\"",country_name, "\")", sep="" )))

  rt_Systrom_df <- data.frame(Date=rt_Systrom$date, R_hat=rt_Systrom$r_t_most_likely, upper=rt_Systrom$r_t_hi, lower=rt_Systrom$r_t_lo, Method="Systrom")
  

  ####   JRC method    #######
  
  
  
  rt_jrc <- est_rt_JRC(out$y  %>% setNames(out$Date), window_interval=7)
  
  rt_jrc_df <- data.frame(Date=as.Date(names(rt_jrc)), R_hat=rt_jrc, Method="JRC, GT=7")
  
  
  #############  RKI method ################
  
  
  # RKI method as specified in 2020-04-24 version of the article
  rt_rki <- est_rt_rki(out$y  %>% setNames(out$Date), GT=4L)
  
  rt_rki_df <- data.frame(Date=as.Date(names(rt_rki)), R_hat=rt_rki, Method="RKI, GT=4")
  
  
  
  
  ###### Wallinga & Lipsitch method    #############
  
  
  # Exponential growth with discrete Laplace transform of a GT \equiv = 4 distribution
  rt_exp2 <- est_rt_exp( out$y %>% setNames(out$Date), GT_obj=R0::generation.time("empirical", c(0,0,0,0,0,0,0,1)), half_window_width=3)
  
  # Exponential growth with discrete Laplace transform of the correct GT distribution
  rt_exp <- est_rt_exp( out$y %>% setNames(out$Date), GT_obj=GT_obj, half_window_width=3)
  
  rt_exp_df <- data.frame(Date=as.Date(names(rt_exp)), R_hat=rt_exp, Method="Stochastic generation time")
  #rt_exp_df <- data.frame(Date=as.Date(names(rt_exp_int[1:dim(rt_exp_int)[2]])), R_hat=rt_exp_int[1,], lower= rt_exp_int[2,], upper=rt_exp_int[3,], Method="Exp-Growth, random GT")
  
  rt_exp2_df <- data.frame(Date=as.Date(names(rt_exp2)), R_hat=rt_exp2, Method="Generation time 7 days (constant)")
  
  
  ### Wallinga and Teuniss############
  
  rt_wt <- est_rt_wt( out$y, GT=GT_obj)
  
  rt_wt_df <- cbind(Date=out$Date[as.numeric(names(rt_wt$R))], R_hat=rt_wt$R, rt_wt$conf.int, Method="W & T, ramdom GT")
  
  ### Create dataset with all estmates
  
  
  
  out_data <- data.frame(Date = as.Date(out$Date),Cases = out$y) 

  # Merge the data

names(rt_Cislaghi_3_df)[2] <- "Cislaghi_3"  
out_data <-  merge(x=out_data, y=rt_Cislaghi_3_df[,1:2], by = "Date", all.x=TRUE)
names(rt_Cislaghi_4_df)[2] <- "Cislaghi_4"  
out_data <-  merge(x=out_data, y=rt_Cislaghi_4_df[,1:2], by = "Date", all.x=TRUE)
names(rt_Cislaghi_5_df)[2] <- "Cislaghi_4"  
out_data <-  merge(x=out_data, y=rt_Cislaghi_5_df[,1:2], by = "Date", all.x=TRUE)
names(rt_Cislaghi_6_df)[2] <- "Cislaghi_4"  
out_data <-  merge(x=out_data, y=rt_Cislaghi_6_df[,1:2], by = "Date", all.x=TRUE)
names(rt_Systrom_df)[2:4] <- c("Systrom", "Systrom_upper", "Systrom_lower")  
out_data <-  merge(x=out_data, y=rt_Systrom_df[,1:4], by = "Date", all.x=TRUE)
names(rt_jrc_df)[2] <- "JRC"  
out_data <-  merge(x=out_data, y=rt_jrc_df[,1:2], by = "Date", all.x=TRUE)
names(rt_rki_df)[2] <- "RKI"  
out_data <-  merge(x=out_data, y=rt_rki_df[,1:2], by = "Date", all.x=TRUE)
names(rt_exp_df)[2] <- "Exp_growth_random_GT"  
out_data <-  merge(x=out_data, y=rt_exp_df[,1:2], by = "Date", all.x=TRUE)
names(rt_exp2_df)[2] <- "Exp_growth_constant_GT"  
out_data <-  merge(x=out_data, y=rt_exp2_df[,1:2], by = "Date", all.x=TRUE)
names(rt_wt_df)[2:4] <- c("Wallinga_Teunis", "Wallinga_Teunis_Upper", "Wallinga_Teunis_lower")  
out_data <-  merge(x=out_data, y=rt_wt_df[,1:4], by = "Date", all.x=TRUE)

return(out_data)

}


R_estimates <- Est_R0(start_date=starting_date, country_name=Country_Name)

write.csv(R_estimates, file.path(Homedir, "Reproduction_number_estimates.csv"))
