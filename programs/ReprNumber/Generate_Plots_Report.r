##### Code used for generation graphs used in JRC article on effective reproduction number estimating   ####
#####  parts of the code reutilizes already partly written code by othe researcrhes, the sources for these are:
#####https://www.datacamp.com/community/tutorials/replicating-in-r-covid19
##### https://www.r-bloggers.com/effective-reproduction-number-estimation/
#####
######
######
##### 

### instructions to run the code  ####
### 1. set the homedirectory in your computer where to save the results
### 2. Uncomment the package installations if they are not installed
### 2. Execute the code, all output will be saved in the homedirectory

##   ADD HOMEDIRECTORY, for example "c:/DATA/COVID/RESULTS_CALCULATION"
Homedir<-"ADD DIRECTORY TO SAVE RESULTS HERE"

# uncomment if the packages are not installed
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


suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(kableExtra))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(accelerometry))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(IRdisplay))
suppressPackageStartupMessages(library(smoother))
suppressPackageStartupMessages(library(HDInterval))


country_startdate <- function(Country_Name) {
if(Country_Name == "Italy") return("2020-02-20")  
if(Country_Name == "Germany") return("2020-02-26") 
  if(Country_Name == "Sweden") return("2020-02-27") 
  if(Country_Name == "Hungary") return("2020-03-05") 
  if(Country_Name == "Poland") return("2020-03-04") 
}

### Set up distribution on generation time
GT_pmf <- structure( c(0, 0.1, 0.1, 0.2, 0.2, 0.2, 0.1, 0.1), names=0:7)
GT_obj <- R0::generation.time("empirical", val=GT_pmf)

GT_obj<-R0::generation.time("gamma", c(6.6, 1.5))



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






#### cislaghi method  ######

data2<-data_ECDC[,c("dateRep", "cases", "countriesAndTerritories")]
# Set cases with minus value to 0
data2$cases[which(data2$cases <0)]<-0
# change format of Dates 
data2$dateRep <- as.Date(as.character(data2$dateRep), format="%d/%m/%Y")
names(data2) <- c("Date", "y", "CountryName")
out <- data2[which(data2$CountryName == "Italy"),]
# Get starting date for the country
out<-out[dim(out)[1]:1,]
start_date <- as.Date(country_startdate("Italy"))
out <- out[out$Date %in% seq(from=start_date , to=out$Date[length(out$Date)] , by=1),]




rt_Cislaghi_3 <- est_rt_Cislaghi(out$y  %>% setNames(out$Date), window_interval=3)
rt_Cislaghi_4 <- est_rt_Cislaghi(out$y  %>% setNames(out$Date), window_interval=4)
rt_Cislaghi_5 <- est_rt_Cislaghi(out$y  %>% setNames(out$Date), window_interval=5)
rt_Cislaghi_6 <- est_rt_Cislaghi(out$y  %>% setNames(out$Date), window_interval=6)

rt_Cislaghi_3_df <- data.frame(Date=as.Date(names(rt_Cislaghi_3)), R_hat=rt_Cislaghi_3, Method=" 3")
rt_Cislaghi_4_df <- data.frame(Date=as.Date(names(rt_Cislaghi_4)), R_hat=rt_Cislaghi_4, Method=" 4")
rt_Cislaghi_5_df <- data.frame(Date=as.Date(names(rt_Cislaghi_5)), R_hat=rt_Cislaghi_5, Method=" 5")
rt_Cislaghi_6_df <- data.frame(Date=as.Date(names(rt_Cislaghi_6)), R_hat=rt_Cislaghi_6, Method=" 6")
max_R <- 5
Normalizer <-  max(out$y) / max_R 
bar_data <- data.frame(Date = as.Date(out$Date),y_upd <- out$y/max(out$y) * max_R) 
out$Method <- "New cases"
cols = c("red", "green", "blue", "orange", "black")
### produce the plots   ####

png(file.path(Homedir, "Casliaghi_Italy.png"), width = 800, height = 480)
ggplot(rt_Cislaghi_4_df, aes(x=Date, y=R_hat, color=Method)) +  
  geom_line(data=rt_Cislaghi_3_df , lwd=1.5) + 
  geom_line(data=rt_Cislaghi_4_df , lwd=1.5) + 
  geom_line(data=rt_Cislaghi_5_df , lwd=1.5) + 
  geom_line(data=rt_Cislaghi_6_df , lwd=1.5) + 
  geom_line(data=out,
            aes(x = as.Date(Date),
                y = y / Normalizer), lwd=2)+
  
  scale_colour_manual(values = cols)+

  geom_hline(yintercept =1, lwd=1.5) + 
  ylab(expression(R(t))) +
  ggtitle("Cislaghi (on Italian data)")+
  scale_x_date(date_labels = "%b %d", breaks = "2 week")+

  scale_y_continuous(expand = c(0, 0), limits=c(0,max_R), sec.axis = sec_axis(trans= ~.*Normalizer,
                                                                              name = 'Number of new confirmed cases')) +
  theme(axis.text.x = element_text(angle=45, hjust = 1),legend.position="bottom",plot.title = element_text(size = 20),axis.line = element_line(colour = "black", size=1.5), axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 20),panel.background = element_blank())
dev.off()




#### Kevins Systrom method    ######

# r_t_range is a vector of possible values for R_t
R_T_MAX = 12
r_t_range = seq(0, R_T_MAX, length = R_T_MAX*100 + 1)

# Gamma is 1/serial interval
# https://wwwnc.cdc.gov/eid/article/26/6/20-0357_article
GAMMA = 1/7


data2<-data_ECDC[,c("dateRep", "cases", "countriesAndTerritories")]
# Set cases with minus value to 0
data2$cases[which(data2$cases <0)]<-0
# change format of Dates 
data2$dateRep <- as.Date(as.character(data2$dateRep), format="%d/%m/%Y")
names(data2) <- c("Date", "y", "CountryName")
out <- data2[which(data2$CountryName == "Germany"),]
# Get starting date for the country
out<-out[dim(out)[1]:1,]
start_date <- as.Date(country_startdate("Germany"))
out <- out[out$Date %in% seq(from=start_date , to=out$Date[length(out$Date)] , by=1),]
names(out)<-c("date", "cases1", "state")
out2<-out
out$cases[1] <-out$cases1[1]
for (j in 2:dim(out)[1]) {
 out$cases[j] <-   out$cases[j-1]+  out$cases1[j]
}



rt_Systrom <- est_rt_Systrom(covid_cases=out, state_selected="Germany") 



rt_Systrom_df <- data.frame(Date=rt_Systrom$date, R_hat=rt_Systrom$r_t_most_likely, upper=rt_Systrom$r_t_hi, lower=rt_Systrom$r_t_lo, Method="Systrom")

max_R <- 5
Normalizer <- max(out$cases1)   / max_R 
bar_data <- data.frame(Date = as.Date(out$date),y_upd <- out$cases1/max(out$cases1) * max_R) 
out$Method <- "Daily reported new cases"
cols = c("black", "blue")

### produce the plots   ####

png(file.path(Homedir, "Systrom_Germanyy.png"), width = 800, height = 480)
ggplot(rt_Systrom_df, aes(x=Date, y=R_hat, color=Method)) +  
  geom_line(lwd=1.5)+
  geom_ribbon(data=rt_Systrom_df, aes( ymin=lower, ymax=upper),inherit.aes=TRUE, fill="lightgrey", alpha=0.3) +
  geom_line(data=out,
            aes(x = as.Date(date),
                y= cases1 / Normalizer, color="New cases"), lwd=2)+
  
  scale_colour_manual(values = cols)+
  
  geom_hline(yintercept =1, lwd=1.5) + 
  ylab(expression(R(t))) +
  ggtitle("Systrom (on German data)")+
  scale_x_date(date_labels = "%b %d", breaks = "2 week")+
  
  scale_y_continuous(expand = c(0, 0), limits=c(0,max_R), sec.axis = sec_axis(trans= ~.*Normalizer,
                                                                              name = 'Number of new confirmed cases')) +
  theme(axis.text.x = element_text(angle=45, hjust = 1),legend.position="bottom",plot.title = element_text(size = 20),axis.line = element_line(colour = "black", size=1.5), axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 20),panel.background = element_blank())
dev.off()



####   JRC method    #######

data2<-data_ECDC[,c("dateRep", "cases", "countriesAndTerritories")]
# Set cases with minus value to 0
data2$cases[which(data2$cases <0)]<-0
# change format of Dates 
data2$dateRep <- as.Date(as.character(data2$dateRep), format="%d/%m/%Y")
names(data2) <- c("Date", "y", "CountryName")
out <- data2[which(data2$CountryName == "Germany"),]
# Get starting date for the country
out<-out[dim(out)[1]:1,]
start_date <- as.Date(country_startdate("Germany"))
out <- out[out$Date %in% seq(from=start_date , to=out$Date[length(out$Date)] , by=1),]




rt_jrc <- est_rt_JRC(out$y  %>% setNames(out$Date), window_interval=7)

rt_jrc_df <- data.frame(Date=as.Date(names(rt_jrc)), R_hat=rt_jrc, Method="JRC, GT=7")






max_R <- 5
Normalizer <-  max(out$y) / max_R 
bar_data <- data.frame(Date = as.Date(out$Date),y_upd <- out$y/max(out$y) * max_R) 
out$Method <- "New cases"

cols = c("red", "black")
### Make a plot of JRC method and epidemic curve  ###
png(file.path(Homedir, "JRC_Germany.png"), width = 800, height = 480)
ggplot(rt_jrc_df, aes(x=Date, y=R_hat, color=Method)) +  
  geom_line(data=rt_jrc_df,lwd=1.5)+ 
  geom_line(data=out,
            aes(x = as.Date(Date),
                y = y / Normalizer), lwd=2)+
  # coord_cartesian(ylim=c(0, ,max_R)) +
  geom_hline(yintercept =1, lwd=1.5) + 
  ylab(expression(R(t))) +
  scale_colour_manual(values = cols)+
  ggtitle("JRC method (on German data)")+
  scale_x_date(date_labels = "%b %d", breaks = "2 week")+
  
  scale_y_continuous(expand = c(0, 0), limits=c(0,max_R), sec.axis = sec_axis(trans= ~.*Normalizer,
                                                                              name = 'Number of new confirmed cases')) +
  theme(axis.text.x = element_text(angle=45, hjust = 1),legend.position="bottom",plot.title = element_text(size = 20),axis.line = element_line(colour = "black", size=1.5), axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 20),panel.background = element_blank())
dev.off()



#############  RKI method ################


data2<-data_ECDC[,c("dateRep", "cases", "countriesAndTerritories")]
# Set cases with minus value to 0
data2$cases[which(data2$cases <0)]<-0
# change format of Dates 
data2$dateRep <- as.Date(as.character(data2$dateRep), format="%d/%m/%Y")
names(data2) <- c("Date", "y", "CountryName")
out <- data2[which(data2$CountryName == "Germany"),]
# Get starting date for the country
out<-out[dim(out)[1]:1,]
start_date <- as.Date(country_startdate("Germany"))
out <- out[out$Date %in% seq(from=start_date , to=out$Date[length(out$Date)] , by=1),]



# RKI method as specified in 2020-04-24 version of the article
rt_rki <- est_rt_rki(out$y  %>% setNames(out$Date), GT=4L)

rt_rki_df <- data.frame(Date=as.Date(names(rt_rki)), R_hat=rt_rki, Method="RKI, GT=4")


max_R <- 5
Normalizer <-  max(out$y) / max_R 
bar_data <- data.frame(Date = as.Date(out$Date),y_upd <- out$y/max(out$y) * max_R) 
out$Method <- "New cases"

cols = c("black", "red")
### Make a plot of JRC method and epidemic curve  ###
png(file.path(Homedir, "RKI_method_Germany.png"), width = 800, height = 480)
ggplot(rt_rki_df, aes(x=Date, y=R_hat, color=Method)) +  
  geom_line(data=rt_rki_df,lwd=1.5)+ 
  geom_line(data=out,
            aes(x = as.Date(Date),
                y = y / Normalizer), lwd=2)+
  # coord_cartesian(ylim=c(0, ,max_R)) +
  geom_hline(yintercept =1, lwd=1.5) + 
  ylab(expression(R(t))) +
  scale_colour_manual(values = cols)+
  ggtitle("Robert Koch Institute method (on German data)")+
  scale_x_date(date_labels = "%b %d", breaks = "2 week")+
  
  scale_y_continuous(expand = c(0, 0), limits=c(0,max_R), sec.axis = sec_axis(trans= ~.*Normalizer,
                                                                              name = 'Number of new confirmed cases')) +
  theme(axis.text.x = element_text(angle=45, hjust = 1),legend.position="bottom",plot.title = element_text(size = 20),axis.line = element_line(colour = "black", size=1.5), axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 20),panel.background = element_blank())
dev.off()





###### Wallinga & Lipsitch method    #############

data2<-data_ECDC[,c("dateRep", "cases", "countriesAndTerritories")]
# Set cases with minus value to 0
data2$cases[which(data2$cases <0)]<-0
# change format of Dates 
data2$dateRep <- as.Date(as.character(data2$dateRep), format="%d/%m/%Y")
names(data2) <- c("Date", "y", "CountryName")
out <- data2[which(data2$CountryName == "Germany"),]
# Get starting date for the country
out<-out[dim(out)[1]:1,]
start_date <- as.Date(country_startdate("Germany"))
out <- out[out$Date %in% seq(from=start_date , to=out$Date[length(out$Date)] , by=1),]



# Exponential growth with discrete Laplace transform of a GT \equiv = 4 distribution
rt_exp2 <- est_rt_exp( out$y %>% setNames(out$Date), GT_obj=R0::generation.time("empirical", c(0,0,0,0,0,0,0,1)), half_window_width=3)

# Exponential growth with discrete Laplace transform of the correct GT distribution
rt_exp <- est_rt_exp( out$y %>% setNames(out$Date), GT_obj=GT_obj, half_window_width=3)

rt_exp_df <- data.frame(Date=as.Date(names(rt_exp)), R_hat=rt_exp, Method="Stochastic generation time")
#rt_exp_df <- data.frame(Date=as.Date(names(rt_exp_int[1:dim(rt_exp_int)[2]])), R_hat=rt_exp_int[1,], lower= rt_exp_int[2,], upper=rt_exp_int[3,], Method="Exp-Growth, random GT")

rt_exp2_df <- data.frame(Date=as.Date(names(rt_exp2)), R_hat=rt_exp2, Method="Generation time 7 days (constant)")


max_R <- 5
Normalizer <-  max(out$y) / max_R 
bar_data <- data.frame(Date = as.Date(out$Date),y_upd <- out$y/max(out$y) * max_R) 
out$Method <- "New cases"

cols = c("red", "black", "blue")
### Make a plot of JRC method and epidemic curve  ###
png(file.path(Homedir, "Wallinga_Lipsitch_germany.png"), width = 800, height = 480)
ggplot(rt_exp_df, aes(x=Date, y=R_hat, color=Method)) +  
  geom_line()+
  geom_line(data=rt_exp2_df , lwd=1.5) + 
  geom_line(data=out,
            aes(x = as.Date(Date),
                y = y / Normalizer), lwd=2)+
  # coord_cartesian(ylim=c(0, ,max_R)) +
  geom_hline(yintercept =1, lwd=1.5) + 
  ylab(expression(R(t))) +
  scale_colour_manual(values = cols)+
  ggtitle("Wallinga and Lipsitch (exponential growth method) (on German data)")+
  scale_x_date(date_labels = "%b %d", breaks = "2 week")+
  
  scale_y_continuous(expand = c(0, 0), limits=c(0,10), sec.axis = sec_axis(trans= ~.*Normalizer,
                                                                              name = 'Number of new confirmed cases')) +
  theme(axis.text.x = element_text(angle=45, hjust = 1),legend.position="bottom",plot.title = element_text(size = 20),axis.line = element_line(colour = "black", size=1.5), axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 20),panel.background = element_blank())
dev.off()



### Wallinga and Teuniss############


data2<-data_ECDC[,c("dateRep", "cases", "countriesAndTerritories")]
# Set cases with minus value to 0
data2$cases[which(data2$cases <0)]<-0
# change format of Dates 
data2$dateRep <- as.Date(as.character(data2$dateRep), format="%d/%m/%Y")
names(data2) <- c("Date", "y", "CountryName")
out <- data2[which(data2$CountryName == "Germany"),]
# Get starting date for the country
out<-out[dim(out)[1]:1,]
start_date <- as.Date(country_startdate("Germany"))
out <- out[out$Date %in% seq(from=start_date , to=out$Date[length(out$Date)] , by=1),]



rt_wt <- est_rt_wt( out$y, GT=GT_obj)

rt_wt_df <- cbind(Date=out$Date[as.numeric(names(rt_wt$R))], R_hat=rt_wt$R, rt_wt$conf.int, Method="W & T, ramdom GT")

max_R <- 5
Normalizer <-  max(out$y) / max_R 
bar_data <- data.frame(Date = as.Date(out$Date),y_upd <- out$y/max(out$y) * max_R) 
out$Method <- "New cases"

cols = c("black", "blue")

png(file.path(Homedir, "Wallinga_Teunis_Germany.png"), width = 800, height = 480)
ggplot(rt_wt_df, aes(x=Date, y=R_hat, color=Method)) +  
  geom_ribbon(data=rt_wt_df, aes( ymin=lower, ymax=upper),inherit.aes=TRUE, fill="lightgrey", alpha=0.6) +
    geom_line(lwd=0.5)+
  geom_line(data=out,
            aes(x = as.Date(Date),
                y= y / Normalizer, color="New cases"), lwd=2)+
  
scale_colour_manual(values = cols)+
  
  geom_hline(yintercept =1, lwd=1.5) + 
  ylab(expression(R(t))) +
  ggtitle("Wallinga and Teunis (on German data")+
  scale_x_date(date_labels = "%b %d", breaks = "2 week")+
  
  scale_y_continuous(expand = c(0, 0), limits=c(0,7), sec.axis = sec_axis(trans= ~.*Normalizer,
                                                                              name = 'Number of new confirmed cases')) +
  theme(axis.text.x = element_text(angle=45, hjust = 1),legend.position="bottom",plot.title = element_text(size = 20),axis.line = element_line(colour = "black", size=1.5), axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 20),panel.background = element_blank())
dev.off()













###### Apply different methods for Italy #######

data2<-data_ECDC[,c("dateRep", "cases", "countriesAndTerritories")]
# Set cases with minus value to 0
data2$cases[which(data2$cases <0)]<-0
# change format of Dates 
data2$dateRep <- as.Date(as.character(data2$dateRep), format="%d/%m/%Y")
names(data2) <- c("Date", "y", "CountryName")
out_sys <- data2[which(data2$CountryName == "Italy"),]
# Get starting date for the country
out_sys<-out_sys[dim(out_sys)[1]:1,]
start_date <- as.Date(country_startdate("Italy"))
out_sys <- out_sys[out_sys$Date %in% seq(from=start_date , to=out_sys$Date[length(out_sys$Date)] , by=1),]
names(out_sys)<-c("date", "cases1", "state")
out2<-out_sys
out_sys$cases[1] <-out_sys$cases1[1]
for (j in 2:dim(out_sys)[1]) {
  out_sys$cases[j] <-   out_sys$cases[j-1]+  out_sys$cases1[j]
}

data2<-data_ECDC[,c("dateRep", "cases", "countriesAndTerritories")]
# Set cases with minus value to 0
data2$cases[which(data2$cases <0)]<-0
# change format of Dates 
data2$dateRep <- as.Date(as.character(data2$dateRep), format="%d/%m/%Y")
names(data2) <- c("Date", "y", "CountryName")
out <- data2[which(data2$CountryName == "Italy"),]
# Get starting date for the country
out<-out[dim(out)[1]:1,]
start_date <- as.Date(country_startdate("Italy"))
out <- out[out$Date %in% seq(from=start_date , to=out$Date[length(out$Date)] , by=1),]

rt_Cislaghi_5 <- est_rt_Cislaghi(out$y  %>% setNames(out$Date), window_interval=5)
rt_Cislaghi_5_df <- data.frame(Date=as.Date(names(rt_Cislaghi_5)), R_hat=rt_Cislaghi_5, Method=" 5")
rt_Systrom <- est_rt_Systrom(covid_cases=out_sys, state_selected="Italy") 
rt_Systrom_df <- data.frame(Date=rt_Systrom$date, R_hat=rt_Systrom$r_t_most_likely, upper=rt_Systrom$r_t_hi, lower=rt_Systrom$r_t_lo, Method="Systrom")
rt_jrc <- est_rt_JRC(out$y  %>% setNames(out$Date), window_interval=7)
rt_jrc_df <- data.frame(Date=as.Date(names(rt_jrc)), R_hat=rt_jrc, Method="JRC, GT=7")
rt_rki <- est_rt_rki(out$y  %>% setNames(out$Date), GT=4L)
rt_rki_df <- data.frame(Date=as.Date(names(rt_rki)), R_hat=rt_rki, Method="RKI, GT=4")
rt_exp2 <- est_rt_exp( out$y %>% setNames(out$Date), GT_obj=R0::generation.time("empirical", c(0,0,0,0,0,0,0,1)), half_window_width=3)
rt_exp <- est_rt_exp( out$y %>% setNames(out$Date), GT_obj=GT_obj, half_window_width=3)
rt_exp_df <- data.frame(Date=as.Date(names(rt_exp)), R_hat=rt_exp, Method="Stochastic generation time")
rt_exp2_df <- data.frame(Date=as.Date(names(rt_exp2)), R_hat=rt_exp2, Method="Generation time 7 days (constant)")
rt_wt <- est_rt_wt( out$y, GT=GT_obj)
rt_wt_df <- cbind(Date=out$Date[as.numeric(names(rt_wt$R))], R_hat=rt_wt$R, rt_wt$conf.int, Method="W & T, ramdom GT")


max_R <- 5
Normalizer <-  max(out$y) / max_R 
bar_data <- data.frame(Date = as.Date(out$Date),y_upd <- out$y/max(out$y) * max_R) 
out$Method <- "New cases"



date_col <- c(as.Date(names(rt_exp2)), out$Date[as.numeric(names(rt_wt$R))],as.Date(names(rt_jrc)),  as.Date(names(rt_rki)), as.Date(names(rt_exp)), as.Date(names(rt_Cislaghi_5)), as.Date(rt_Systrom$date) )

R_col <- as.numeric(c(rt_exp2, rt_wt$R, rt_jrc, rt_rki, rt_exp, rt_Cislaghi_5, rt_Systrom_df$R_hat))

Method_col <- c( rep("Exp-Growth, GT = 7", length= length(rt_exp2)),  rep("W & T, ramdom GT", length= length(rt_wt$R)), rep("JRC, GT=7", length= length(rt_jrc)), rep("RKI, GT=4", length= length(rt_rki)), rep("Exp-Growth, random GT", length= length(rt_exp)), rep("Cislaghi, Window=5", length= length(rt_Cislaghi_5)),  rep("Systrom", length= dim(rt_Systrom)[1]))

Data_Total <- data.frame(Date=date_col, R_hat=R_col, Method=Method_col)

Data_total1 <- Data_Total[Data_Total$Method %in% c("Exp-Growth, GT = 7", "W & T, ramdom GT", "JRC, GT=7", "RKI, GT=4", "Exp-Growth, random GT", "Cislaghi, Window=5","Systrom"),]


### Plot the methods   ######
png(file.path(Homedir, "All_methods_Italy.png"), width = 800, height = 480)
ggplot(Data_total1, aes(x=Date, y=R_hat, group=Method)) +  
  
  geom_ribbon(data= rt_wt_df, aes( ymin=lower, ymax=upper),inherit.aes=TRUE, fill="lightgrey", alpha=0.5) +
  #geom_line(data=Ret_true, size=2, color="red") +  
  geom_line(aes(linetype=Method, color=Method, size=Method)) +
  scale_size_manual(values=c(1,1,1,1,1,1,1))+
    #scale_linetype_manual(values=c("twodash","longdash","dotdash","dashed","twodash","longdash","solid"))+
  coord_cartesian(ylim=c(0, 5)) +
  ylab(expression(R(t))) +

  scale_x_date(date_labels = "%b %d", breaks = "2 week")+
  geom_hline(yintercept=1, linetype="solid", color = "black")+
  ggtitle("R(t) estimate for Italy, grey area = 95% Wallinga & Teunis confidence interval")+
  theme(axis.text.x = element_text(angle=45, hjust = 1),legend.position="bottom",plot.title = element_text(size = 20),axis.line = element_line(colour = "black", size=1.5), axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 20),panel.background = element_blank())

dev.off()

png(file.path(Homedir, "Epidemic_curve_Italy.png"), width = 800, height = 480)
ggplot(out, aes(x=Date, y=y)) +
  geom_line(data=out , lwd=1.5) +
     ylab("Number of new confirmed cases") +
  ggtitle("Epidemic curve for Italy")+
  scale_x_date(date_labels = "%b %d")+
  theme(axis.text.x = element_text(angle=45, hjust = 1),legend.position="bottom",plot.title = element_text(size = 20),axis.line = element_line(colour = "black", size=1.5), axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 10),panel.background = element_blank())

dev.off()


###### Apply different methods for Germany #######

data2<-data_ECDC[,c("dateRep", "cases", "countriesAndTerritories")]
# Set cases with minus value to 0
data2$cases[which(data2$cases <0)]<-0
# change format of Dates 
data2$dateRep <- as.Date(as.character(data2$dateRep), format="%d/%m/%Y")
names(data2) <- c("Date", "y", "CountryName")
out_sys <- data2[which(data2$CountryName == "Germany"),]
# Get starting date for the country
out_sys<-out_sys[dim(out_sys)[1]:1,]
start_date <- as.Date(country_startdate("Germany"))
out_sys <- out_sys[out_sys$Date %in% seq(from=start_date , to=out_sys$Date[length(out_sys$Date)] , by=1),]
names(out_sys)<-c("date", "cases1", "state")
out2<-out_sys
out_sys$cases[1] <-out_sys$cases1[1]
for (j in 2:dim(out_sys)[1]) {
  out_sys$cases[j] <-   out_sys$cases[j-1]+  out_sys$cases1[j]
}

data2<-data_ECDC[,c("dateRep", "cases", "countriesAndTerritories")]
# Set cases with minus value to 0
data2$cases[which(data2$cases <0)]<-0
# change format of Dates 
data2$dateRep <- as.Date(as.character(data2$dateRep), format="%d/%m/%Y")
names(data2) <- c("Date", "y", "CountryName")
out <- data2[which(data2$CountryName == "Germany"),]
# Get starting date for the country
out<-out[dim(out)[1]:1,]
start_date <- as.Date(country_startdate("Germany"))
out <- out[out$Date %in% seq(from=start_date , to=out$Date[length(out$Date)] , by=1),]

rt_Cislaghi_5 <- est_rt_Cislaghi(out$y  %>% setNames(out$Date), window_interval=5)
rt_Cislaghi_5_df <- data.frame(Date=as.Date(names(rt_Cislaghi_5)), R_hat=rt_Cislaghi_5, Method=" 5")
rt_Systrom <- est_rt_Systrom(covid_cases=out_sys, state_selected="Germany") 
rt_Systrom_df <- data.frame(Date=rt_Systrom$date, R_hat=rt_Systrom$r_t_most_likely, upper=rt_Systrom$r_t_hi, lower=rt_Systrom$r_t_lo, Method="Systrom")
rt_jrc <- est_rt_JRC(out$y  %>% setNames(out$Date), window_interval=7)
rt_jrc_df <- data.frame(Date=as.Date(names(rt_jrc)), R_hat=rt_jrc, Method="JRC, GT=7")
rt_rki <- est_rt_rki(out$y  %>% setNames(out$Date), GT=4L)
rt_rki_df <- data.frame(Date=as.Date(names(rt_rki)), R_hat=rt_rki, Method="RKI, GT=4")
rt_exp2 <- est_rt_exp( out$y %>% setNames(out$Date), GT_obj=R0::generation.time("empirical", c(0,0,0,0,0,0,0,1)), half_window_width=3)
rt_exp <- est_rt_exp( out$y %>% setNames(out$Date), GT_obj=GT_obj, half_window_width=3)
rt_exp_df <- data.frame(Date=as.Date(names(rt_exp)), R_hat=rt_exp, Method="Stochastic generation time")
rt_exp2_df <- data.frame(Date=as.Date(names(rt_exp2)), R_hat=rt_exp2, Method="Generation time 7 days (constant)")
rt_wt <- est_rt_wt( out$y, GT=GT_obj)
rt_wt_df <- cbind(Date=out$Date[as.numeric(names(rt_wt$R))], R_hat=rt_wt$R, rt_wt$conf.int, Method="W & T, ramdom GT")


max_R <- 5
Normalizer <-  max(out$y) / max_R 
bar_data <- data.frame(Date = as.Date(out$Date),y_upd <- out$y/max(out$y) * max_R) 
out$Method <- "New cases"



date_col <- c(as.Date(names(rt_exp2)), out$Date[as.numeric(names(rt_wt$R))],as.Date(names(rt_jrc)),  as.Date(names(rt_rki)), as.Date(names(rt_exp)), as.Date(names(rt_Cislaghi_5)), as.Date(rt_Systrom$date) )

R_col <- as.numeric(c(rt_exp2, rt_wt$R, rt_jrc, rt_rki, rt_exp, rt_Cislaghi_5, rt_Systrom_df$R_hat))

Method_col <- c( rep("Exp-Growth, GT = 7", length= length(rt_exp2)),  rep("W & T, ramdom GT", length= length(rt_wt$R)), rep("JRC, GT=7", length= length(rt_jrc)), rep("RKI, GT=4", length= length(rt_rki)), rep("Exp-Growth, random GT", length= length(rt_exp)), rep("Cislaghi, Window=5", length= length(rt_Cislaghi_5)),  rep("Systrom", length= dim(rt_Systrom)[1]))

Data_Total <- data.frame(Date=date_col, R_hat=R_col, Method=Method_col)

Data_total1 <- Data_Total[Data_Total$Method %in% c("Exp-Growth, GT = 7", "W & T, ramdom GT", "JRC, GT=7", "RKI, GT=4", "Exp-Growth, random GT", "Cislaghi, Window=5","Systrom"),]


### Plot the methods   ######
png(file.path(Homedir, "All_methods_Germany.png"), width = 800, height = 480)
ggplot(Data_total1, aes(x=Date, y=R_hat, group=Method)) +  
  
  geom_ribbon(data= rt_wt_df, aes( ymin=lower, ymax=upper),inherit.aes=TRUE, fill="lightgrey", alpha=0.5) +
  #geom_line(data=Ret_true, size=2, color="red") +  
  geom_line(aes(linetype=Method, color=Method, size=Method)) +
  scale_size_manual(values=c(1,1,1,1,1,1,1))+
  #scale_linetype_manual(values=c("twodash","longdash","dotdash","dashed","twodash","longdash","solid"))+
  coord_cartesian(ylim=c(0, 5)) +
  ylab(expression(R(t))) +
  
  scale_x_date(date_labels = "%b %d", breaks = "2 week")+
  geom_hline(yintercept=1, linetype="solid", color = "black")+
  ggtitle("R(t) estimate for Germany, grey area = 95% Wallinga & Teunis confidence interval")+
  theme(axis.text.x = element_text(angle=45, hjust = 1),legend.position="bottom",plot.title = element_text(size = 20),axis.line = element_line(colour = "black", size=1.5), axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 20),panel.background = element_blank())

dev.off()

png(file.path(Homedir, "Epidemic_curve_Germany.png"), width = 800, height = 480)
ggplot(out, aes(x=Date, y=y)) +
  geom_line(data=out , lwd=1.5) +
  ylab("Number of new confirmed cases") +
  ggtitle("Epidemic curve for Germany")+
  scale_x_date(date_labels = "%b %d")+
  theme(axis.text.x = element_text(angle=45, hjust = 1),legend.position="bottom",plot.title = element_text(size = 20),axis.line = element_line(colour = "black", size=1.5), axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 10),panel.background = element_blank())

dev.off()



###### Apply different methods for Sweden #######

data2<-data_ECDC[,c("dateRep", "cases", "countriesAndTerritories")]
# Set cases with minus value to 0
data2$cases[which(data2$cases <0)]<-0
# change format of Dates 
data2$dateRep <- as.Date(as.character(data2$dateRep), format="%d/%m/%Y")
names(data2) <- c("Date", "y", "CountryName")
out_sys <- data2[which(data2$CountryName == "Sweden"),]
# Get starting date for the country
out_sys<-out_sys[dim(out_sys)[1]:1,]
start_date <- as.Date(country_startdate("Sweden"))
out_sys <- out_sys[out_sys$Date %in% seq(from=start_date , to=out_sys$Date[length(out_sys$Date)] , by=1),]
names(out_sys)<-c("date", "cases1", "state")
out2<-out_sys
out_sys$cases[1] <-out_sys$cases1[1]
for (j in 2:dim(out_sys)[1]) {
  out_sys$cases[j] <-   out_sys$cases[j-1]+  out_sys$cases1[j]
}

data2<-data_ECDC[,c("dateRep", "cases", "countriesAndTerritories")]
# Set cases with minus value to 0
data2$cases[which(data2$cases <0)]<-0
# change format of Dates 
data2$dateRep <- as.Date(as.character(data2$dateRep), format="%d/%m/%Y")
names(data2) <- c("Date", "y", "CountryName")
out <- data2[which(data2$CountryName == "Sweden"),]
# Get starting date for the country
out<-out[dim(out)[1]:1,]
start_date <- as.Date(country_startdate("Sweden"))
out <- out[out$Date %in% seq(from=start_date , to=out$Date[length(out$Date)] , by=1),]

rt_Cislaghi_5 <- est_rt_Cislaghi(out$y  %>% setNames(out$Date), window_interval=5)
rt_Cislaghi_5_df <- data.frame(Date=as.Date(names(rt_Cislaghi_5)), R_hat=rt_Cislaghi_5, Method=" 5")
rt_Systrom <- est_rt_Systrom(covid_cases=out_sys, state_selected="Sweden") 
rt_Systrom_df <- data.frame(Date=rt_Systrom$date, R_hat=rt_Systrom$r_t_most_likely, upper=rt_Systrom$r_t_hi, lower=rt_Systrom$r_t_lo, Method="Systrom")
rt_jrc <- est_rt_JRC(out$y  %>% setNames(out$Date), window_interval=7)
rt_jrc_df <- data.frame(Date=as.Date(names(rt_jrc)), R_hat=rt_jrc, Method="JRC, GT=7")
rt_rki <- est_rt_rki(out$y  %>% setNames(out$Date), GT=4L)
rt_rki_df <- data.frame(Date=as.Date(names(rt_rki)), R_hat=rt_rki, Method="RKI, GT=4")
rt_exp2 <- est_rt_exp( out$y %>% setNames(out$Date), GT_obj=R0::generation.time("empirical", c(0,0,0,0,0,0,0,1)), half_window_width=3)
rt_exp <- est_rt_exp( out$y %>% setNames(out$Date), GT_obj=GT_obj, half_window_width=3)
rt_exp_df <- data.frame(Date=as.Date(names(rt_exp)), R_hat=rt_exp, Method="Stochastic generation time")
rt_exp2_df <- data.frame(Date=as.Date(names(rt_exp2)), R_hat=rt_exp2, Method="Generation time 7 days (constant)")
rt_wt <- est_rt_wt( out$y, GT=GT_obj)
rt_wt_df <- cbind(Date=out$Date[as.numeric(names(rt_wt$R))], R_hat=rt_wt$R, rt_wt$conf.int, Method="W & T, ramdom GT")


max_R <- 5
Normalizer <-  max(out$y) / max_R 
bar_data <- data.frame(Date = as.Date(out$Date),y_upd <- out$y/max(out$y) * max_R) 
out$Method <- "New cases"



date_col <- c(as.Date(names(rt_exp2)), out$Date[as.numeric(names(rt_wt$R))],as.Date(names(rt_jrc)),  as.Date(names(rt_rki)), as.Date(names(rt_exp)), as.Date(names(rt_Cislaghi_5)), as.Date(rt_Systrom$date) )

R_col <- as.numeric(c(rt_exp2, rt_wt$R, rt_jrc, rt_rki, rt_exp, rt_Cislaghi_5, rt_Systrom_df$R_hat))

Method_col <- c( rep("Exp-Growth, GT = 7", length= length(rt_exp2)),  rep("W & T, ramdom GT", length= length(rt_wt$R)), rep("JRC, GT=7", length= length(rt_jrc)), rep("RKI, GT=4", length= length(rt_rki)), rep("Exp-Growth, random GT", length= length(rt_exp)), rep("Cislaghi, Window=5", length= length(rt_Cislaghi_5)),  rep("Systrom", length= dim(rt_Systrom)[1]))

Data_Total <- data.frame(Date=date_col, R_hat=R_col, Method=Method_col)

Data_total1 <- Data_Total[Data_Total$Method %in% c("Exp-Growth, GT = 7", "W & T, ramdom GT", "JRC, GT=7", "RKI, GT=4", "Exp-Growth, random GT", "Cislaghi, Window=5","Systrom"),]


### Plot the methods   ######
png(file.path(Homedir, "All_methods_Sweden.png"), width = 800, height = 480)
ggplot(Data_total1, aes(x=Date, y=R_hat, group=Method)) +  
  
  geom_ribbon(data= rt_wt_df, aes( ymin=lower, ymax=upper),inherit.aes=TRUE, fill="lightgrey", alpha=0.5) +
  #geom_line(data=Ret_true, size=2, color="red") +  
  geom_line(aes(linetype=Method, color=Method, size=Method)) +
  scale_size_manual(values=c(1,1,1,1,1,1,1))+
  #scale_linetype_manual(values=c("twodash","longdash","dotdash","dashed","twodash","longdash","solid"))+
  coord_cartesian(ylim=c(0, 5)) +
  ylab(expression(R(t))) +
  
  scale_x_date(date_labels = "%b %d", breaks = "2 week")+
  geom_hline(yintercept=1, linetype="solid", color = "black")+
  ggtitle("R(t) estimate for Sweden, grey area = 95% Wallinga & Teunis confidence interval")+
  theme(axis.text.x = element_text(angle=45, hjust = 1),legend.position="bottom",plot.title = element_text(size = 20),axis.line = element_line(colour = "black", size=1.5), axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 20),panel.background = element_blank())

dev.off()

png(file.path(Homedir, "Epidemic_curve_Sweden.png"), width = 800, height = 480)
ggplot(out, aes(x=Date, y=y)) +
  geom_line(data=out , lwd=1.5) +
  ylab("Number of new confirmed cases") +
  ggtitle("Epidemic curve for Sweden")+
  scale_x_date(date_labels = "%b %d")+
  theme(axis.text.x = element_text(angle=45, hjust = 1),legend.position="bottom",plot.title = element_text(size = 20),axis.line = element_line(colour = "black", size=1.5), axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 10),panel.background = element_blank())

dev.off()






###### Apply different methods for Hungary #######

data2<-data_ECDC[,c("dateRep", "cases", "countriesAndTerritories")]
# Set cases with minus value to 0
data2$cases[which(data2$cases <0)]<-0
# change format of Dates 
data2$dateRep <- as.Date(as.character(data2$dateRep), format="%d/%m/%Y")
names(data2) <- c("Date", "y", "CountryName")
out_sys <- data2[which(data2$CountryName == "Hungary"),]
# Get starting date for the country
out_sys<-out_sys[dim(out_sys)[1]:1,]
start_date <- as.Date(country_startdate("Hungary"))
out_sys <- out_sys[out_sys$Date %in% seq(from=start_date , to=out_sys$Date[length(out_sys$Date)] , by=1),]
names(out_sys)<-c("date", "cases1", "state")
out2<-out_sys
out_sys$cases[1] <-out_sys$cases1[1]
for (j in 2:dim(out_sys)[1]) {
  out_sys$cases[j] <-   out_sys$cases[j-1]+  out_sys$cases1[j]
}

data2<-data_ECDC[,c("dateRep", "cases", "countriesAndTerritories")]
# Set cases with minus value to 0
data2$cases[which(data2$cases <0)]<-0
# change format of Dates 
data2$dateRep <- as.Date(as.character(data2$dateRep), format="%d/%m/%Y")
names(data2) <- c("Date", "y", "CountryName")
out <- data2[which(data2$CountryName == "Hungary"),]
# Get starting date for the country
out<-out[dim(out)[1]:1,]
start_date <- as.Date(country_startdate("Hungary"))
out <- out[out$Date %in% seq(from=start_date , to=out$Date[length(out$Date)] , by=1),]

rt_Cislaghi_5 <- est_rt_Cislaghi(out$y  %>% setNames(out$Date), window_interval=5)
rt_Cislaghi_5_df <- data.frame(Date=as.Date(names(rt_Cislaghi_5)), R_hat=rt_Cislaghi_5, Method=" 5")
rt_Systrom <- est_rt_Systrom(covid_cases=out_sys, state_selected="Hungary") 
rt_Systrom_df <- data.frame(Date=rt_Systrom$date, R_hat=rt_Systrom$r_t_most_likely, upper=rt_Systrom$r_t_hi, lower=rt_Systrom$r_t_lo, Method="Systrom")
rt_jrc <- est_rt_JRC(out$y  %>% setNames(out$Date), window_interval=7)
rt_jrc_df <- data.frame(Date=as.Date(names(rt_jrc)), R_hat=rt_jrc, Method="JRC, GT=7")
rt_rki <- est_rt_rki(out$y  %>% setNames(out$Date), GT=4L)
rt_rki_df <- data.frame(Date=as.Date(names(rt_rki)), R_hat=rt_rki, Method="RKI, GT=4")
rt_exp2 <- est_rt_exp( out$y %>% setNames(out$Date), GT_obj=R0::generation.time("empirical", c(0,0,0,0,0,0,0,1)), half_window_width=3)
rt_exp <- est_rt_exp( out$y %>% setNames(out$Date), GT_obj=GT_obj, half_window_width=3)
rt_exp_df <- data.frame(Date=as.Date(names(rt_exp)), R_hat=rt_exp, Method="Stochastic generation time")
rt_exp2_df <- data.frame(Date=as.Date(names(rt_exp2)), R_hat=rt_exp2, Method="Generation time 7 days (constant)")
rt_wt <- est_rt_wt( out$y, GT=GT_obj)
rt_wt_df <- cbind(Date=out$Date[as.numeric(names(rt_wt$R))], R_hat=rt_wt$R, rt_wt$conf.int, Method="W & T, ramdom GT")


max_R <- 5
Normalizer <-  max(out$y) / max_R 
bar_data <- data.frame(Date = as.Date(out$Date),y_upd <- out$y/max(out$y) * max_R) 
out$Method <- "New cases"



date_col <- c(as.Date(names(rt_exp2)), out$Date[as.numeric(names(rt_wt$R))],as.Date(names(rt_jrc)),  as.Date(names(rt_rki)), as.Date(names(rt_exp)), as.Date(names(rt_Cislaghi_5)), as.Date(rt_Systrom$date) )

R_col <- as.numeric(c(rt_exp2, rt_wt$R, rt_jrc, rt_rki, rt_exp, rt_Cislaghi_5, rt_Systrom_df$R_hat))

Method_col <- c( rep("Exp-Growth, GT = 7", length= length(rt_exp2)),  rep("W & T, ramdom GT", length= length(rt_wt$R)), rep("JRC, GT=7", length= length(rt_jrc)), rep("RKI, GT=4", length= length(rt_rki)), rep("Exp-Growth, random GT", length= length(rt_exp)), rep("Cislaghi, Window=5", length= length(rt_Cislaghi_5)),  rep("Systrom", length= dim(rt_Systrom)[1]))

Data_Total <- data.frame(Date=date_col, R_hat=R_col, Method=Method_col)

Data_total1 <- Data_Total[Data_Total$Method %in% c("Exp-Growth, GT = 7", "W & T, ramdom GT", "JRC, GT=7", "RKI, GT=4", "Exp-Growth, random GT", "Cislaghi, Window=5","Systrom"),]


### Plot the methods   ######
png(file.path(Homedir, "All_methods_Hungary.png"), width = 800, height = 480)
ggplot(Data_total1, aes(x=Date, y=R_hat, group=Method)) +  
  
  geom_ribbon(data= rt_wt_df, aes( ymin=lower, ymax=upper),inherit.aes=TRUE, fill="lightgrey", alpha=0.5) +
  #geom_line(data=Ret_true, size=2, color="red") +  
  geom_line(aes(linetype=Method, color=Method, size=Method)) +
  scale_size_manual(values=c(1,1,1,1,1,1,1))+
  #scale_linetype_manual(values=c("twodash","longdash","dotdash","dashed","twodash","longdash","solid"))+
  coord_cartesian(ylim=c(0, 5)) +
  ylab(expression(R(t))) +
  
  scale_x_date(date_labels = "%b %d", breaks = "2 week")+
  geom_hline(yintercept=1, linetype="solid", color = "black")+
  ggtitle("R(t) estimate for Hungary, grey area = 95% Wallinga & Teunis confidence interval")+
  theme(axis.text.x = element_text(angle=45, hjust = 1),legend.position="bottom",plot.title = element_text(size = 20),axis.line = element_line(colour = "black", size=1.5), axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 20),panel.background = element_blank())

dev.off()

png(file.path(Homedir, "Epidemic_curve_Hungary.png"), width = 800, height = 480)
ggplot(out, aes(x=Date, y=y)) +
  geom_line(data=out , lwd=1.5) +
  ylab("Number of new confirmed cases") +
  ggtitle("Epidemic curve for Hungary")+
  scale_x_date(date_labels = "%b %d")+
  theme(axis.text.x = element_text(angle=45, hjust = 1),legend.position="bottom",plot.title = element_text(size = 20),axis.line = element_line(colour = "black", size=1.5), axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 10),panel.background = element_blank())

dev.off()



###### Apply different methods for Poland #######

data2<-data_ECDC[,c("dateRep", "cases", "countriesAndTerritories")]
# Set cases with minus value to 0
data2$cases[which(data2$cases <0)]<-0
# change format of Dates 
data2$dateRep <- as.Date(as.character(data2$dateRep), format="%d/%m/%Y")
names(data2) <- c("Date", "y", "CountryName")
out_sys <- data2[which(data2$CountryName == "Poland"),]
# Get starting date for the country
out_sys<-out_sys[dim(out_sys)[1]:1,]
start_date <- as.Date(country_startdate("Poland"))
out_sys <- out_sys[out_sys$Date %in% seq(from=start_date , to=out_sys$Date[length(out_sys$Date)] , by=1),]
names(out_sys)<-c("date", "cases1", "state")
out2<-out_sys
out_sys$cases[1] <-out_sys$cases1[1]
for (j in 2:dim(out_sys)[1]) {
  out_sys$cases[j] <-   out_sys$cases[j-1]+  out_sys$cases1[j]
}

data2<-data_ECDC[,c("dateRep", "cases", "countriesAndTerritories")]
# Set cases with minus value to 0
data2$cases[which(data2$cases <0)]<-0
# change format of Dates 
data2$dateRep <- as.Date(as.character(data2$dateRep), format="%d/%m/%Y")
names(data2) <- c("Date", "y", "CountryName")
out <- data2[which(data2$CountryName == "Poland"),]
# Get starting date for the country
out<-out[dim(out)[1]:1,]
start_date <- as.Date(country_startdate("Poland"))
out <- out[out$Date %in% seq(from=start_date , to=out$Date[length(out$Date)] , by=1),]

rt_Cislaghi_5 <- est_rt_Cislaghi(out$y  %>% setNames(out$Date), window_interval=5)
rt_Cislaghi_5_df <- data.frame(Date=as.Date(names(rt_Cislaghi_5)), R_hat=rt_Cislaghi_5, Method=" 5")
rt_Systrom <- est_rt_Systrom(covid_cases=out_sys, state_selected="Poland") 
rt_Systrom_df <- data.frame(Date=rt_Systrom$date, R_hat=rt_Systrom$r_t_most_likely, upper=rt_Systrom$r_t_hi, lower=rt_Systrom$r_t_lo, Method="Systrom")
rt_jrc <- est_rt_JRC(out$y  %>% setNames(out$Date), window_interval=7)
rt_jrc_df <- data.frame(Date=as.Date(names(rt_jrc)), R_hat=rt_jrc, Method="JRC, GT=7")
rt_rki <- est_rt_rki(out$y  %>% setNames(out$Date), GT=4L)
rt_rki_df <- data.frame(Date=as.Date(names(rt_rki)), R_hat=rt_rki, Method="RKI, GT=4")
rt_exp2 <- est_rt_exp( out$y %>% setNames(out$Date), GT_obj=R0::generation.time("empirical", c(0,0,0,0,0,0,0,1)), half_window_width=3)
rt_exp <- est_rt_exp( out$y %>% setNames(out$Date), GT_obj=GT_obj, half_window_width=3)
rt_exp_df <- data.frame(Date=as.Date(names(rt_exp)), R_hat=rt_exp, Method="Stochastic generation time")
rt_exp2_df <- data.frame(Date=as.Date(names(rt_exp2)), R_hat=rt_exp2, Method="Generation time 7 days (constant)")
rt_wt <- est_rt_wt( out$y, GT=GT_obj)
rt_wt_df <- cbind(Date=out$Date[as.numeric(names(rt_wt$R))], R_hat=rt_wt$R, rt_wt$conf.int, Method="W & T, ramdom GT")


max_R <- 5
Normalizer <-  max(out$y) / max_R 
bar_data <- data.frame(Date = as.Date(out$Date),y_upd <- out$y/max(out$y) * max_R) 
out$Method <- "New cases"



date_col <- c(as.Date(names(rt_exp2)), out$Date[as.numeric(names(rt_wt$R))],as.Date(names(rt_jrc)),  as.Date(names(rt_rki)), as.Date(names(rt_exp)), as.Date(names(rt_Cislaghi_5)), as.Date(rt_Systrom$date) )

R_col <- as.numeric(c(rt_exp2, rt_wt$R, rt_jrc, rt_rki, rt_exp, rt_Cislaghi_5, rt_Systrom_df$R_hat))

Method_col <- c( rep("Exp-Growth, GT = 7", length= length(rt_exp2)),  rep("W & T, ramdom GT", length= length(rt_wt$R)), rep("JRC, GT=7", length= length(rt_jrc)), rep("RKI, GT=4", length= length(rt_rki)), rep("Exp-Growth, random GT", length= length(rt_exp)), rep("Cislaghi, Window=5", length= length(rt_Cislaghi_5)),  rep("Systrom", length= dim(rt_Systrom)[1]))

Data_Total <- data.frame(Date=date_col, R_hat=R_col, Method=Method_col)

Data_total1 <- Data_Total[Data_Total$Method %in% c("Exp-Growth, GT = 7", "W & T, ramdom GT", "JRC, GT=7", "RKI, GT=4", "Exp-Growth, random GT", "Cislaghi, Window=5","Systrom"),]


### Plot the methods   ######
png(file.path(Homedir, "All_methods_Poland.png"), width = 800, height = 480)
ggplot(Data_total1, aes(x=Date, y=R_hat, group=Method)) +  
  
  geom_ribbon(data= rt_wt_df, aes( ymin=lower, ymax=upper),inherit.aes=TRUE, fill="lightgrey", alpha=0.5) +
  #geom_line(data=Ret_true, size=2, color="red") +  
  geom_line(aes(linetype=Method, color=Method, size=Method)) +
  scale_size_manual(values=c(1,1,1,1,1,1,1))+
  #scale_linetype_manual(values=c("twodash","longdash","dotdash","dashed","twodash","longdash","solid"))+
  coord_cartesian(ylim=c(0, 5)) +
  ylab(expression(R(t))) +
  
  scale_x_date(date_labels = "%b %d", breaks = "2 week")+
  geom_hline(yintercept=1, linetype="solid", color = "black")+
  ggtitle("R(t) estimate for Poland, grey area = 95% Wallinga & Teunis confidence interval")+
  theme(axis.text.x = element_text(angle=45, hjust = 1),legend.position="bottom",plot.title = element_text(size = 20),axis.line = element_line(colour = "black", size=1.5), axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 20),panel.background = element_blank())

dev.off()

png(file.path(Homedir, "Epidemic_curve_Poland.png"), width = 800, height = 480)
ggplot(out, aes(x=Date, y=y)) +
  geom_line(data=out , lwd=1.5) +
  ylab("Number of new confirmed cases") +
  ggtitle("Epidemic curve for Poland")+
  scale_x_date(date_labels = "%b %d")+
  theme(axis.text.x = element_text(angle=45, hjust = 1),legend.position="bottom",plot.title = element_text(size = 20),axis.line = element_line(colour = "black", size=1.5), axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 10),panel.background = element_blank())

dev.off()








#################    TEST MODELLED DATA   ####################



################### Test modelled data   ##############################


# Define time varying effective reproduction number
Ret1 <- function(date1) {
  if(date1 <= as.Date("2020-03-15")) {
    R_val <- 2.5
  }
  if (date1 > as.Date("2020-03-15") & date1 <= as.Date("2020-04-15")) {
    R_val <- 1.5
  }
  if (date1 > as.Date("2020-04-15"))  {
    R_val <- 0.75
  }
  return(R_val)
}
Ret<-Ret1

routbreak <- function(n=100, Ret, GT_obj, initial_cases = 10) {
  # Set up time series of incident cases
  y <- rep(0, n + length(GT_pmf))
  y[seq_len(length(initial_cases))] <- initial_cases
  # Outbreak starts on 2020-02-15
  dates <- as.Date("2020-02-15") + 0:(n-1)
  # Extract serial interval PMF, ignore support at 0.
  GT_pmf <- GT_obj$GT[-1]
  
  # Loop over all time points
  for (i in 1:n) {
    date <- dates[i]
    y[i + 1:length(GT_pmf)] <- y[i] * (Ret(date) * GT_pmf) + y[i + 1:length(GT_pmf)]
  }
  
  # Data frame with the result. Assume we start on 15th of Feb
  res <- data.frame(Date=dates, y=y[1:n])
  
  #Done
  return(res)
}


# Generate an outbreak (no stochasticity, just the difference equation)
out <- routbreak(n=100, Ret=Ret, GT_obj=GT_obj)
out <- out %>% mutate(ratio = y/lag(y))
# Data frame with the true values
True_val <- NULL
for (i in 1:length(out$Date)) {
  date <- out$Date[i]
  True_val <- c(True_val, Ret(date))
}
n1<-100
names(True_val) <- as.Date("2020-02-15") + 0:(n1-1)
Data_R_sim <- data.frame(Date=as.Date("2020-02-15") + 0:(n1-1), R_sim=True_val, Method="Simulated R(t)")

max_R <- 2.5
Normalizer <-  max(out$y) / max_R 
bar_data <- data.frame(Date = as.Date(out$Date),y_upd <- out$y/max(out$y) * max_R) 
out$Method <- "Daily reported new cases"
cols = c("black", "red")

png(file.path(Homedir, "Model_outbreak.png"), width = 800, height = 480)
ggplot(Data_R_sim, aes(x=Date, y=R_sim, color=Method)) + geom_line(lwd=1.5) +
  ylab("R(t)")+
  scale_colour_manual(values = cols)+
  geom_line(data=out,
            aes(x = as.Date(Date),
                y= y / Normalizer, color="New cases"), lwd=2)+
  scale_y_continuous(expand = c(0, 0), limits=c(0,max_R), sec.axis = sec_axis(trans= ~.*Normalizer,
                                                                              name = 'Daily new cases (simulated)')) +
  scale_x_date(date_labels = "%b %d", breaks = "2 week")+    theme(axis.text.x = element_text(angle=45, hjust = 1),legend.position="bottom",plot.title = element_text(size = 20),axis.line = element_line(colour = "black", size=1.5), axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 20),panel.background = element_blank())
dev.off()






# RKI method as specified in 2020-04-24 version of the article
rt_rki <- est_rt_rki(out$y  %>% setNames(out$Date), GT=4L)

# Exponential growth with discrete Laplace transform of a GT \equiv = 4 distribution
rt_exp2 <- est_rt_exp( out$y %>% setNames(out$Date), GT_obj=R0::generation.time("empirical", c(0,0,0,0,0,0,0,1)), half_window_width=3)

# Exponential growth with discrete Laplace transform of the correct GT distribution
rt_exp <- est_rt_exp( out$y %>% setNames(out$Date), GT_obj=GT_obj, half_window_width=3)


# Wallinga Teunis approach with correct GT distribution
rt_wt <- est_rt_wt( out$y, GT=GT_obj)

# Data frame with the true values
True_val <- NULL
for (i in 1:length(out$Date)) {
  date <- out$Date[i]
  True_val <- c(True_val, Ret(date))
}

Ret_true <- data.frame(Date=out$Date, R_hat = True_val, Method="True R(t)")


rt_jrc <- est_rt_JRC(out$y  %>% setNames(out$Date), window_interval=7)

# Cislaghi method
rt_Cislaghi_3 <- est_rt_Cislaghi(out$y  %>% setNames(out$Date), window_interval=3)
rt_Cislaghi_4 <- est_rt_Cislaghi(out$y  %>% setNames(out$Date), window_interval=4)
rt_Cislaghi_5 <- est_rt_Cislaghi(out$y  %>% setNames(out$Date), window_interval=5)
rt_Cislaghi_6 <- est_rt_Cislaghi(out$y  %>% setNames(out$Date), window_interval=6)

rt_jrc <- est_rt_JRC(out$y  %>% setNames(out$Date), window_interval=7)

out2 <- out
out2$state <- "simulated"
out2$date <- as.character(out2$Date)
out2$cases[1] <-out2$y[1]
for (j in 2:dim(out2)[1]) {
  out2$cases[j] <-   out2$cases[j-1]+  out2$y[j]
}
out2$cases <- round(out2$cases,0)

rt_Systrom <- est_rt_Systrom(covid_cases=out2, state_selected="simulated") 





#########
# Convert fits to unified data.frames
##########
rt_wt_df <- cbind(Date=out$Date[as.numeric(names(rt_wt$R))], R_hat=rt_wt$R, rt_wt$conf.int, Method="W & T, ramdom GT")

rt_rki_df <- data.frame(Date=as.Date(names(rt_rki)), R_hat=rt_rki, Method="RKI, GT=4")

rt_exp_df <- data.frame(Date=as.Date(names(rt_exp)), R_hat=rt_exp, Method="Exp-Growth, correct GT")

rt_exp2_df <- data.frame(Date=as.Date(names(rt_exp2)), R_hat=rt_exp2, Method="Exp-Growth, GT=7")

rt_Cislaghi_3_df <- data.frame(Date=as.Date(names(rt_Cislaghi_3)), R_hat=rt_Cislaghi_3, Method="Cislaghi, Window=3")
rt_Cislaghi_4_df <- data.frame(Date=as.Date(names(rt_Cislaghi_4)), R_hat=rt_Cislaghi_4, Method="Cislaghi, Window=4")
rt_Cislaghi_5_df <- data.frame(Date=as.Date(names(rt_Cislaghi_5)), R_hat=rt_Cislaghi_5, Method="Cislaghi, Window=5")
rt_Cislaghi_6_df <- data.frame(Date=as.Date(names(rt_Cislaghi_6)), R_hat=rt_Cislaghi_6, Method="Cislaghi, Window=6")


rt_jrc_df <- data.frame(Date=as.Date(names(rt_jrc)), R_hat=rt_jrc, Method="JRC, GT=7")
rt_Systrom_df <- data.frame(Date=as.Date(rt_Systrom$date), R_hat=rt_Systrom$r_t_most_likely,Method="Systrom")

########## Make a major dataframe
Ret_true <- data.frame(Date=out$Date, R_hat = True_val, Method="True R(t)")

date_col <- c(as.Date(names(rt_exp2)),out$Date, out$Date[as.numeric(names(rt_wt$R))],as.Date(names(rt_jrc)),  as.Date(names(rt_rki)), as.Date(names(rt_exp)), as.Date(names(rt_Cislaghi_4)), as.Date(names(rt_Cislaghi_5)),as.Date(names(rt_Cislaghi_6)))
R_col <- as.numeric(c(rt_exp2, True_val, rt_wt$R, rt_jrc, rt_rki, rt_exp, rt_Cislaghi_4, rt_Cislaghi_5, rt_Cislaghi_6))

Method_col <- c( rep("Exp-Growath, GT = 7", length= length(rt_exp2)), rep("True R(t)", length= length(True_val)), rep("W & T, ramdom GT", length= length(rt_wt$R)), rep("JRC, GT=7", length= length(rt_jrc)), rep("RKI, GT=4", length= length(rt_rki)), rep("Exp-Growth, random GT", length= length(rt_exp)), rep("Cislaghni, Window=4", length= length(rt_Cislaghi_4)), rep("Cislaghni, Window=5", length= length(rt_Cislaghi_5)), rep("Cislaghni, Window=6", length= length(rt_Cislaghi_6)))

Data_total <- data.frame(Date=date_col, R_hat=R_col, Method=Method_col)



date_col <- c(as.Date(names(rt_exp2)),out$Date, out$Date[as.numeric(names(rt_wt$R))],as.Date(names(rt_jrc)),  as.Date(names(rt_rki)), as.Date(names(rt_exp)), as.Date(names(rt_Cislaghi_3)),as.Date(names(rt_Cislaghi_4)), as.Date(names(rt_Cislaghi_5)),as.Date(names(rt_Cislaghi_6)), as.Date(rt_Systrom$date) )

R_col <- as.numeric(c(rt_exp2, True_val, rt_wt$R, rt_jrc, rt_rki, rt_exp, rt_Cislaghi_3, rt_Cislaghi_4, rt_Cislaghi_5, rt_Cislaghi_6,rt_Systrom_df$R_hat))

Method_col <- c( rep("Exp-Growth, GT = 7", length= length(rt_exp2)), rep("True R(t)", length= length(True_val)), rep("W & T, ramdom GT", length= length(rt_wt$R)), rep("JRC, GT=7", length= length(rt_jrc)), rep("RKI, GT=4", length= length(rt_rki)), rep("Exp-Growth, random GT", length= length(rt_exp)), rep("Cislaghi, Window=3", length= length(rt_Cislaghi_3)),rep("Cislaghi, Window=4", length= length(rt_Cislaghi_4)), rep("Cislaghi, Window=5", length= length(rt_Cislaghi_5)), rep("Cislaghi, Window=6", length= length(rt_Cislaghi_6)), rep("Systrom", length= dim(rt_Systrom)[1]))

Data_Total <- data.frame(Date=date_col, R_hat=R_col, Method=Method_col)

Data_total <- Data_Total[Data_Total$Method %in% c("Exp-Growth, GT = 7", "True R(t)", "W & T, ramdom GT", "JRC, GT=7", "RKI, GT=4", "Exp-Growth, random GT", "Cislaghi, Window=5","Systrom"),]


### Plot the methods   ######
png(file.path(Homedir, "R_e_Simulated.png"), width = 800, height = 480)
ggplot(Data_total, aes(x=Date, y=R_hat, group=Method)) +  
  
  geom_ribbon(data= rt_wt_df, aes( ymin=lower, ymax=upper),inherit.aes=TRUE, fill="lightgrey", alpha=0.5) +
  #geom_line(data=Ret_true, size=2, color="red") +  
  geom_line(aes(linetype=Method, color=Method, size=Method)) +
  scale_size_manual(values=c(1,1,1,1,1,1,2,1))+
  scale_linetype_manual(values=c("twodash","longdash","dotdash","dashed","twodash","longdash","solid","dotted"))+
  coord_cartesian(ylim=c(0, 5)) +
  ylab(expression(R(t))) +
  scale_x_date(date_labels = "%b %d", breaks = "2 week")+
  geom_hline(yintercept=1, linetype="solid", color = "black")+
  ggtitle("R(t) estimate with simulated data, grey area = 95% Wallinga & Teunis confidence interval")+
  theme(axis.text.x = element_text(angle=45, hjust = 1),legend.position="bottom",plot.title = element_text(size = 20),axis.line = element_line(colour = "black", size=1.5), axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 20),panel.background = element_blank())

dev.off()

######  Routint to plot the simulated data    ######

### Cislaghi


Data_cislaghi <- Data_Total[Data_Total$Method %in% c("True R(t)", "Cislaghi, Window=3", "Cislaghi, Window=4", "Cislaghi, Window=5", "Cislaghi, Window=6"),]



##### Plot the methods   ######
png(file.path(Homedir, "cislaghi, simulated.png"), width = 800, height = 480)
ggplot(Data_cislaghi, aes(x=Date, y=R_hat, group=Method)) +  
  
  
  geom_line(aes(linetype=Method, color=Method, size=Method)) +
  scale_size_manual(values=c(1,1,1,1,2))+
  scale_linetype_manual(values=c("twodash","longdash","dotdash","dashed","solid"))+
  coord_cartesian(ylim=c(0, 5)) +
  ylab(expression(R(t))) +
  scale_x_date(date_labels = "%b %d", breaks = "2 week")+
  geom_hline(yintercept=1, linetype="solid", color = "black")+
  ggtitle("Cislaghi method with simulated data")+
  
  theme(axis.text.x = element_text(angle=45, hjust = 1),legend.position="bottom",plot.title = element_text(size = 20),axis.line = element_line(colour = "black", size=1.5), axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 10),panel.background = element_blank())


dev.off()


### Systrom


Data_Systrom <- Data_Total[Data_Total$Method %in% c("True R(t)", "Systrom"),]







### Plot the methods   ######
png(file.path(Homedir, "Systrom_simulated.png"), width = 800, height = 480)
ggplot(Data_Systrom, aes(x=Date, y=R_hat, group=Method)) +  
  
  
  geom_line(aes(linetype=Method, color=Method, size=Method)) +
  scale_size_manual(values=c(1,2))+
  scale_linetype_manual(values=c("dashed","solid"))+
  coord_cartesian(ylim=c(0, 5)) +
  ylab(expression(R(t))) +
  scale_x_date(date_labels = "%b %d", breaks = "2 week")+
  geom_hline(yintercept=1, linetype="solid", color = "black")+
  ggtitle("Systrom method with simulated data")+
  theme(axis.text.x = element_text(angle=45, hjust = 1),legend.position="bottom",plot.title = element_text(size = 20),axis.line = element_line(colour = "black", size=1.5), axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 20),panel.background = element_blank())




dev.off()

### JRC


Data_JRC <- Data_Total[Data_Total$Method %in% c("True R(t)", "JRC, GT=7"),]



### Plot the methods   ######
png(file.path(Homedir, "JRC_simulated.png"), width = 800, height = 480)
ggplot(Data_JRC, aes(x=Date, y=R_hat, group=Method)) +  
  
  
  geom_line(aes(linetype=Method, color=Method, size=Method)) +
  scale_size_manual(values=c(1,2))+
  scale_linetype_manual(values=c("dashed","solid"))+
  coord_cartesian(ylim=c(0, 5)) +
  ylab(expression(R(t))) +
  scale_x_date(date_labels = "%b %d", breaks = "2 week")+
  geom_hline(yintercept=1, linetype="solid", color = "black")+
  ggtitle("JRC method with simulated data")+
  theme(axis.text.x = element_text(angle=45, hjust = 1),legend.position="bottom",plot.title = element_text(size = 20),axis.line = element_line(colour = "black", size=1.5), axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 20),panel.background = element_blank())


dev.off()


### RKI


Data_RKI <- Data_Total[Data_Total$Method %in% c("True R(t)", "RKI, GT=4"),]



### Plot the methods   ######
png(file.path(Homedir, "RKI_simulated.png"), width = 800, height = 480)
ggplot(Data_RKI, aes(x=Date, y=R_hat, group=Method)) +  
  
  
  geom_line(aes(linetype=Method, color=Method, size=Method)) +
  scale_size_manual(values=c(1,2))+
  scale_linetype_manual(values=c("dashed","solid"))+
  coord_cartesian(ylim=c(0, 5)) +
  ylab(expression(R(t))) +
  scale_x_date(date_labels = "%b %d", breaks = "2 week")+
  geom_hline(yintercept=1, linetype="solid", color = "black")+
  ggtitle("RKI method with simulated data")+
  theme(axis.text.x = element_text(angle=45, hjust = 1),legend.position="bottom",plot.title = element_text(size = 20),axis.line = element_line(colour = "black", size=1.5), axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 20),panel.background = element_blank())


dev.off()

### Wallinga & Lipsitch


Data_Expo <- Data_Total[Data_Total$Method %in% c("True R(t)", "Exp-Growth, random GT", "Exp-Growth, GT = 7"),]



### Plot the methods   ######
png(file.path(Homedir, "Exp_simulated.png"), width = 800, height = 480)
ggplot(Data_Expo, aes(x=Date, y=R_hat, group=Method)) +  
  
  
  geom_line(aes(linetype=Method, color=Method, size=Method)) +
  scale_size_manual(values=c(1,1, 2))+
  scale_linetype_manual(values=c("dashed","dotdash","solid"))+
  coord_cartesian(ylim=c(0, 5)) +
  ylab(expression(R(t))) +
  scale_x_date(date_labels = "%b %d", breaks = "2 week")+
  geom_hline(yintercept=1, linetype="solid", color = "black")+
  ggtitle("Exponential growth method with simulated data")+
  theme(axis.text.x = element_text(angle=45, hjust = 1),legend.position="bottom",plot.title = element_text(size = 20),axis.line = element_line(colour = "black", size=1.5), axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 20),panel.background = element_blank())


dev.off()

### Wallinga & Teunis


Data_Teunis <- Data_Total[Data_Total$Method %in% c("True R(t)", "W & T, ramdom GT"),]



### Plot the methods   ######
png(file.path(Homedir, "WT_simulated.png"), width = 800, height = 480)
ggplot(Data_Teunis, aes(x=Date, y=R_hat, group=Method)) +  
  
  
  geom_line(aes(linetype=Method, color=Method, size=Method)) +
  geom_ribbon(data= rt_wt_df, aes( ymin=lower, ymax=upper),inherit.aes=TRUE, fill="lightgrey", alpha=0.5) +
  scale_size_manual(values=c(2,1))+
  scale_linetype_manual(values=c("solid", "dashed"))+
  coord_cartesian(ylim=c(0, 5)) +
  ylab(expression(R(t))) +
  scale_x_date(date_labels = "%b %d", breaks = "2 week")+
  geom_hline(yintercept=1, linetype="solid", color = "black")+
  ggtitle("Wallinga & Teunis method with simulated data (shaded area 95% confidence interval)")+
  theme(axis.text.x = element_text(angle=45, hjust = 1),legend.position="bottom",plot.title = element_text(size = 20),axis.line = element_line(colour = "black", size=1.5), axis.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 20),panel.background = element_blank())



dev.off()







