library(terra)
library(tidyverse)
library(itsmr)
library(zoo)
library(Rssa)
library(splines)
library(minpack.lm)
file.create(".here")
import_global.ocean <- function() {
  # Load historical dataset from external project
  Ocean_data_update <- readRDS("~/projects/Ocean_Analysis/data/Ocean_data_update.rds")
  ocean_stack <- Ocean_data_update$sst_stack
  sst_stack <- terra::unwrap(ocean_stack)
  date <- time(sst_stack)

  # Load Berkeley Earth mask and derive ocean fractions/masks
  berk_mask <- terra::rast("~/projects/Ocean_Analysis/data/Land_and_Ocean_1x1_SST_at_sea_ice.nc")$land_mask

  ocean_mask <- terra::mask(berk_mask, berk_mask < 0.1, maskvalues = FALSE)
  ocean_frac <- 1 - berk_mask

  # Calculate spatial area weights (relative to maximum equatorial cell area)
  weights_abs <- terra::cellSize(sst_stack[[1]], unit = "km")
  weights_rel <- weights_abs / terra::global(weights_abs, "max", na.rm = TRUE)[1, 1]

  # Resample and apply NA masking to match target geometry
  ocean_mask_fixed <- terra::resample(ocean_frac, sst_stack, method = "near")
  ocean_mask_na <- ocean_mask_fixed
  ocean_mask_na[ocean_mask_na == 0] <- NA

  # Return core objects as a list for downstream analysis
  return(list(
    sst_stack     = sst_stack,
    date          = date,
    weights_rel   = weights_rel,
    ocean_mask_na = ocean_mask_na,
    ocean_mask_fixed = ocean_mask_fixed
  ))
}
region_ts<-function(regext){
  global.data<-import_global.ocean()
  sst_stack<-global.data$sst_stack
  ocean_mask_na<-global.data$ocean_mask_na
  weights_rel<-global.data$weights_rel
  extd=regext
  # crop data to regext
  reg_stack  <- crop(sst_stack, extd)
  reg_mask  <- crop(ocean_mask_na, extd)
  reg_weights_rel <- crop(weights_rel, extd)
  # correct coord
  reg_stack        <- rotate(reg_stack)
  reg_mask    <- rotate(reg_mask)
  reg_weights_rel      <- rotate(reg_weights_rel)
  # 3. Combine mask and weights for the region
  reg_effective_weights <- reg_weights_rel * reg_mask
  # 4. Find the total sum of weights specifically within the Mediterranean basin
  reg_total_weight <- global(reg_effective_weights, fun = "sum", na.rm = TRUE)[1, 1]
  regional_ts <- global(reg_stack * reg_effective_weights, fun = "sum", na.rm = TRUE) / reg_total_weight
  regional_ts<-unlist(regional_ts)
  ts=tibble(dt.mnth=as.yearmon(time(reg_stack))%>%as.numeric(),
            Anomaly=regional_ts)

  return(ts)
}



global.data<-import_global.ocean()
med_extent <- ext(-6, 36, 30, 46)
Medi_anomaly=region_ts(med_extent)
M.ssn=c("season",12,"season",6)
Medi_anomaly<-Medi_anomaly%>%
  mutate(anomaly=Resid(Anomaly,M.ssn),smth.anom=predict(loess(anomaly~ dt.mnth,span = 0.15) ))
# combine solar_monthly with ocean anomaly
solar_mnthly=readRDS("~/projects/Ocean_Analysis/data/S_power.rds")%>%
  mutate(SI=TSI-mean(TSI))
# anchor  at time-series solar power

# find zero crossings of Solar Irradiation SI
find_nodes <- function(x) {
  which(diff(sign(diff(x))) == 2) + 1
}
# nodes <- find_nodes(solar_mnthly$SI)
# time coordinates of zero-crossings
zero_crossing_times <- solar_mnthly$dt.mnth[find_nodes(solar_mnthly$SI)]

# 2. Fit a continuous spline locked @ exact nodes
# degree = 1 creates  sawtooth a continuous line that bends at the zero-crossings
# spline: degree = 3 creates a smooth, continuous curve through the zero-crossings
trd.Medisolar <- lm(anomaly ~ bs(dt.mnth, knots = zero_crossing_times, degree = 3),
                    data = Medi_anomaly)

Medi_anomaly<-Medi_anomaly%>%mutate(trend.solar=predict(trd.Medisolar))
Medi_residual<-Medi_anomaly%>%mutate(residual=anomaly-trend.solar)

#========
Medi_anomaly%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=anomaly-trend.solar),col="grey")+
  geom_line(aes(y=trend.solar))
# Spectrum Isolation on a. Trend and b. Residuals
 solar.trend<-Medi_anomaly$trend.solar
 solar.residuals <-Medi_residual$residual
N<-length(solar.residuals) #2121
L<-floor(N/3) # 707
s.trend<-ssa(solar.trend,L)
s.resds<-ssa(solar.residuals,L)
#=====
w.trend<-wcor(s.trend, groups = 1:20)%>%as.matrix()
w.resds<-wcor(s.resds,groups=1:20)%>%as.matrix()
w_df <- expand.grid(F1 = 1:20, F2 = 1:20) %>%
  mutate(value.trend = abs(as.vector(w.trend)),
         value.resds = abs(as.vector(w.resds)))
# plot the two cor matrices
ggplot(w_df, aes(x = F1, y = F2, fill = value.trend)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black", limits = c(0, 1)) +
  scale_x_continuous(breaks = 1:20, expand = c(0, 0)) +
  scale_y_continuous(breaks = 1:20, expand = c(0, 0)) +
  coord_fixed()+
  labs(title = "trend cor matrix")
ggplot(w_df, aes(x = F1, y = F2, fill = value.resds)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue", limits = c(0, 1)) +
  scale_x_continuous(breaks = 1:20, expand = c(0, 0)) +
  scale_y_continuous(breaks = 1:20, expand = c(0, 0)) +
  coord_fixed()+
  labs(title = "resds cor matrix")
# reconstruct sum dominant periods
rc.Medi.trend <- reconstruct(s.trend, groups = list(domin=1:2,smth.trd=1:10,full.trd=1:20,harmon=3:10)
)
tibble(dt.mnth=Medi_anomaly$dt.mnth,
       smth.trd=unlist(rc.Medi.trend$smth.trd),
       full.trd=unlist(rc.Medi.trend$full.trd),
       domin=unlist(rc.Medi.trend$domin),
       harmon=unlist(rc.Medi.trend$harmon))%>%
  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=full.trd),lwd=2,col="grey")+
  geom_line(aes(y=smth.trd),col=2)+
  geom_line(aes(y=domin),linetype=2)+
  geom_line(aes(y=harmon))
tibble(dt.mnth=Medi_anomaly$dt.mnth,
       smth.trd=unlist(rc.Medi.trend$smth.trd),
       full.trd=unlist(rc.Medi.trend$full.trd),
       domin=unlist(rc.Medi.trend$domin),
       harmon=unlist(rc.Medi.trend$harmon))%>%
  pivot_longer(cols=-dt.mnth,names_to = "Componts",
               values_to= "trends")%>%
  ggplot(aes(x=dt.mnth))+geom_line(aes(y=trends,col=Componts))
# parestimate
pe.trd <- parestimate(s.trend, groups = list(c(1:2),c(3:10), c(1:10)))
3450/12 # 287 yrs
252/12 # 21 yrs
586/12 # 49 yrs
370/12 # 30.8 yrs
1577/12 # 131 yrs
# reconstruct residuals
rc.Medi.resds <- reconstruct(s.resds, groups = list(domin=1:2,smth.trd=1:6,full.trd=1:10,harmon=3:10)
)
tibble(dt.mnth=Medi_anomaly$dt.mnth,
       smth.trd=unlist(rc.Medi.resds$smth.trd),
       full.trd=unlist(rc.Medi.resds$full.trd),
       domin=unlist(rc.Medi.resds$domin),
       harmon=unlist(rc.Medi.resds$harmon))%>%
  pivot_longer(cols=-dt.mnth,names_to = "Componts",
               values_to= "resds")%>%
  ggplot(aes(x=dt.mnth))+geom_line(aes(y=resds,col=Componts))
# rollingcor to show str break on dominant harmonic
library(slider)
Medi_residual<-Medi_residual%>%mutate(lng.anoma=hr(anomaly,N/1:3),lng.resd=hr(residual,N/1:3))
vector_obs<-Medi_residual$lng.anoma
vector_res<-Medi_residual$lng.resd
rolling_corr.Medi <- slide_index_dbl(
  .x = 1:length(vector_obs),
  .i = 1:length(vector_obs),
  .f = function(idx) {
    # If the window doesn't have a full 11 years of data, return NA
    if(length(idx) < 132) return(NA)

    # Calculate simple vector correlation on the current slice of indices
    cor(vector_obs[idx], vector_res[idx])
  },
  .before = 132,   # Looks 131 steps backward + current step = 132 months
  .complete = TRUE # Automatically returns NA for incomplete windows at the start
)
tibble(dt.mnth=Medi_anomaly$dt.mnth,
       Rcorr=rolling_corr.Medi)%>%
  ggplot(aes(x=dt.mnth,y=Rcorr))+geom_line()

library(mgcv)
library(broom)
N<-NROW(Medi_anomaly)
N_train<- floor(N/2)
df_full<-Medi_anomaly
df_train<-Medi_anomaly[1:N_train,]
gam_train <- gam(anomaly ~ s(dt.mnth, k = 90, bs = "cc"), data = df_train)
df_train<-augment(gam_train)%>%dplyr::select(dt.mnth,"train"=.fitted)
gam_full<- gam(anomaly~ s(dt.mnth,k=180,bs="cc"),data=Medi_anomaly)
df_train<-df_full%>%mutate(ext.preind=predict(gam_train, newdata = df_full))
df_lng.anoma<-gam_full%>%augment()%>%
  dplyr::select(dt.mnth,"full"=.fitted)
df_full<-df_full%>%left_join(df_lng.anoma)%>%left_join(df_train)
mean.histor=mean(df_full$ext.preind)
df_full%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=anomaly),col="grey")+
  geom_line(aes(y=ext.preind),col=2)+
  geom_line(aes(y=full),col=3)+
  geom_hline(yintercept = mean.histor,col=4)

library(slider)
vector_obs<-df_full$full
vector_res<-df_full$ext.preind
rolling_gam.Medi <- slide_index_dbl(
  .x = 1:length(vector_obs),
  .i = 1:length(vector_obs),
  .f = function(idx) {
    # If the window doesn't have a full 11 years of data, return NA
    if(length(idx) < 144) return(NA)

    # Calculate simple vector correlation on the current slice of indices
    cor(vector_obs[idx], vector_res[idx])
  },
  .before = 144,   # Looks 131 steps backward + current step = 132 months
  .complete = TRUE # Automatically returns NA for incomplete windows at the start
)
length(rolling_gam.Medi)
Roll.gam.tbl=tibble(dt.mnth=Medi_anomaly$dt.mnth,
       Rcor=rolling_gam.Medi)
Roll.gam.tbl%>%
  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Rcor))
Medi_anomaly$dt.mnth[which.min(rolling_gam.Medi)] # 1957.5
T.cor<-Roll.gam.tbl%>%arrange(Rcor)%>%subset(dt.mnth>1960)%>%slice(1)
TP=T.cor[[1]]
Medi_anomaly%>%subset(dt.mnth==TP)


fit_logistic <- nlsLM(
  anomaly ~ b+A*(1+exp(-B*(dt.mnth-TP)))^-(1/nu),
  data = df_full,
  start=list(b=0,A=2,B=0.12,nu=1),
  control = nls.control(maxiter = 200)
)
coef(fit_logistic) # b= -0.16; A= 2.9; B= 0.064; nu= 0.35
summary(fit_logistic)
Logst.tbl=tibble(dt.mnth=Medi_anomaly$dt.mnth,
       logstc=predict(fit_logistic))
df_full%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=anomaly),col="grey")+
  geom_line(aes(y=ext.preind),col=2)+
  geom_line(aes(y=full),col=3)+
  geom_hline(yintercept = mean.histor,col=4)+
  geom_vline(xintercept = TP,linetype = 3)+
  geom_line(data=Logst.tbl,aes(y=logstc))+
  labs(x=NULL,title = "Mediterranean Decomposed",
       subtitle="Richards trend symmetry point: Rollcor ")
# separating the human caused trend from observed
df_full<-df_full%>% left_join(Logst.tbl)
DF_natrl<-df_full%>%dplyr::select(dt.mnth,anomaly,logstc)%>%
  mutate(natr.resds=anomaly-logstc)
DF_natrl%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=natr.resds),col="grey")+
  geom_smooth(aes(y=natr.resds))+
  labs(x=NULL,title = "Natural Harmonics Part of Mediterranean\n Anomaly")
ggsave("Basin_figs/Natural_anomaly.png")
