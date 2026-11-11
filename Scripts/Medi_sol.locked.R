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
library(itsmr)
M.ssn=c("season",12,"season",6)
Medi_anomaly<-Medi_anomaly%>%
  mutate(anomaly=Resid(Anomaly,M.ssn),smth.anom=predict(loess(anomaly~ dt.mnth,span = 0.15) ))
dt.mnth<-Medi_anomaly$dt.mnth
Medi_sign=Medi_anomaly$anomaly
fit.loess<-loess(anomaly ~ dt.mnth, span=0.25,data = Medi_anomaly)
Medi_anomaly<-Medi_anomaly%>%
  mutate(secular=predict(fit.loess))
# combine solar_monthly with ocean anomaly
solar_mnthly=readRDS("~/projects/Ocean_Analysis/data/S_power.rds")%>%
  mutate(SI=TSI-mean(TSI))
# anchor  at time-series solar power
library(splines)
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
Medi_anomaly=Medi_anomaly%>%
  left_join(solar_mnthly,by="dt.mnth")%>%
  dplyr::select(dt.mnth,anomaly,"Solar_Variation"=SI,secular,trend.solar)%>%
  mutate(trd_hr.res=trend.solar-secular)

plt.medi.trd=Medi_anomaly%>%ggplot(aes(x=dt.mnth))+
  geom_point(aes(y=anomaly),col="grey",size=0.3)+
  geom_line(aes(y=secular,col="trend"))+
  labs(x="",y="anomaly",title = "Mediterranean Anomaly",
       subtitle="secular = loess.smth span=0.25 ")
Medi_anomaly%>%colnames()


plt_Medi.solar<-Medi_anomaly%>%
  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=trend.solar,col="trd.sol."))+
  geom_line(aes(y=secular,col="secular"))+
  geom_line(aes(y=trd_hr.res,col="har.part"))+
  labs(x="",title="Mediterranean: Solar-locked-trend\nSecular Trend")

ggsave("Basin_figs/sol.trds.png")
saveRDS(Medi_anomaly,"data/Medi_anomaly.rds")
#--------
# 3. Extract the stable Baseline Trend and visualise
plt.decomp_solar.lckd<-Medi_anomaly%>%
  mutate(solar.res=anomaly-trend.solar)%>%
  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=solar.res),col="grey")+
  geom_line(aes(y=trend.solar,col="trend.solar"))+
  labs(title="Mediterranean-locked Trend",
       subtitle = "extracted with spline\n fit to solar power")
plt.decomp_solar.lckd
#===========
# periods of trend.solar estimate
N<-NROW(Medi_anomaly) # 2121
trd_hr.res <-Medi_anomaly$trd_hr.res
FFT.har.trd=tibble(idx=1:N-1,
                   spc=fft(trd_hr.res),
                   amp=Mod(spc))
FFT.har.trd%>%subset(idx<20)%>%
  ggplot(aes(x=idx))+
  geom_point(aes(y=amp))+
  labs(x="rel.frequency",
       title="Fourier Coef. Harmonic Trend")
FFT.max<-tibble(idx=1:N-1,spc=rep(0+0i,N),amp=Mod(spc))
which.max(FFT.har.trd$amp) #8 22yrs 530 mnth
idx.mx=FFT.har.trd%>%subset(amp>25)%>% pull(idx)
FFT.max[idx.mx,]<-FFT.har.trd[idx.mx,]
library(gsignal)
plt.har.comp=tibble(dt.mnth=Medi_anomaly$dt.mnth,
                    y.mx=Re(ifft(FFT.max$spc)))%>%
  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=y.mx))+
  labs(x=NULL,title="Harmonic Component of\n Solar Locked Trend")
library(itsmr)
# 6 max harmonics
Medi_anomaly%>%
  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=trend.solar,col="trd.solar"))+
  geom_line(aes(y=Re(ifft(FFT.max$spc)),col="har.part"))+
  geom_line(aes(y=trd_hr.res,col="trd_hr.res"))+
  labs(title="Mediterranean: Solar-locked-trend",
       subtitle = "low-pass-filtered extracted periods")
#======== Decompose SSA
library(Rssa)
medi.trd<-Medi_anomaly$trend.solar
L = floor(N/3)
Medi.s=ssa(medi.trd,L=L)
# 1. Compute wcor
w <- wcor(s, groups = 1:20)
# 2. Extract matrix and convert to long-format tibble
w_mat <- as.matrix(w)

w_df <- expand.grid(F1 = 1:20, F2 = 1:20) %>%
  mutate(value = abs(as.vector(w_mat)))

# 3. Plot with clear axes
ggplot(w_df, aes(x = F1, y = F2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black", limits = c(0, 1)) +
  scale_x_continuous(breaks = 1:20, expand = c(0, 0)) +
  scale_y_continuous(breaks = 1:20, expand = c(0, 0)) +
  coord_fixed() +
  labs(
    title = "W-Correlation Matrix",
    x = "Eigenvector Index",
    y = "Eigenvector Index",
    fill = "W-Cor"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 8),
    panel.grid = element_blank()
  )
plot(s,type="values")
