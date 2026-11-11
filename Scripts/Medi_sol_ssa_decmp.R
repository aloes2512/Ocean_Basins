# from Medi_sol.locked.R line 89
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
  mutate(smth.25=predict(fit.loess))
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
  dplyr::select(dt.mnth,anomaly,"Solar_Variation"=SI,smth.25,trend.solar)%>%
  mutate(difference=trend.solar-smth.25)

#_______

N=NROW(Medi_anomaly)
# decompose a.trend.solar & b.difference
Trend.solar=Medi_anomaly$trend.solar
difference=Medi_anomaly$difference
library(Rssa)
L=floor(N/3)
s.trd<-ssa(Trend.solar,L)
s.resds<-ssa(difference,L)
# 1. Compute wcor
w.trd <- wcor(s.trd, groups = 1:20)
w.resds <-wcor(s.resds, groups = 1:20)
# 2. Extract matrix and convert to long-format tibble
w_mat.trd <- as.matrix(w.trd)
w_mat.resds <- as.matrix(w.resds)

w_df.trd <- expand.grid(F1 = 1:20, F2 = 1:20) %>%
  mutate(value = abs(as.vector(w_mat.trd)))
w_df.resds <- expand.grid(F1 = 1:20, F2 = 1:20) %>%
  mutate(value = abs(as.vector(w_mat.resds)))

# 3. Plot with clear axes
ggplot(w_df.trd, aes(x = F1, y = F2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black", limits = c(0, 1)) +
  scale_x_continuous(breaks = 1:20, expand = c(0, 0)) +
  scale_y_continuous(breaks = 1:20, expand = c(0, 0)) +
  coord_fixed() +
  labs(
    title = "W.trd-Correlation Matrix",
    x = "Eigenvector Index",
    y = "Eigenvector Index",
    fill = "W-Cor"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 8),
    panel.grid = element_blank()
  )
ggplot(w_df.resds, aes(x = F1, y = F2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue", limits = c(0, 1)) +
  scale_x_continuous(breaks = 1:20, expand = c(0, 0)) +
  scale_y_continuous(breaks = 1:20, expand = c(0, 0)) +
  coord_fixed() +
  labs(
    title = "W.resds-Correlation Matrix",
    x = "Eigenvector Index",
    y = "Eigenvector Index",
    fill = "W-Cor"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 8),
    panel.grid = element_blank()
  )
# 1. Extract the reconstructed trend component vector
rc.trd <- reconstruct(s.trd, groups = list(Trend = 1:3,
                                           Har2=4:7,
                                           Har3=7:8,
                                           Har4=9:10,
                                           trd.full=1:40))
rc.resds<- reconstruct(s.resds, groups = list( Trend.rsd=1:2,
                                               Har.rsd2=3:4,
                                               Har.rsd3=5:6,
                                               Har.rsd4=8:9,
                                               rsds.full=1:40))
Trend.rc=tibble(dt.mnth=Medi_anomaly$dt.mnth,
               Trend=unlist(rc.trd$Trend),
       Har2=unlist(rc.trd$Har2),
       Har3=unlist(rc.trd$Har3))
Trend.rc%>%
  ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Trend,col="Trend"))+
  geom_line(aes(y=Har2,col="Har2"))+
  geom_line(aes(y=Har3,col="Har3"))+
  labs(x=NULL,title="SSA decomposed Solar locked Trend",
       subtitle= "includes quasi harmonics periods (yrs):91; 41&24; 32")

Resds.rc= tibble(dt.mnth=Medi_anomaly$dt.mnth,
                 Har1=unlist(rc.resds$Trend.rsd),
                 Har2=unlist(rc.resds$Har.rsd2),
                 Har3=unlist(rc.resds$Har.rsd3),
                 Har4= unlist(rc.resds$Har.rsd4))
Resds.rc%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=Har4,col="Har4"))+
  geom_line(aes(y=Har1,col="Har1"))+
  geom_line(aes(y=Har2,col="Har2"))+
  geom_line(aes(y=Har3,col="Har3"))+
  labs(x=NULL, y= "harmonic",title="SSA decomposed Residual",
       subtitle=" periods(years): 22.67; 21.67;26.6;24.9")
ggsave("Basin_figs/resds.sol_decomposed.png")
#======= parestimate
pe.trd <- parestimate(s.trd, groups = list(c(1:3),c(4:7), c(7,8), c(9,10)))
          # Trend =91 yr;Har2= 41 yr & 24;Har3 =32 yrs;
pe.rsds= parestimate(s.resds, groups=list(c(1:2),(3:4),(5:6),(8:9)))
          # 22.67 yr ;21.67; 26.6; 24.9

