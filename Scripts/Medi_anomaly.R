library(terra)
library(tidyverse)
library(itsmr)
library(zoo)
library(Rssa)
library(splines)
library(minpack.lm)
library(knitr)
library(kableExtra)
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
# Solar locked trend
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
Medi_anomaly<-Medi_anomaly%>%mutate(residual=anomaly-trend.solar)
Medi_anomaly%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=trend.solar),col=2,lwd= 1.3)+
  geom_point(aes(y=residual),size=0.2,col="grey")+
  labs(x=NULL,title="Solar Locked Medterranean Anomaly",
       subtitle="decomposed: trend harmonics ")
ggsave("Basin_figs/Medi.anoma.decomp.png")
#===================

# SSA decompose
#medi_residual<-Medi_anomaly$residual
Medi_anomaly<-Medi_anomaly%>% mutate(smth.res=smooth.fft(residual,f=0.1))
medi.res.smth<-Medi_anomaly$smth.res
N<- NROW(Medi_anomaly)
L <- floor(N/3)
s.res <- ssa(medi.res.smth, L = L)
# 1. Compute wcor
w.res <- wcor(s.res, groups = 1:20)
# 2. Extract matrix and convert to long-format tibble
w.res_mat <- as.matrix(w.res)

w.res_df <- expand.grid(F1 = 1:20, F2 = 1:20) %>%
  mutate(value = abs(as.vector(w.res_mat)))
# 3. Plot with clear axes
ggplot(w.res_df, aes(x = F1, y = F2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black", limits = c(0, 1)) +
  scale_x_continuous(breaks = 1:20, expand = c(0, 0)) +
  scale_y_continuous(breaks = 1:20, expand = c(0, 0)) +
  coord_fixed() +
  labs(
    title = "W.res-Correlation Matrix",
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
rc.Medi.res <- reconstruct(s.res, groups = list(base =1:2,Har2=3:4,Har3=5:6,Har4=7:10,full=1:20)
                      )

# Estimate parameters (frequencies, periods, amplitudes) for component groups
pe <- parestimate(s.res, groups = list(c(1:2),c(3,4), c(5, 6), c(7, 10),c(1:20)))


Medi.res_rec=tibble(dt.mnth=Medi_anomaly$dt.mnth,
                  base=unlist(rc.Medi.res$base),   # 12.1 /12
                  Har2=unlist(rc.Medi.res$Har2),  #  3.9 yrs
                  Har3=unlist(rc.Medi.res$Har3),  #  5.8 yrs
                  Har4=unlist(rc.Medi.res$Har4),   #  4.6 yrs
                  full=unlist(rc.Medi.res$full))
Medi.res_rec%>%
  pivot_longer(cols=-dt.mnth,
               names_to = "comp",
               values_to = "resid")%>%
  ggplot(aes(x=dt.mnth,y=resid,col=comp))+
  geom_line()+labs(x=NULL,title="Mediterranean Residuals Decomposed")
ggsave("Basin_figs/SSA_dominant_residuals.png")
# harmonics
Medi.res_rec%>%ggplot(aes(x=dt.mnth))+
  geom_line(aes(y=base),col=2)+
  geom_line(aes(y=Har2),col=3)+
  geom_line(aes(y=Har3),col=4)
Medi.res_rec%>%dplyr::select(-c(full,Har4))%>%
  pivot_longer(cols=-dt.mnth,names_to = "compnts",values_to = "resids")%>%
  ggplot(aes(x=dt.mnth,y=resids,col=compnts))+
  geom_line()+
  labs(x=NULL,title="Medi Residual SSA decomposed",
       subtitle = " dominant perriods (yrs): base= 4; Har2= 6.6; Har3 = 1.8 ")
ggsave("Basin_figs/Medi_sol.residuals.ssa.png")
pe.res <- parestimate(s.res, groups = list(base=c(1,2),Har2=c(3, 4),Har3= c(5, 6),Har4= c(7, 8)))

pe.res_tibble <- map_dfr(names(pe.res), function(nm) {
  item <- pe.res[[nm]]

  # Convert the fdimpars object to a data frame by extracting its internal structure or unclassing
  # Usually fdimpars objects contain a matrix or data frame of the estimates
  df<- tibble()

  df <- as.data.frame(unclass(item))
  df$component <- nm
  df<-df%>%dplyr::select(component,periods,moduli,rates)

  df
}, .id = NULL) %>% as_tibble()%>%subset(periods>0)
pe.res_tibble
#======================
trd.signal<-Medi_anomaly$trend.solar
N<- length(trd.signal) #2121
L <- floor(N/3)
s.trd <- ssa(trd.signal, L = L)
# 1. Compute wcor
w.trd <- wcor(s.trd, groups = 1:20)
# 2. Extract matrix and convert to long-format tibble
w.trd_mat <- as.matrix(w.trd)

w.trd_df <- expand.grid(F1 = 1:20, F2 = 1:20) %>%
  mutate(trd.value = abs(as.vector(w.trd_mat)))
# 3. Plot with clear axes
ggplot(w.trd_df, aes(x = F1, y = F2, fill = trd.value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black", limits = c(0, 1)) +
  scale_x_continuous(breaks = 1:20, expand = c(0, 0)) +
  scale_y_continuous(breaks = 1:20, expand = c(0, 0)) +
  coord_fixed() +
  labs(
    title = "W_trend-Correlation Matrix",
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
rc.Medi.trd <- reconstruct(s.trd, groups = list(base =1:2,Har2=3:4,Har3=5:6,Har4=7:10,full=1:20)
)
pe.trd<-parestimate(s.trd, groups = list(base=c(1,2),Har2=c(3, 4),Har3= c(5, 6),Har4= c(7, 8)))

pe.trd_tibble <- map_dfr(names(pe.trd), function(nm) {
  item <- pe.trd[[nm]]

  # Convert the fdimpars object to a data frame by extracting its internal structure or unclassing
  # Usually fdimpars objects contain a matrix or data frame of the estimates
  df<- tibble()

  df <- as.data.frame(unclass(item))
  df$component <- nm
  df<-df%>%dplyr::select(component,periods,moduli,rates)

  df
}, .id = NULL) %>%
  as_tibble()%>%subset(periods>0)

Medi.trd_rec=tibble(dt.mnth=Medi_anomaly$dt.mnth,
                base=unlist(rc.Medi.trd$base),   # 1.01=12.1 /12
                Har2=unlist(rc.Medi.trd$Har2),  #  3.9 yrs
                Har3=unlist(rc.Medi.trd$Har3),  #  5.8 yrs
                Har4=unlist(rc.Medi.trd$Har4),   #  0.9 yrs
                full=unlist(rc.Medi.trd$full))
pe.trd<-parestimate(s.trd, groups = list(base=c(1,2),Har2=c(3, 4),Har3= c(5, 6),Har4= c(7, 8)))
trd.har=pe.trd_tibble%>%subset(component %in% c("Har2","Har3"))
trd.har$component<-c("trd.Har2","trd.Har3")
plt.sol.trd.har=Medi.trd_rec%>%dplyr::select(dt.mnth,base,Har2,Har3)%>%
  pivot_longer(cols=-dt.mnth,names_to = "components",values_to = "trd.sol")%>%
  ggplot(aes(x=dt.mnth,y=trd.sol,col=components))+
  geom_line()+
  labs(x=NULL,title="Dominant Trend (Long Periods)\nMediterranean Solar Locked ")

ggsave("Basin_figs/Sol.trd.harmon.pdf")

natr.harmonics=pe_tibble %>%
  mutate(
    periods = round(periods, 1),
    moduli  = round(moduli, 3),
    rates   = round(rates, 5)
  ) %>% bind_rows(trd.har)%>%
  kable(
    format = "markdown",
    col.names = c("Component", "Period (mo)", "Modulus", "Rate"),
    caption = "Estimated parameters for SSA oscillatory components."
  )

natr.harmonics
