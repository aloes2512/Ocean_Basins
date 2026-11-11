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
  reg_regional_ts <- global(reg_stack * reg_effective_weights, fun = "sum", na.rm = TRUE) / reg_total_weight
  return(reg_regional_ts)
}