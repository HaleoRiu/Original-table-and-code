# ===========================
# Option 1: SAR residualization -> PLS-PM
# ===========================

rm(list = ls())
setwd("your/path")
# loading packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(spdep)
  library(spatialreg)
  library(plspm)
})

# 1) Reading data --------------------------------------------------------------
dat <- read.table("env.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
stopifnot(all(req_cols %in% names(dat)))

# 2) Creat martix----------------------------------------
coords <- as.matrix(dat[, c("long","lat")])

k <- 4
nb  <- knn2nb(knearneigh(coords, k = k))

lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

# 3) Data list ---------------------------------------
vars_to_resid <- c("HeS_B","DL_B","DR_B","DL_F","DR_F","AN","AP","TK",
                   "ASV3502", "ASV2170", "ASV1216", "ASV1029", 
                   "ASV5069", "ASV30", "ASV108", "ASV26", 
                   "ASV3")

# 4) SAR -----------------------------
sar_residualize <- function(yname, dat, lw) {
  frm <- as.formula(paste(yname, "~ elev + HFI")) 
  fit <- lagsarlm(frm, data = dat, listw = lw, zero.policy = TRUE, method = "eigen")
  res <- residuals(fit) 
  mi  <- moran.test(res, lw, zero.policy = TRUE)
  cat("\n---", yname, "SAR model ---\n")
  print(summary(fit))
  cat("Moran's I of residuals:\n")
  print(mi)
  return(res)
}
for(v in vars_to_resid) {
  dat[[paste0(v, "_res")]] <- sar_residualize(v, dat, lw)
}


# 5) blocks & paths  -----------------------------------------

dat$DR_B_res = -1 * dat$DR_B_res
dat$DR_F_res = -1 * dat$DR_F_res
dat$AN_res = -1 * dat$AN_res
dat$DR_B = -1 * dat$DR_B
dat$DR_F = -1 * dat$DR_F

fix_blocks <- list(
  Elev      = c("elev"),
  HFI     = c("HFI"),
  Keystone_Bac = c("ASV2170", "ASV5069", "ASV26", "ASV3"),
  AsseP_F = c("DL_F", "DR_F"),
  AsseP_B    = c("DL_B", "DR_B")
  
)

fix_modes <- c("A","A","B","B","B")


Elev      <- c(0,0,0,0,0)
HFI     <- c(1,0,0,0,0)
Keystone_Bac <- c(1,1,0,0,0)
AsseP_F <- c(1,1,1,0,0)
AsseP_B <- c(1,1,1,1,0)
fix_path <- rbind(Elev, HFI, Keystone_Bac, AsseP_F, AsseP_B)
colnames(fix_path) <- rownames(fix_path) <- c("Elev","HFI","Keystone_Bac", 
                                              "AsseP_F", "AsseP_B")

set.seed(123)
pls_fit <- plspm(dat, fix_path, fix_blocks, modes = fix_modes)
print(summary(pls_fit))
cat("\nGoodness-of-Fit:\n")
print(pls_fit$gof)

try({
  innerplot(pls_fit, colpos = 'red3', colneg = 'navy', show.values = TRUE,
            lcol = 'grey30', box.lwd = 1.2)
}, silent = TRUE)
