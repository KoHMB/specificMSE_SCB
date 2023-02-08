###############################################################################
###############################################################################
#
# OM作成のためのコンディショニングコード
#
###############################################################################
###############################################################################
library(knitr)
library(tidyverse)
library(TMB)

###############################################################################
##---------------------------------------------------------------------------##
##        Conditionig
##---------------------------------------------------------------------------##

dat <- read.csv("SC_Jap-B4.csv")
TT <- length(dat$year[21:71])
data <- list(y1=dat$CPUE_nominal[21:71], y2=dat$survey_B[21:71], Catch=dat$catch_B[21:71], eps=0.001)

#compile("sspm2.cpp")
dyn.load(dynlib("sspm2"))

par_init <- list(
  logit_D0 = 1,
  log_r = log(0.3),
  log_K = log(10000),
  log_z = log(1),
  log_q1 = log(0.0008),
  log_q2 = log(2.),
  log_sig1 = log(0.2),
  log_sig2 = log(0.3),
  log_tau = log(0.2),
  log_Dep = c(log(readRDS("Dep_int2.rds")[21:70]), log(0.8))
)
#map <- list(log_z = factor(NA))
map <- list(log_q2 = factor(NA), log_z = factor(NA))

# Estimation by TMB
obj.sspm <- MakeADFun(data, par_init, random="log_Dep", map=map, DLL="sspm2")
opt<-nlminb(obj.sspm$par, obj.sspm$fn, obj.sspm$gr)#, 
            #lower = obj.sspm$par/10#,
            #upper = obj.sspm$par*20
            #)
obj <- obj.sspm
LL_F <- opt$objective

# Extracting results
sdr.sspm <- sdreport(obj.sspm)
res1_est <- sdr.sspm$value
res1_sd <- sdr.sspm$sd
#prof.L <- tmbprofile(obj.sspm, "log_q2", ytol=5, parm.range=c(-Inf,log(1.25)))
#plot(exp(prof.L), xlim=c(0,1.25), cex=1.5, lwd = 3)
points(exp(prof.L[which(min(prof.L$value,na.rm=T)==prof.L$value)[1],]),
       pch=18, col="blue", cex=5)

r_est <- res1_est[1] %>% as.numeric()
K_est <- res1_est[2] %>% as.numeric()
MSY_est <- res1_est[3] %>% as.numeric()
Bmsy_est <- res1_est[4] %>% as.numeric()
q1_est <- res1_est[105] %>% as.numeric()
q2_est <- res1_est[106] %>% as.numeric()
sig1_est <- res1_est[107] %>% as.numeric()
sig2_est <- res1_est[108] %>% as.numeric()
tau_est <- res1_est[109] %>% as.numeric()

B.est <- sdr.sspm$value[names(sdr.sspm$value)=="B"]
B.se <- sdr.sspm$sd[names(sdr.sspm$value)=="B"]
B.LB <- B.est*exp(-1.96*B.se/B.est)
B.UB <- B.est*exp(+1.96*B.se/B.est)

# Plotting results
Year <- dat$year[21:70]
plot(Year, B.est, type="l", col="red", ylab="Biomass (tons)", 
     ylim=c(0, max(B.UB,na.rm = TRUE)), main = "Dynamics (1970年-)")
     #ylim=c(0, 50000), main = "Dynamics (1970年-)", lwd=2)
#points(Year, B.LB, type="l", lty=2, col="red")
#points(Year, B.UB, type="l", lty=2, col="red")
polygon(c(Year,rev(Year)), c(B.LB,rev(B.UB)), 
        border = gray(.5,alpha=.5), col = gray(.5,alpha=.2))
#points(Year, data$Catch, type = "h", lwd=5, col="darkgray", ylab="Tons")
abline(h=K_est, col="blue", lty=2)

summary_res <- data.frame(Model = "PT(z=0.5)",
                          #Model = "Schaefer",
                          r = r_est,
                          K = K_est,
                          MSY = MSY_est,
                          Bmsy = Bmsy_est, 
                          q1 = q1_est,
                          q2 = q2_est,
                          obs1_error = sig1_est,
                          obs2_error = sig2_est,
                          process_error = tau_est
)
knitr::kable(summary_res)
saveRDS(summary_res, "res_SSPM-PT100_TMB.RDS")
#saveRDS(summary_res, "res_SSPM-S_TMB.RDS")

curve(-r_est/K_est*x^2+r_est*x, xlim = c(0,K_est), cex=2, lwd=2)
abline(h=250, col="red", lwd = 1.5)
points(last(B.est),-r_est/K_est*last(B.est)^2+r_est*last(B.est), col="blue", cex=2.5, pch=16)
segments(last(B.est),-r_est/K_est*last(B.est)^2+r_est*last(B.est),last(B.est),0, lty=2, col="blue")




################################################################################
# Grid approach

z_sets <- c(0.1,0.5,1)
q2_sets <- c(0.35, 1, 1.5)
par_sets <- expand.grid(z_sets, q2_sets)
colnames(par_sets) <- c("z", "q2")

NLL <- numeric(nrow(par_sets))
for(i in 1:nrow(par_sets)){
  par_init$log_z <- log(par_sets[i,1])
  par_init$log_q2 <- log(par_sets[i,2])
  
  obj.sspm <- MakeADFun(data, par_init, random="log_Dep", map=map, DLL="sspm2")
  opt<-nlminb(obj.sspm$par, obj.sspm$fn, obj.sspm$gr)
  obj <- obj.sspm
  NLL[i] <- opt$objective
  sdr.sspm <- sdreport(obj.sspm)
  eval(parse(text = paste0('saveRDS(sdr.sspm, 
                           "Res_conditioning/20210803/z-', par_sets[i,1],
                           '_q2-',par_sets[i,2],
                           '.rds")')))
  rm(obj.sspm); rm(sdr.sspm) 
}
saveRDS(cbind.data.frame(par_sets, NLL = NLL),
        "Res_conditioning/20210803/OMs_Likelihood.rds")

# robustness set 1
z_sets <- c(0.1,0.5,1)
q2_sets <- c(2)
par_sets <- expand.grid(z_sets, q2_sets)
colnames(par_sets) <- c("z", "q2")
NLL <- numeric(nrow(par_sets))
for(i in 1:nrow(par_sets)){
  par_init$log_z <- log(par_sets[i,1])
  par_init$log_q2 <- log(par_sets[i,2])
  obj.sspm <- MakeADFun(data, par_init, random="log_Dep", map=map, DLL="sspm2")
  opt<-nlminb(obj.sspm$par, obj.sspm$fn, obj.sspm$gr)
  obj <- obj.sspm
  NLL[i] <- opt$objective
  sdr.sspm <- sdreport(obj.sspm)
  eval(parse(text = paste0('saveRDS(sdr.sspm, 
                           "Res_conditioning/20210803/z-', par_sets[i,1],
                           '_q2-',par_sets[i,2],
                           '.rds")')))
  rm(obj.sspm); rm(sdr.sspm) 
}
saveRDS(cbind.data.frame(par_sets, NLL = NLL),
        "Res_conditioning/20210803/OMs_Likelihood.rds")











max(data$Catch, na.rm = TRUE)
max(data$y2, na.rm = TRUE)


# profile likelihood
m_sets <- seq(0.01, 1.5, length=100)+1
q2_sets <- (seq(0.1, 2, length=100))
par_sets <- expand.grid(z_sets, q2_sets)
#invq2_sets <- ((seq(1/0.2, 1/10, length=100))) %>% sort(decreasing = FALSE)
#par_sets <- expand.grid(z_sets, invq2_sets)
colnames(par_sets) <- c("m", "q2")

NLL <- numeric(nrow(par_sets))
for(i in 1:nrow(par_sets)){
  par_init$log_z <- log(par_sets[i,1]-1)
  par_init$log_q2 <- log(par_sets[i,2])
#  par_init$log_q2 <- log(1/par_sets[i,2])
  
  obj.sspm <- MakeADFun(data, par_init, random="log_Dep", map=map, DLL="sspm2", silent = TRUE)
  opt<-nlminb(obj.sspm$par, obj.sspm$fn, obj.sspm$gr) %>% try(silent = TRUE)
  if(class(opt)=="try-error") {
    NLL[i] <- NA
  } else if(!opt$convergence==0) {
      NLL[i] <- NA
  } else {
    obj <- obj.sspm
    NLL[i] <- opt$objective
    #sdr.sspm <- sdreport(obj.sspm)
  }
  rm(opt); rm(obj.sspm)#; rm(sdr.sspm) 
}
#saveRDS(cbind.data.frame(par_sets, NLL = NLL), "Res_conditioning/R3-SA/OMs_Likelihood.rds")
NLL <- readRDS("Res_conditioning/R3-SA/OMs_Likelihood.rds")

## q2のスケール
data1 <- readRDS("Res_conditioning/R3-SA/OMs_Likelihood.rds")
data1$NLL <- data1$NLL - min(data1$NLL, na.rm = TRUE)
data1$deltaNLL <- rep(NA, nrow(data1))
for(i in 1:nrow(data1)){
  if(is.na(data1$NLL[i])){
    data1$deltaNLL[i] <- NA
  } else if(data1$NLL[i]<0.1){
    data1$deltaNLL[i] <- "0-"
  } else if(data1$NLL[i]<0.2){
    data1$deltaNLL[i] <- "0.1-"
  } else if(data1$NLL[i]<0.5){
    data1$deltaNLL[i] <- "0.2-"
  } else if(data1$NLL[i]<2){
    data1$deltaNLL[i] <- "0.5-"
  } else if(data1$NLL[i]<4){
    data1$deltaNLL[i] <- "2-"
  } else if(data1$NLL[i]<10){
    data1$deltaNLL[i] <- "4-"
  } else {
    data1$deltaNLL[i] <- "10-"
  } 
}
g_profile <- ggplot(data = data1) +
  geom_tile(aes(x=m, y=q2, fill=deltaNLL))+
  geom_point(data = data.frame(x = rep(c(1.1,1.5,2),3), 
                               y = rep(c(0.35,1,1.5), each = 3)),
             aes(x=x, y=y), col ="red", size=5)+
  geom_point(data = data.frame(x = data1[data1$NLL==min(data1$NLL),1], 
                               y = data1[data1$NLL==min(data1$NLL),2]),
             aes(x=x, y=y), col = "red", size = 10, shape = 4)+
  scale_fill_manual(values = c("gray10", "gray15", "gray20", "gray25",
                               "gray35", "gray70", "gray90"),
                    na.value = "white")+
  coord_cartesian(xlim = c(1,2.5), ylim = c(0,2), expand = 0)+
  ylab("q2")+
  theme_classic()+
  theme(text = element_text(size = 20))


## 逆数のスケール
data2 <- data.frame(z=par_sets[,1], q2=par_sets[,2], NLL=NLL-min(NLL, na.rm = TRUE))
data2$deltaNLL <- rep(NA, nrow(data2))
for(i in 1:nrow(data2)){
  if(is.na(data2$NLL[i])){
    data2$deltaNLL[i] <- NA
  } else if(data2$NLL[i]<0.1){
    data2$deltaNLL[i] <- "0-"
  } else if(data2$NLL[i]<0.2){
    data2$deltaNLL[i] <- "0.1-"
  } else if(data2$NLL[i]<0.5){
    data2$deltaNLL[i] <- "0.2-"
  } else if(data2$NLL[i]<2){
    data2$deltaNLL[i] <- "0.5-"
  } else if(data2$NLL[i]<4){
    data2$deltaNLL[i] <- "2-"
  } else if(data2$NLL[i]<10){
    data2$deltaNLL[i] <- "4-"
  } else {
    data2$deltaNLL[i] <- "10-"
  } 
}
data2$deltaNLL <- factor(data2$deltaNLL, levels=rev(unique(data2$deltaNLL)))
g_profile <- 
  ggplot(data = data2) +
  geom_tile(aes(x=z, y=q2, fill=deltaNLL))+
  #geom_contour(aes(x=z, y=q2, z=deltaNLL, color=..level..), breaks = 0.5, color=1)+
  #geom_contour(aes(x=z, y=q2, z=deltaNLL, color=..level..), breaks = 2, color=1)+
  #geom_contour(aes(x=z, y=q2, z=deltaNLL, color=..level..), breaks = 4, color=1)+
  geom_point(data = data.frame(x = rep(c(0.1,0.5,1),3), 
                               y = rep(1/c(0.35,1,1.5), each = 3)),
             aes(x=x, y=y), col ="red", size=5)+
  geom_point(data = data.frame(x = data2[data2$NLL==min(data2$NLL),1], 
                               y = data2[data2$NLL==min(data2$NLL),2]),
             aes(x=x, y=y), col = "red", size = 10, shape = 4)+
  #scale_fill_grey(na.value = "white")+
  scale_fill_manual(values = c("gray10", "gray15", "gray20", "gray25",
                               "gray35", "gray70", "gray90"),
                    na.value = "white")+
  #scale_x_continuous(limits = c(0,2))+
  #scale_y_continuous(limits = c(0,5))+
  coord_cartesian(xlim = c(0,1.5), ylim = c(0,5), expand = 0)+
  #xlim(0,2)+ylim(0,5)+
  ylab("1/q2")+
  theme_classic()+
  theme(text = element_text(size = 20))

saveRDS(g_profile, "Res_Graph/profile_z-q2.rds")







# sensitivity case
z_sets <- c(0.01, 0.1,0.5,1)
NLL <- numeric(length(z_sets))
for(i in 1:length(z_sets)){
  par_init$log_z <- log(z_sets[i])
  par_init$log_q2 <- log(2)
  
  obj.sspm <- MakeADFun(data, par_init, random="log_Dep", map=map, DLL="sspm2")
  opt<-nlminb(obj.sspm$par, obj.sspm$fn, obj.sspm$gr)
  obj <- obj.sspm
  NLL[i] <- opt$objective
  sdr.sspm <- sdreport(obj.sspm)
  eval(parse(text = paste0('saveRDS(sdr.sspm, 
                           "Res_conditioning/z-', z_sets[i],
                           '_q2-2.rds")')))
  rm(obj.sspm); rm(sdr.sspm) 
}


# Likelihood
NLL <- numeric(nrow(par_sets))
for(i in 1:nrow(par_sets)){
  par_init$log_z <- log(par_sets[i,1])
  par_init$log_q2 <- log(par_sets[i,2])
  
  obj.sspm <- MakeADFun(data, par_init, random="log_Dep", map=map, DLL="sspm2")
  opt<-nlminb(obj.sspm$par, obj.sspm$fn, obj.sspm$gr)
  obj <- obj.sspm
  NLL[i] <- opt$objective
}
par_sets <-cbind(par_sets, NLL)
colnames(par_sets) <- c("z", "q2", "NLL")
par_sets



# Graphical output of Pop. dunamics and Kobe-plot
pdf("Res_Graph/Res_OMs_dynamics_kobe.pdf")
coltmp <- c("#cc00cc", "#9900cc", "#3300cc", "#0066cc",
            "#ff00ff", "#bf00ff", "#4000ff", "#0080ff",
            "#ff33ff", "#cc33ff", "#6633ff", "#3399ff")
pchtmp <- c(rep(9,4),rep(10,4),rep(12,4))
legend.name <- NULL
plot(c(0,X2), c(0,Y2), type = "n", xlab = "B / Bmsy", ylab = "F / Fmsy")
polygon(c(X1,1,1,X1), c(Y1,Y1,1,1), col = rgb(1,0.5,0,alpha=0.5), border = NA)
polygon(c(1,X2,X2,1), c(Y1,Y1,1,1), col = rgb(0,1,0,alpha=0.5), border = NA)
polygon(c(X1,1,1,X1), c(1,1,Y2,Y2), col = rgb(1,0,0,alpha=0.5), border = NA)
polygon(c(1,X2,X2,1), c(1,1,Y2,Y2), col = rgb(1,1,0,alpha=0.5), border = NA)
title(main = "Summary: kobe-plot of OMs")
for(i in 1:12){
  file <- paste0("Res_conditioning/z", par_sets[i,1]*100, "q2", par_sets[i,2]*100, ".rds")
  res <- readRDS(file)
  res1_est <- res$value
  MSY_est <- res1_est[3] %>% as.numeric()
  Bmsy_est <- res1_est[4] %>% as.numeric()
  B.est <- res$value[names(res$value)=="B"]
  points(last(B.est)/Bmsy_est,
         (last(dat$catch_B %>% na.omit())/last(B.est))/(MSY_est/Bmsy_est),
         col = coltmp[i], cex = 2, pch = pchtmp[i])
  legend.name <- c(legend.name, paste0("OM",i,": z=",par_sets[i,1]," q2=",par_sets[i,2]))
}
legend("topleft", legend.name, col = coltmp, pch = pchtmp, ncol=3)

X1 <- 0 ; X2 <- 3 ;  Y1 <- 0 ; Y2 <- 1.5
par(mfrow=c(1,2), oma=c(.5,.5,4,.5))
for(i in 1:12){
  file <- paste0("Res_conditioning/z", par_sets[i,1]*100, "q2", par_sets[i,2]*100, ".rds")
  res <- readRDS(file)
  res1_est <- res$value
  res1_sd <- res$sd
  r_est <- res1_est[1] %>% as.numeric()
  K_est <- res1_est[2] %>% as.numeric()
  MSY_est <- res1_est[3] %>% as.numeric()
  Bmsy_est <- res1_est[4] %>% as.numeric()
  q1_est <- res1_est[105] %>% as.numeric()
  q2_est <- res1_est[106] %>% as.numeric()
  sig1_est <- res1_est[107] %>% as.numeric()
  sig2_est <- res1_est[108] %>% as.numeric()
  tau_est <- res1_est[109] %>% as.numeric()
  B.est <- res$value[names(res$value)=="B"]
  B.se <- res$sd[names(res$value)=="B"]
  B.LB <- B.est*exp(-1.96*B.se/B.est)
  B.UB <- B.est*exp(+1.96*B.se/B.est)
  
  Year <- dat$year[21:70]
  Title_name <- paste0("OM",i,": z=",par_sets[i,1]," q2=",par_sets[i,2])
  plot(Year, B.est, type="l", col="red", ylab="Biomass (tons)", 
       ylim=c(0, max(B.UB)), main = "Pop. dynamics")
  points(Year, B.LB, type="l", lty=2, col="red")
  points(Year, B.UB, type="l", lty=2, col="red")
  polygon(c(Year,rev(Year)), c(B.LB,rev(B.UB)), 
          border = gray(.5,alpha=.5), col = gray(.5,alpha=.2))
  points(Year, data$Catch, type = "h", lwd=5, col="darkgray", ylab="Tons")
  abline(h=K_est, col="blue")
  
  plot(c(0,X2), c(0,Y2), type = "n", xlab = "B / Bmsy", ylab = "F / Fmsy")
  title(main = "Kobe-plot")
  polygon(c(X1,1,1,X1), c(Y1,Y1,1,1), col = rgb(1,0.5,0,alpha=0.5), border = NA)
  polygon(c(1,X2,X2,1), c(Y1,Y1,1,1), col = rgb(0,1,0,alpha=0.5), border = NA)
  polygon(c(X1,1,1,X1), c(1,1,Y2,Y2), col = rgb(1,0,0,alpha=0.5), border = NA)
  polygon(c(1,X2,X2,1), c(1,1,Y2,Y2), col = rgb(1,1,0,alpha=0.5), border = NA)
  points(last(B.est)/Bmsy_est,
         (last(dat$catch_B %>% na.omit())/last(B.est))/(MSY_est/Bmsy_est),
         col = "blue", cex = 2.5, pch = 10)
  mtext(Title_name, side = 3, line = 1, cex = 1.5, outer = TRUE)
}

dev.off()








#############################################################################
#調査漁具能率の違い
plot(0, 0, type = "n", xlim = range(Year), ylim = c(0,7500), xlab = "Year", ylab = "Abundance")
lty.tmp <- 1:3
for(i in 1:3){
  file <- paste0("Res_conditioning/z", 100, "q2", par_sets[i*4,2]*100, ".rds")
  res <- readRDS(file)
  res1_est <- res$value
  B.est <- sdr.sspm$value[names(sdr.sspm$value)=="B"]
  points(Year, B.est, type="l", lty=lty.tmp[i], col="red")
}


plot(z.tmp, LL.tmp, type = "l", cex=1.5, lwd = 3)



#############################################################################
#plot log-likelihood
z.tmp <- seq(0.001,1.5,length=20)
LL.tmp <- numeric(length(z.tmp))
for(i in 1:length(z.tmp)){
  par_init$log_z <- log(z.tmp[i])
  obj.sspm <- MakeADFun(data, par_init, random="log_Dep", map=map, DLL="sspm2")
  opt<-nlminb(obj.sspm$par, obj.sspm$fn, obj.sspm$gr)
  
}
plot(z.tmp, LL.tmp, type = "l", cex=1.5, lwd = 3)


tmp <- c(1:10,seq(10,5,length=10))
tmp2 <- numeric(length(tmp))
for(i in 1:20)tmp2[i] <- rlnorm(1,log(tmp[i]),1)
plot(tmp2, type = "l")

D <- pnorm(scale(tmp2), 0, 1)
last(D)
pnorm(mean(tail(D,3)),mean(D),sd(D))



#############################################################################
# surplus yield curve
surplus.func <- function(x, z, r=0.5, K=1) (r/z)*x*(1-(x/K)^z)

curve(surplus.func(x,z=0.01), xlim=c(0,1),lwd=3, col="black",yaxt="n",ylab="", xlab="", cex=3)
par(new=T)
curve(surplus.func(x,z=0.1), xlim=c(0,1),lwd=3, col="red",yaxt="n",ylab="", xlab="", cex=3)
par(new=T)
curve(surplus.func(x,z=0.5), xlim=c(0,1),lwd=3, col="orange",yaxt="n",ylab="", xlab="", cex=3)
par(new=T)
curve(surplus.func(x,z=1), xlim=c(0,1),lwd=3, col="pink",yaxt="n",ylab="", xlab="", cex=3)
par(new=T)
curve(surplus.func(x,z=1.5), xlim=c(0,1),lwd=3, col="purple",yaxt="n",ylab="", xlab="", cex=3)
abline(v=0.5, col="pink", lty=2)






#############################################################################
# population dynamics trajectories

z_sets <- c(0.1,0.5,1)
q2_sets <- c(0.35, 1, 1.5)
par_sets <- expand.grid(z = z_sets, q2 = q2_sets)
dep_mat <- biomass_mat <- matrix(NA, nrow = 51, ncol = nrow(par_sets))
for(i in 1:nrow(par_sets)){
  file <- paste0("MCMCsamples/20210813/PM_PTz-", par_sets[i,1], "_q2-", par_sets[i,2], "_mcmc.RDS")
  res <- readRDS(file)
  res <- res$samples[201:300,,]
  for(j in 1:51){
    biomass_mat[j,i] <- median(as.numeric(exp(res[,,3]))*as.numeric(exp(res[,,7+j])))
    dep_mat[j,i] <- median(as.numeric(exp(res[,,7+j])))
  }
}
biomass_mat <- cbind.data.frame(Year = dat$year[21:71], biomass_mat)
colnames(biomass_mat)[2:10] <- str_c("OM", 1:9)
dep_mat <- cbind.data.frame(Year = dat$year[21:71], dep_mat)
colnames(dep_mat)[2:10] <- str_c("OM", 1:9)

biomass_longer <- pivot_longer(biomass_mat,
                               cols = -"Year",
                               names_to = "OM",
                               values_to = "Biomass")
biomass_longer <- mutate(biomass_longer,
                         m = ifelse(OM=="OM1",1.1,
                                    ifelse(OM=="OM2",1.5,
                                           ifelse(OM=="OM3",2.0,
                                                  ifelse(OM=="OM4",1.1,
                                                         ifelse(OM=="OM5",1.5,
                                                                ifelse(OM=="OM6",2.0,
                                                                       ifelse(OM=="OM7",1.1,
                                                                              ifelse(OM=="OM8",1.5,2.0)
                                                                       ))))))),
                         q2 = ifelse(OM=="OM1",0.35,
                                    ifelse(OM=="OM2",0.35,
                                           ifelse(OM=="OM3",0.35,
                                                  ifelse(OM=="OM4",1.0,
                                                         ifelse(OM=="OM5",1.0,
                                                                ifelse(OM=="OM6",1.0,1.5)
                                                                       )))))
                         )
dep_longer <- pivot_longer(dep_mat,
                               cols = -"Year",
                               names_to = "OM",
                               values_to = "Depletion")
dep_longer <- mutate(dep_longer,
                     m = ifelse(OM=="OM1",1.1,
                                ifelse(OM=="OM2",1.5,
                                       ifelse(OM=="OM3",2.0,
                                              ifelse(OM=="OM4",1.1,
                                                     ifelse(OM=="OM5",1.5,
                                                            ifelse(OM=="OM6",2.0,
                                                                   ifelse(OM=="OM7",1.1,
                                                                          ifelse(OM=="OM8",1.5,2.0)
                                                                   ))))))),
                     q2 = ifelse(OM=="OM1",0.35,
                                 ifelse(OM=="OM2",0.35,
                                        ifelse(OM=="OM3",0.35,
                                               ifelse(OM=="OM4",1.0,
                                                      ifelse(OM=="OM5",1.0,
                                                             ifelse(OM=="OM6",1.0,1.5)
                                                      )))))
)

ggplot(biomass_longer, 
       aes(x=Year, y=Biomass, color=as.character(q2)))+
  geom_line(aes(group=OM), alpha=0.8, size=3)+
  #geom_point(aes(group=OM), alpha=0.8, size=5)+
  labs(color="q2")+
  ylim(0,11000)+
  facet_wrap(~m)+
  ylab("資源量(トン)")+xlab("時系列")+
  theme_bw()+
  theme(text = element_text(size=40))

ggplot(dep_longer, 
       aes(x=Year, y=Depletion, color=as.character(m), linetype=as.character(q2)))+
  geom_line(aes(group=OM), size=2)+
  scale_linetype_manual(values=c("0.35"=3,"1"=2,"1.5"=1))

