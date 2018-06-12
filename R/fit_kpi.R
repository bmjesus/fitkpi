#' @title Function for determining Kpi from FRRF data
#' @description  Function for determining Kpi from FRRF data
#' @param file_linc path to the file with the lincomycin data
#' @param file_control path to the file with the control data (without inhibitor)
#' @param file_calib path to the calibration file (specific for each machine)
#' @param fit_start light step where to start the fit
#' @param fit_finish light step where to finish the fit
#' @param recovery_start light step where the recovery period starts
#' @param recovery_finish light step where the recovery period finishes
#' @param light_steps number of light steps inside the file, default=10
#' @return prints the results in the console and produces a text file with them also
#' @export
fit_kpi<-function(file_linc,file_control,file_calib,
                  fit_start = 2,
                  fit_finish = 8,
                  recovery_start = 8,
                  recovery_finish = 10,
                  light_steps= 10
)
{

#this imports the data from a lincomycin treated sample
#it uses the fit_model_fofm from the psifluo package that processes
#data from the FL3500 PSI fluorometers

psifluo::fit_model_fofm(file_name = file_linc,
               num_steps = light_steps,
               calib_file = file_calib,
               out_name = "marlene_linco"
  )

################################################################################
#GOAL 1 - plot Fv/Fm 2s dark + lincomycin vs. time
#single phase exponential decay = kpi
#(first order rate constant for photoinactivation at the given irradiance level)

#using the number of light steps to extract data corresponding to the 2s dark measurements
select_i <- light_steps + 1
select_f <- length(etr_fit$eff)

#extracting the Fv/Fm data of the 2s dark measurements
fv_fm_linco<-etr_fit$eff[select_i:select_f]


#extract the time
time_simple<-strptime(as.character(time.data$V3), format="%m/%d/%Y %I:%M:%S %p")

#convert it to minutes from starting point

t_init<-time_simple[1]
t_abs<-numeric()
a<-numeric()
for (i in 2:length(time_simple))
{
  a<-time_simple[i]-t_init
  t_abs<-c(t_abs,a)
}
t_abs<-c(0,t_abs)

#convert time to seconds
t_abs<-t_abs*60


#data subselection between 2 and 8 light steps (456 ppfd)
#T0 for high light treatment = step 2, then 6 steps at 456 takes to step 8
x<-t_abs[fit_start:fit_finish] - t_abs[fit_start]
y<-fv_fm_linco[fit_start:fit_finish]

#correct y (fv_fm) for influence of low light recovery periods 9 & 10 after high light periods 3-8.
y_corr <- c(y[1],(y[fit_start:(fit_finish-1)] + (fv_fm_linco[recovery_finish]-fv_fm_linco[recovery_start ])))

#generate a shorter set of y_corr points to cover the initial decline over periods 2,3,4,5
fit_final<-(fit_finish-fit_start)-1
y_trim <-y_corr[1:fit_final]
x_trim <-x[1:fit_final]

#fitting the single phase exponential decay and estimating Kpi, over both y_corr & y_trim
#fit_kpi<-nlsLM(y_corr ~ y_corr[1]*exp(-kpi*x),  start=list(kpi=0.0001))

fit_kpi_trim<-minpack.lm::nlsLM(y_trim ~ y_trim[1]*exp(-kpi*x_trim),  start=list(kpi=0.0001))

#extracting the value from the fitted model
# p_kpi<-summary(fit_kpi)$parameters[1]
p_kpi_trim<-summary(fit_kpi_trim)$parameters[1]

#predicting values using the kpi from the model,
#this is useful info to overlay with original data and check goodness of fit
# pred_y<-y_corr[1]*exp(-p_kpi*x)
pred_y_trim<-y_trim[1]*exp(-p_kpi_trim*x_trim)

#it's a big number, better to reformat
kpi_legend<-formatC(p_kpi_trim, format = "e", digits = 2)

#storing the various axis so that I can use them later

x_trim1<-x_trim
y_trim1<-y_trim

#correcting the starting point for the trimmed data by re-adding the time
x_trim1<-x_trim1+t_abs[fit_start]


#storing the plot inside a function so that I can build a multiplot figure at the end
plot1<-function(){
  plot(t_abs,fv_fm_linco,xlab="",ylab = "",las=1,pch=21,bg=1)
  #overlay with predicted values and y_trim values
  points(x_trim1,pred_y_trim,col=2,lty=2,type='l')
  points(x_trim1, y_trim1,pch=21,bg=2)

#legend('topright',legend = c('original','trim-corrected'),pch=21,col=c(1,2),pt.bg = c(1,2),bty='n',cex=0.5)

  mtext(side = 2, "Fv/Fm",las=3,cex=0.8,line=2.5)
  mtext(side = 1, "time (seconds)",las=1,cex=0.8,line=2)
  mtext(side =3, paste('Kpi = ',kpi_legend),line=-1.5)
}



################################################################################
#GOAL 2
#plot Fv/Fm 2s dark + lincomycin vs. cumulative photons m-2  (= photons m-2 s-1 x elapsed time)
#single phase exponential decay = sigma i
#(the probability of photoinactivation per incident photon)


#extract light intensities to a separate object
light_levels<-volt$par

#calculating light accumulation at each step, assuming time step does not change
#between steps which I think cannot be done with the PSI script anyway,
#so we are ok with this assumption
t_step<-t_abs[fit_start]
light_steps<-light_levels*t_step*6.022*10^17

#calculating accumulated light
l_cum<-cumsum(light_steps)

####trimming and correct the xaxis as above
#data subselection between 2 and 8 light steps (456 ppfd)
#T0 for high light treatment = step 2, then 6 steps at 456 takes to step 8
x<-l_cum[fit_start:fit_finish] - l_cum[fit_start]
x_trim <-x[1:fit_final]

#single phase exponential decay = sigma i
#(the probability of photoinactivation per incident photon)
fit_sigma<-minpack.lm::nlsLM(y_trim ~ y_trim[1]*exp(-k_sigma*x_trim),  start=list(k_sigma=1*10^-25))

p_sigma<-summary(fit_sigma)$parameters[1]

pred_y_sigma<-y_trim[1]*exp(-p_sigma*x_trim)

#it is a huge number so best to express it scientifically
sigma_legend<-formatC(p_sigma, format = "e", digits = 2)

#storing the various axis so that I can use them later
x2<-x
y2<-y
x_trim2<-x_trim
y_trim2<-y_trim

#correcting the starting point for the trimmed data by re-adding light
x_trim2<-x_trim2+l_cum[fit_start]


#storing the plot inside a function so that I  can build a multiplot figure at the end
plot2<-function(){
plot(l_cum,fv_fm_linco,xlab="",ylab = "Fv/Fm",las=1,pch=21,bg=1,yaxt="n")
points(x_trim2,pred_y_sigma,type='l',lty=2,col=2)
points(x_trim2,pred_y_sigma,pch=21,bg=2)

#legend('topright',legend = c('original','trim-corrected'),pch=21,col=c(1,2),pt.bg = c(1,2),bty='n',cex=0.5)

axis(4)
mtext(side = 4, "Fv/Fm",las=3,cex=0.8,line=2.5)
mtext(side = 1, "cumulative photons m-2",las=1,cex=0.8,line=2)
mtext(side =3, paste('sigma = ',sigma_legend),line=-1.5)

}


################################################################################
#GOAL 3
#plot Fv/Fm 2s dark + lincomycin vs. cumulative photons PSII-1  (= photons m-2 s-1 x elapsed time x sigmaPSII)
#single phase exponential decay = sigma psii i
#(the probability of photoinactivation per photon delivered to PSII)

#storing the sigma PSII in a separate object
sigmaPSII<-etr_fit$sigma[select_i:select_f]

#calculating photons PSII-1 at each light level
light_steps_sigmaPSII<-light_levels*t_step*sigmaPSII*6.022*10^17*1*10^(-20)

#calculating cumulative photons PSII-1
l_cum_sigmaPSII<-cumsum(light_steps_sigmaPSII)

####trimming and correct the xaxis as above
#data subselection between 2 and 8 light steps (456 ppfd)
#T0 for high light treatment = step 2, then 6 steps at 456 takes to step 8
x<-l_cum_sigmaPSII[fit_start:fit_finish] - l_cum_sigmaPSII[fit_start]
x_trim <-x[1:fit_final]

#single phase exponential decay = sigma i
#(the probability of photoinactivation per incident photon)
fit_sigmaPSII<-minpack.lm::nlsLM(y_trim ~ y_trim[1]*exp(-k_sigmaPSII*x_trim),  start=list(k_sigmaPSII=1*10^-08))

p_sigmaPSII<-summary(fit_sigmaPSII)$parameters[1]

pred_y_sigmaPSII<-y_trim[1]*exp(-p_sigmaPSII*x_trim)

sigmaPSII_legend<-formatC(p_sigmaPSII, format = "e", digits = 2)


#storing the various axis so that I can use them later
x3<-x
y3<-y
x_trim3<-x_trim
y_trim3<-y_trim

#correcting the starting point for the trimmed data by re-adding light
x_trim3<-x_trim3+l_cum_sigmaPSII[fit_start]

#storing the plot inside a function so that I build a multiplot figure at the end
plot3<-function(){
plot(l_cum_sigmaPSII,fv_fm_linco,xlab="",ylab = "",las=1,pch=21,bg=1)

points(x_trim3,pred_y_sigmaPSII,type='l',lty=2,col=2)
points(x_trim3,pred_y_sigmaPSII,pch=21,bg=2)

#legend('topright',legend = c('original','trim-corrected'),pch=21,col=c(1,2),pt.bg = c(1,2),bty='n',cex=0.5)


mtext(side = 2, "Fv/Fm",las=3,cex=0.8,line=2.5)
mtext(side = 1, "cumulative photons m-2.photons PSII-1",las=1,cex=0.8,line=2)
mtext(side =3, paste('sigmaPSII i = ',sigmaPSII_legend),line=-1.5)


}

################################################################################
#GOAL 4
#Estimate krec
#plot Fv/Fm2s dark no lincomycin vs. time
#Use Kok equation to fit for krec, using  kpi as an input.
################################################################################
#change file to the control file
psifluo::fit_model_fofm(file_name = file_control,
               num_steps = light_steps,
               calib_file = file_calib,
               out_name = "marlene_control"
)


#extracting the Fv/Fm data of the 2s dark measurements
fv_fm_control<-etr_fit$eff[select_i:select_f]



#extract the time
time_simple<-strptime(as.character(time.data$V3), format="%m/%d/%Y %I:%M:%S %p")

#convert it to minutes from starting point
t_init<-time_simple[1]
t_abs<-numeric()
a<-numeric()
for (i in 2:length(time_simple))
{
  a<-time_simple[i]-t_init
  t_abs<-c(t_abs,a)
}
t_abs<-c(0,t_abs)

#convert time to seconds
t_abs<-t_abs*60

#data subselection between 2 and 8 light steps (456 ppfd)
#T0 for high light treatment = step 2, then 6 steps at 456 takes to step 8
x<-t_abs[fit_start:fit_finish] - t_abs[fit_start]
y<-fv_fm_control[fit_start:fit_finish]

#correct y (fv_fm) for influence of low light recovery periods 9 & 10 after high light periods 3-8..
y_corr <- c(y[1],(y[fit_start:(fit_finish-1)] + (fv_fm_linco[recovery_finish]-fv_fm_linco[recovery_start])))

#generate a shorter set of y_corr points cover the initial decline over periods 2,3,4,5
y_trim <-y_corr[1:fit_final]
x_trim <-x[1:fit_final]


fit_krec<-minpack.lm::nlsLM(y_trim~y_trim[1]*((krec+(p_kpi_trim*((exp(-(p_kpi_trim+krec)*x_trim)))))/(p_kpi_trim+krec)),start=list(krec=0.0001))

p_krec<-summary(fit_krec)$parameters[1]

pred_krec<-y_trim[1]*((p_krec+(p_kpi_trim*((exp(-(p_kpi_trim+p_krec)*x_trim)))))/(p_kpi_trim+p_krec))

krec_legend<-formatC(p_krec, format = "e", digits = 2)


#storing the various axis so that I can use them later
x4<-x
y4<-y
x_trim4<-x_trim
y_trim4<-y_trim

#correcting the starting point for the trimmed data by re-adding time
x_trim4 <- x_trim4 + t_abs[fit_start]


#storing the plot inside a function so that I can build a multiplot figure at the end
plot4<-function(){
plot(t_abs,fv_fm_control,xlab="",ylab='',las=1,
     main = ,pch=21,bg=1,yaxt='n')

points(x_trim4,pred_krec,type = 'l',lty=2,col=2)
points(x_trim4,pred_krec,pch=21,bg=2)

legend(2000,0.54,legend = c('original','trim-corrected'),pch=21
       ,col=c(1,2),pt.bg = c(1,2),bty='n',cex=0.8,adj=0,
       y.intersp=0.4)


axis(4)
mtext(side = 4, "Fv/Fm",las=3,cex=0.8,line=2.5)
mtext(side = 1, "time (seconds)",las=1,cex=0.8,line=2)
mtext(side =3, paste('Krec = ',krec_legend),line=-1.5)


}

################################################################################
#####SECTION FOR PLOTTING
################################################################################

nf<-layout(matrix(c(1,2,3,4),ncol=2,byrow=TRUE))
#par(oma=c(4,4,3,4),mar=c(2,2,2,2), xpd=NA,tcl=-0.3,bg="white",cex=0.8,
#    cex.axis=0.9,cex.lab=0.9,bty="o",las=1,mgp=c(3,0.5,0),adj=0.5)

par(oma=c(2,4,2,4),mar=c(2,0,2,0),cex=0.8,cex.axis=0.9,cex.lab=0.9,bty="o"
    ,las=1)

plot1()
plot2()
plot3()
plot4()


################################################################################
#####SECTION FOR OUTPUT RESULTS
################################################################################

#printing result on the console
print(paste('Kpi = ',kpi_legend))
print(paste('sigma = ',sigma_legend))
print(paste('sigmaPSII_i = ',sigmaPSII_legend))
print(paste('Krec = ',krec_legend))

#write a csv with the results
output<-as.data.frame(cbind(c('Kpi','sigma','sigmaPSII_i','Krec'),
              c(kpi_legend,sigma_legend,sigmaPSII_legend,krec_legend)))

names(output)<-c('parameter','value')

write.table(output,file = 'Kpi_output.txt',sep=",",col.names = TRUE,row.names = FALSE)



}
