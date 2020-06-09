######################
### LOAD FUNCTIONS ###
######################
source('fun.R')

############################################
### GET PRIOR FROM CODIACS DATA ANALYSIS ###
############################################
#source('codiacs analysis semiparametric.R')
#pfoo0 = getRandProb(foo,base=2)
#print(pfoo0)
#pi10 = pfoo0$PI1[1]
#pi1 = pi11 = pfoo0$PI1[2]
#pi20 = pfoo0$PI2[,3]
#pi2 = pi21 = pfoo0$PI2[,4]
### above was commented out
pi1 = pi11 = 0.6749323
pi10 = 1 - pi11
pi2 = pi21 = c(0.7016099, 0.3799304, 0.4030912, 0.2576907)
pi20 = 1 - pi21

##############################################
#### TRUE REGRESSION AND MODEL PARAMETERS ####
##############################################
Theta1 = c(2.206, 5.594, 8.294, 7.745, -12.102, 6.455, -13.045, 6.59)     # Scene 1
Theta2 = c(2.206, 5.594, 8.294, 7.745, -6.051, 6.455, -6.5225, 0.0675)    # Scene 2
Theta3 = c(2.206, 5.594, 8.294, 7.745, -12.102, -6.455, -13.045, 19.5)    # Scene 3
Theta4 = c(2.206, 5.594, 8.294, 7.745, -6.051, -6.455, -6.5225, 12.9775)  # Scene 4
Theta5 = c(1.320, 6.480, 9.180, 9.555, -11.822, 4.645, -14.855, 6.382)    # Scene 5 (fitted saturated model)
Theta6 = c(2.206, 5.594, 8.294, 7.745, -12.102, -6.455, -13.045, -19.5)    # Scene 6 (modified scene 3)


theta = Theta3


Q2sat = function(a1,resp,a2, theta0 = theta) {
	dim(theta0) = c(8,1)
	x = cbind(1, a1, a2, resp, a1*a2, a1*resp, a2*resp, a1*a2*resp)
	return( x %*% theta0)
}
p0 = 29/56
p1 = 28/52
sigma = 6.7

#############################################
#### DERIVE VALUE OF ALL DTRs AND        ####
#### EVALUATE THE OPTIMAL DTR			 ####
#############################################
V0 = matrix(rep(NA,4*8),nrow=8)
V0[,1] = c(rep(0,4),rep(1,4))
V0[,2] = c(0,0,1,1,0,0,1,1)
V0[,3] = c(0,1,0,1,0,1,0,1)
for (i in 1:8) {
	d1 = V0[i,1]
	if (d1==1) w1 = p1
	if (d1==0) w1 = p0
	w0 = 1 - w1
	V0[i,4] = w0 * Q2sat(d1,0,V0[i,2]) + w1 * Q2sat(d1,1,V0[i,3])
}
Vest = V0[,4]
pos = which(Vest==max(Vest))
dtrBest = V0[pos,1:3]

cat("REGRESSION PARAMETERS OF Saturated model:\n")
thetaprn = matrix(theta,nrow=1)
colnames(thetaprn) = c("beta0","beta1","beta2","gamma1","beta3","gamma3","gamma2","gamma4")
print(thetaprn)
cat("\nAll possible DTRs and their values","\n")
colnames(V0) = c("d1","d2(0)","d2(1)","Value")
V0row = rep("",nrow(V0))
V0row[pos] = "OPTIMAL"
pos = which(Vest==min(Vest))
V0row[pos] = "WORST"
rownames(V0) = V0row
print(V0,digits=3)

###########################
#### DESIGN PARAMETERS ####
###########################
n <- 100
Nmin = n1 = 20
n2 = n - n1
base = 10
tau = 0.75
tau = n1 * tau^{1/(base-1)}  # slightly different from the manuscript notation (to speed up simulation a little bit)


############################################
#### PARAMETERS FOR MONITORING SCHEDULE ####
############################################
accrate = 4   # number of patients per month
obswin = 6  # observations will become available at calendar time (arrival + obswin) = Final outcome evaluation


set.seed(0202)
nsim = 5
DTRg = DTRhat = matrix(rep(NA,nsim*3),nrow=nsim)
valg = val = rep(NA,nsim)
Y = matrix(rep(NA,nsim*n),nrow=nsim)

for (r in 1:nsim) {
	if (round(r/100)==(r/100)) print(r)
	
	### Patient outcomes and arrival times
	eps = rnorm(n,0,sigma)
	arrival = c(0, cumsum(rexp(n-1,rate=accrate)))

	### Sequential experiment (a1, a2)
	a1 = a2 = R = y = rep(NA,n)
	for (i in 1:n) {
		indcomp = which(arrival < (arrival[i] - obswin))
		ncomp = length(indcomp)
		if (ncomp < n1)  {   # AR does not begin until there are n1 complete observations
			a1[i] = rbinom(1,1,pi1)  # pi1(initial random allocation prob) = 0.67 as determined from CODIACs
			if (a1[i]==0)  R[i] = rbinom(1,1,p0)
			else R[i] = rbinom(1,1,p1)
			if (a1[i]==0 & R[i]==0) a2[i] = rbinom(1,1,pi2[1])
			else if (a1[i]==0 & R[i]==1) a2[i] = rbinom(1,1,pi2[2])
			else if (a1[i]==1 & R[i]==0) a2[i] = rbinom(1,1,pi2[3])
			else a2[i] = rbinom(1,1,pi2[4])
		}
		else {
			y = Q2sat(a1,R,a2) + eps
			ycomp = y[indcomp]
			a1comp = a1[indcomp]
			a2comp = a2[indcomp]
			Rcomp = R[indcomp]
			foo = cbind(a1comp,Rcomp,a2comp,ycomp)
			pfoo = getRandProb(foo,base=base)
			w0 = (tau/ncomp)^{base-1}
			w1 = 1-w0	
			pi1hat = pi11hat = pfoo$PI1[2]
			pi10hat = pfoo$PI1[1]
			rho10 = exp(w0*log(pi10) + w1*log(pi10hat))
			rho11 = exp(w0*log(pi11) + w1*log(pi11hat))
			pi1til = rho11  / (rho10 + rho11)		

			pi20hat = pfoo$PI2[,3]
			pi2hat = pi21hat = pfoo$PI2[,4]
			rho20 = exp(w0*log(pi20) + w1*log(pi20hat))
			rho21 = exp(w0*log(pi21) + w1*log(pi21hat))
			pi2til = rho21 / (rho20+rho21)

			a1[i] = rbinom(1,1,pi1til)
			if (a1[i]==1)  R[i] = rbinom(1,1,p1)
			else R[i] = rbinom(1,1,p0)
			if (a1[i]==0 & R[i]==0) a2[i] = rbinom(1,1,pi2til[1])
			else if (a1[i]==0 & R[i]==1) a2[i] = rbinom(1,1,pi2til[2])
			else if (a1[i]==1 & R[i]==0) a2[i] = rbinom(1,1,pi2til[3])
			else a2[i] = rbinom(1,1,pi2til[4])
		}
	}
	
	### Final fit after complete follow-up of all subjects
	y = Q2sat(a1,R,a2) + eps
	foo = cbind(a1,R,a2,y)
	fit1 = learning(foo)

	Y[r,] = y
	DTRhat[r,] = dtrhat = fit1$odtr
	pos = which(V0[,1]==dtrhat[1] & V0[,2]==dtrhat[2] & V0[,3]==dtrhat[3])
	val[r] = V0[pos,4]
	DTRg[r,] = getGestimate(foo)$odtr	
	pos = which(V0[,1]==DTRg[r,1] & V0[,2]==DTRg[r,2] & V0[,3]==DTRg[r,3])
	valg[r] = V0[pos,4]
}


cat("\nSMART-AR continuous monitoring\n")
cat("N:\t",n,"\n")
cat("Nmin:",n1,"\n")
cat("base:",base,"\n")
cat("tau:",tau,"\n")
cat("pi1:",pi1,"\n")
cat("pi2:",pi2,"\n")
cat("arate:", accrate,"\n")
cat("nsim:",nsim,"\n\n")

v4 = V0[,4]
vmax = which(v4==max(v4))
vmax = V0[vmax,4]
vmin = which(v4==min(v4))
vmin = V0[vmin,4]
pcs = length(which(abs(val-vmax)<=0.0001))/nsim
pcsg = length(which(abs(valg-vmax)<=0.0001))/nsim

cat("Distribution of the values of the selected DTR (dhat):\n")
Value.dhat = val
print(summary(Value.dhat))
print(table(Value.dhat),digits=3)
cat("\n")

cat("Distribution of the average patient outcome per trial\n")
patientval = apply(Y,1,mean)
print(summary(patientval))

#cat("Some key summary statistics\n")
cat("\nProbability of selecting the true optimal DTR:", pcs, "\n\n")
cat("Value of selected dtr (Mean):", mean(val),"\n")
cat("Value of selected dtr (SD):", sd(val),"\n")
cat("Value of selected dtr (Q1,Q2,Q3)", quantile(val,c(0.25,0.5,0.75)),"\n")
cat("Adjusted value of selected dtr (Mean)", (mean(val)-vmin)/(vmax-vmin),"\n")
cat("Adjusted value of selected dtr (SD)", sd(val)/(vmax-vmin),"\n")
cat("Adjusted value of selected dtr (Q1,Q2,Q3)", quantile(val-vmin,c(0.25,0.5,0.75))/(vmax-vmin),"\n\n")


cat("Average patient outcome (Mean):", mean(patientval),"\n")
cat("Average patient outcome (SD):", sd(patientval),"\n")
cat("Average patient outcome (Q1,Q2,Q3):", quantile(patientval,c(0.25,0.5,0.75)),"\n")
cat("Adjusted average patient outcome (Mean):", (mean(patientval)-vmin)/(vmax-vmin),"\n")
cat("Adjusted average patient outcome (SD):", sd(patientval)/(vmax-vmin),"\n")
cat("Adjusted average patient outcome (Q1,Q2,Q3):", quantile(patientval-vmin,c(0.25,0.5,0.75))/(vmax-vmin),"\n")


cat("\nTHE RESULTS BELOW ARE BASED ON G-estimation\n") 


cat("Distribution of the values of the selected DTR based on G-estimation (dhat):\n")
Value.dhatg = valg
print(summary(Value.dhatg))
print(table(Value.dhatg),digits=3)
cat("\n")


cat("Probability of selecting the true optimal DTR:", pcsg, "\n\n")
cat("Value of selected dtr (Mean):", mean(valg),"\n")
cat("Value of selected dtr (SD):", sd(valg),"\n")
cat("Value of selected dtr (Q1,Q2,Q3)", quantile(valg,c(0.25,0.5,0.75)),"\n")
cat("Adjusted value of selected dtr (Mean)", (mean(valg)-vmin)/(vmax-vmin),"\n")
cat("Adjusted value of selected dtr (SD)", sd(valg)/(vmax-vmin),"\n")
cat("Adjusted value of selected dtr (Q1,Q2,Q3)", quantile(valg-vmin,c(0.25,0.5,0.75))/(vmax-vmin),"\n\n")
