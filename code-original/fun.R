### Contains functions to generate, evaluate value, and perform Q-learning.
### Specifics: Each time 2 possible treatments; two stage; a binary intermediate response
### Stage-2 regression model
getRandProb = function(foo, base=2) {
	a1 = z1 = foo[,1]   # a1
	a2 = z2 = foo[,3]   # a2
	z3 = a1*a2     # a1*a2
	resp = z4 = foo[,2]   # resp
	z5 = resp*(1-a1)*a2
	z6 = resp*a1*(1-a2)

	y = foo[,4]
	n = length(y)
	fit = lm(y ~ z1+z2+z3+z4+z5+z6)
	coef = fit$coef
	coef[which(is.na(coef))] = 0
	yhat = rep(NA,n)  # pseudo outcome
	for (i in 1:n) {
		yhat[i] = Q2max(a1[i],resp[i],coef[1],coef[2],coef[3],coef[4],coef[5],coef[6],coef[7])[2]
	}
	fit1 = lm(yhat ~ a1)
	q1hat = c(mean(fit1$fit[a1==0]), mean(fit1$fit[a1==1]))  # values correspond to a1=0, a1=1
	d1 = q1hat - min(q1hat)
	d1 = d1 / summary(fit1)$sigma
	d1 = pmin(d1, log(100000)/log(base))
	PI1 = base^d1 / sum(base^d1)
	PI2 = matrix(rep(NA, 16), nrow=4)
	colnames(PI2) = c("a1","resp","prob(a2=0)","prob(a2=1)")
	iter = 1
	for (acur in 0:1) { 
		for (rcur in 0:1) {
			q2hat = c(Q2(acur,rcur,0,coef[1],coef[2],coef[3],coef[4],coef[5],coef[6],coef[7]),Q2(acur,rcur,1,coef[1],coef[2],coef[3],coef[4],coef[5],coef[6],coef[7]))
			d2 = q2hat - min(q2hat)
			d2 = d2 / summary(fit)$sigma
			d2 = pmin(d2, log(100000)/log(base))
			prand2 = base^d2 / sum(base^d2)
			PI2[iter,] = c(acur, rcur, prand2)
			iter = iter + 1
		}
	}
	return(list(PI1=PI1, PI2=PI2))
}

learning = function(foo) {
	a1 = z1 = foo[,1]   # a1
	a2 = z2 = foo[,3]   # a2
	z3 = a1*a2     # a1*a2
	resp = z4 = foo[,2]   # resp
	z5 = resp*(1-a1)*a2
	z6 = resp*a1*(1-a2)

	y = foo[,4]
	n = length(y)
	fit = lm(y ~ z1+z2+z3+z4+z5+z6)
	coef = fit$coef
	coef[which(is.na(coef))] = 0 
	yhat = rep(NA,n)  # pseudo outcome
	for (i in 1:n) {
		yhat[i] = Q2max(a1[i],resp[i],coef[1],coef[2],coef[3],coef[4],coef[5],coef[6],coef[7])[2]
	}
	fit1 = lm(yhat ~ a1)
	q1hat = c(mean(fit1$fit[a1==0]), mean(fit1$fit[a1==1]))  # values correspond to a1=0, a1=1
	q20hat = c(Q2max(0,0, coef[1],coef[2],coef[3], coef[4],coef[5],coef[6],coef[7])[2], Q2max(0,1,coef[1],coef[2],coef[3],coef[4],coef[5],coef[6],coef[7])[2])  ## values correspond to a1=0; r1=0,1
	
	
	q21hat = c(Q2max(1,0, coef[1],coef[2],coef[3], coef[4],coef[5],coef[6],coef[7])[2], Q2max(1,1,coef[1],coef[2],coef[3],coef[4],coef[5],coef[6],coef[7])[2])  ## values correspond to a1=0; r1=0,1
	
	odtr = rep(NA,3)
	if (q1hat[1] >= q1hat[2])  {
		odtr[1] = 0
		odtr[2:3] = c(Q2max(0,0, coef[1],coef[2],coef[3], coef[4],coef[5],coef[6],coef[7])[1], Q2max(0,1,coef[1],coef[2],coef[3],coef[4],coef[5],coef[6],coef[7])[1]) 
	}
	else { 
		odtr[1] = 1
		odtr[2:3] = c(Q2max(1,0, coef[1],coef[2],coef[3], coef[4],coef[5],coef[6],coef[7])[1], Q2max(1,1,coef[1],coef[2],coef[3],coef[4],coef[5],coef[6],coef[7])[1]) 
	}
	return(list(odtr=odtr, fit=summary(fit), fit1=summary(fit1),q1hat=q1hat))
}

Q2 = function(a1,resp,a2,b0,b1,b2,b3,g1,g2,g3) {
	b0 + b1*a1 + b2*a2 + b3*a1*a2 + g1*resp + g2*resp*(1-a1)*a2 + g3*resp*a1*(1-a2)
}


Q2max = function(a1,resp,b0,b1,b2,b3,g1,g2,g3) {
	Q20 = Q2(a1,resp,0,b0,b1,b2,b3,g1,g2,g3)
	Q21 = Q2(a1,resp,1,b0,b1,b2,b3,g1,g2,g3)
	if (Q20 >= Q21)  return( c(0,Q20) )
	else { return( c(1,Q21) ) }
}


Q1 = function(a1,b0,b1,b2,b3,g1,g2,g3,p) {  #not used in Qlearn because p = response rate under a1 is not available
	val = (1-p)*Q2max(a1,0, b0,b1,b2,b3,g1,g2,g3)[2] + p*Q2max(a1,1,b0,b1,b2,b3,g1,g2,g3)[2]
	return(val)
}

getTrueOptimalDTR = function(b0,b1,b2,b3,g1,g2,g3,p0, p1) {  # (a1,a20,a21) fir given true parameters
	q10 = Q1(0, b0,b1,b2,b3,g1,g2,g3,p0)	
	q11 = Q1(1, b0,b1,b2,b3,g1,g2,g3,p1)
	
	if (q10 >= q11) { 
		q20 = Q2max(0,0,b0,b1,b2,b3,g1,g2,g3)
		q21 = Q2max(0,1,b0,b1,b2,b3,g1,g2,g3)	
		
		return(list(DTR = c(0, q20[1], q21[1]), Q1 = q10, Q20 = q20[2], Q21 = q21[2]))
	}
	else {
		q20 = Q2max(1,0,b0,b1,b2,b3,g1,g2,g3)
		q21 = Q2max(1,1,b0,b1,b2,b3,g1,g2,g3)	
		return(list(DTR = c(1, q20[1], q21[1]), Q1 = q11, Q20 = q20[2], Q21 = q21[2]))
	}
}

getTrueValue = function(dtr,b0,b1,b2,b3,g1,g2,g3,p0,p1) {  # evaluate the value of a DTR for given true parameters
	a1 = dtr[1]
	a20 = dtr[2]  # a2 when r=0
	a21 = dtr[3]  # a2 when r=1
	if (a1==0)  p = p0
	else p  = p1
	
	y0 = Q2(a1,0,a20,b0,b1,b2,b3,g1,g2,g3)
	y1 = Q2(a1,1,a21,b0,b1,b2,b3,g1,g2,g3)
	val = y0*(1-p) + y1*p
	return(val)
}


getGestimate = function(foo) {
	a1 = foo[,1]
	a2 = foo[,3]
	resp = foo[,2]
	y = foo[,4]
	
	Vhat = matrix(rep(NA,5*8),nrow=8)  # Store DTR and the estimated values
	Vhat[,1] = c(rep(0,4),rep(1,4))
	Vhat[,2] = c(0,0,1,1,0,0,1,1)
	Vhat[,3] = c(0,1,0,1,0,1,0,1)
	for (r in 1:8) {
		d1 = Vhat[r,1]
		d20 = Vhat[r,2]
		d21 = Vhat[r,3]
		
		# G-estimation
		p1hat = length(which(a1==d1 & resp==1)) / length(which(a1==d1))
		p0hat = 1- p1hat
		m0 = mean(y[a1==d1 & resp==0 & a2==d20])
		m1 = mean(y[a1==d1 & resp==1 & a2==d21])
		gmu = Vhat[r,4] = p0hat * m0 + p1hat * m1
		# Simple mean
		ind = which(a1==d1 & ((resp==0&a2==d20) | (resp==1&a2==d21)))
		mu = Vhat[r,5] = mean(y[ind])
	}
	colnames(Vhat) = c("d1","d20","d21","Gest","Ave")
	Vest = Vhat[,4]
	pos = which(Vest==max(Vest,na.rm=T))
	return(list(odtr=Vhat[pos,1:3], Values = Vhat))
}

