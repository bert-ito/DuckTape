llsense_full<-function(N1,N2,p,sizefactors,dispfunc)
{
	lb1 <- sum(dbinom(N2,N1+N2,p,log=TRUE))
	SN <- N1 
	AN <- N2 
	TN <- SN + AN
	baseLevelTN <- TN/sizefactors
	TNq <- sum(baseLevelTN)/length(baseLevelTN)
	rawVarTN <- dispfunc(TNq) 
	TNmu <- TNq*sizefactors
	fullVarTN <- TNmu + rawVarTN * TNmu^2
	probTN <- TNmu/fullVarTN
	sizeTN <- TNmu^2/(fullVarTN-TNmu)
	lnb1 <- sum(dnbinom(TN,size=sizeTN,mu=TNmu, log=TRUE))
	totalln<-sum(lb1)+sum(lnb1)
	result<-list(SN=N1,
			AN=N2,
			TN=TN,
			TNmu=TNmu,
			FullVar=fullVarTN,
			probTN = probTN,
			sizeTN = sizeTN,
			lnb = lnb1,
			lb = lb1)
	return(totalln)
}

llanti_full<-function(x,N,p,sizefactors,dispfunc)
{
	# assuming balanced design
	nconds<-length(N)/2
	a<-x[1:nconds]; b<-x[(nconds+1):(nconds*2)]
	N1 <- N[1:nconds]; N2 <- N[(nconds+1):(nconds*2)]
	lb1 <- sum(dbinom(a,N1-b+a,p,log=TRUE));lb2 <- sum(dbinom(b,N2-a+b,p,log=TRUE))
	SN <- N1 - b + a;AN <- N2 - a + b
	baseLevelSN <- SN/sizefactors; baseLevelAN <- AN/sizefactors
	SNq <- sum(baseLevelSN)/length(baseLevelSN);ANq <- sum(baseLevelAN)/length(baseLevelAN) 
	rawVarSN <- max(dispfunc(max(SNq,1e-20)),0);rawVarAN <- max(dispfunc(max(ANq,1e-20)),0)
	SNmu <- SNq*sizefactors;ANmu <- ANq*sizefactors
	
	fullVarSN <- SNmu + rawVarSN * SNmu^2; fullVarAN <- ANmu + rawVarAN * ANmu^2
	probSN <- SNmu/fullVarSN; probAN <- ANmu/fullVarAN
	sizeSN <- SNmu^2/(fullVarSN-SNmu);sizeAN <- ANmu^2/(fullVarAN-ANmu)
	lnb1 <- pmax(dnbinom(SN,size=sizeSN,mu=SNmu, log=TRUE),-1e4);lnb2 <- pmax(dnbinom(AN,size=sizeAN,mu=ANmu, log=TRUE),-1e4)
	totalln<-sum(lb1)+sum(lb2)+sum(lnb1)+sum(lnb2)
	if (is.na(totalln))
	{return(-1e20)}
	return(totalln)
}

llanti_red<-function(x,N,p,sizefactors,dispfunc)
{
	#number of conditions is half of the N
	nconds<-length(N)/2
	a<-x[1:nconds]; b<-x[(nconds+1):(nconds*2)]
	N1 <- N[1:nconds]; N2 <- N[(nconds+1):(nconds*2)]
	lb1 <- sum(dbinom(a,N1-b+a,p,log=TRUE));lb2 <- sum(dbinom(b,N2-a+b,p,log=TRUE))
	SN <- N1 - b + a;AN <- N2 - a + b
	baseLevelSN <- SN/sizefactors; baseLevelAN <- AN/sizefactors
	SNq <- sum(baseLevelSN)/length(baseLevelSN);ANq <- sum(baseLevelAN)/length(baseLevelAN) 
	rawVarSN <- max(dispfunc(SNq),0);rawVarAN <- max(dispfunc(ANq),0)
	SNmu <- SNq*sizefactors;ANmu <- ANq*sizefactors
	
	fullVarSN <- SNmu + rawVarSN * SNmu^2; fullVarAN <- ANmu + rawVarAN * ANmu^2
	probSN <- SNmu/fullVarSN; probAN <- ANmu/fullVarAN
	sizeSN <- SNmu^2/(fullVarSN-SNmu);sizeAN <- ANmu^2/(fullVarAN-ANmu)
	lnb1 <- dnbinom(SN,size=sizeSN,mu=SNmu, log=TRUE);lnb2 <- dnbinom(AN,size=sizeAN,mu=ANmu, log=TRUE)
	lnb11 <- dnbinom(SN,size=sizeSN,prob=probSN, log=TRUE);lnb21 <- dnbinom(AN,size=sizeAN,prob=probAN, log=TRUE)
	
	totalln<-sum(lb1)+sum(lb2)+sum(lnb1)
	if (is.na(totalln))
	{return(-1e20)}
	return(totalln)
}

hillclimb<-function(xstart,lowb,highb,llf,...)
{
	#start are the starting values
	# lowb: lower bound of the values, for binomial and negative binomial, lowb should be 0
	# highb: high bound of the values, for binomial and negative binomial, highb should be N
	#source('~/Scripts/WorkSpace/AntisenseR/AntisenseR.R')
	args <- list(...) 

	N1<-args[[1]];N2<-args[[2]]
	p<-args[[3]]
	SF<-args[[4]]
	rawvarfunc<-args[[5]]

	
	nvar <- length(xstart)
	nlbound <- length(lowb)
	nhbound <- length(highb)
	if ((nvar!=nlbound)||(nvar!=nhbound))
	{
		warning( "number of bounds are different from number of variables." )
		break
	}
	step <-1
	count<-0
	delta<-c(1e-20) # minimum improvement.
	while(1){
		startll<-llf(xstart,c(N1,N2),p,SF,rawvarfunc)
		lastll<-startll
		xlast <- xstart
		for (i in 1:nvar)
		{
			# scan each variable
			higher<-xlast[i]+step
			lower<-xlast[i]-step
			if (higher>highb[i]){higher<-highb[i]}
			if (lower<lowb[i]){lower<-lowb[i]}
			xlasth<-xlast
			xlasth[i]<-higher
			higherll<-llf(xlasth,c(N1,N2),p,SF,rawvarfunc)
			
			xlastl<-xlast; xlastl[i]<-lower
			lowerll<-llf(xlastl,c(N1,N2),p,SF,rawvarfunc)
			if(higherll>=lowerll)
			{   # in case of tie, go with higher
				diff <- higherll - lastll
				if (diff<0)
				{   #no improvement
					next
				}
				else
				{
					xlast<-xlasth
					lastll<-higherll
					next
				}
			}
			else
			{
				diff <- lowerll - lastll
				if (diff<0)
				{   #no improvement
					next
				}
				else
				{
					xlast<-xlastl
					lastll<-lowerll
					next
				}
				
			}
		}
		if (lastll>startll)
		{
			diff <- lastll-startll
			if (diff<delta){return(xstart)}
			else{
				xstart <- xlast
				startll <- lastll
				next
			}
		}
		else
		{
			return(xstart)
		}
	}
	
}

sumll<-function(lowb,highb,llf=llantimodel_new,...)
{
	# this is a function to calculate the sum of likelihood.
	args <- list(...) 
	if (max(highb)>10000)
	{
		print('too slow to compute')
		return(0)
	}
	N<-args[[1]]
	p<-args[[2]]
	SF<-args[[3]]
	dispfunc<-args[[4]]
	out<-0
	out<-matrix(0.0,ncol=3,nrow=highb[1]*highb[2])
	for (i in lowb[1]:highb[1])
	{
		for(j in lowb[2]:highb[2])
		{
			index <- (i-1)*j+j
			if(is.na(llf(c(i,j),N,p,SF,dispfunc)))
			{
				print(c(i,j,N))
			}
			out[index,1]<-i
			out[index,2]<-j
			out[index,3]<-exp(llf(c(i,j),N,p,SF,dispfunc))
		}
	}
	return(out)
}

BICanti<-function(x,N,p,SF,rawvarfunc,llanti,llsense)
{
	nconds<-length(N)/2
	llvalueanti<-llanti(x,N,p,SF,rawvarfunc)
	llvaluesense<-llsense(N[1:nconds],N[(nconds+1):(2*nconds)],p,SF,rawvarfunc)
	d1 <- 2
	d2 <- 0
	n <- length(x)
	deltaBIC = 2*(llvalueanti-llvaluesense)-(d1-d2)*nconds*log(n)
	return(deltaBIC)
}



