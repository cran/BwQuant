bwRT <- function(x,y,tau){

  if(!is.numeric(tau)|!is.vector(tau)|any(!is.finite(tau))) stop("The parameter 'tau' must be a single number between 0 and 1")
  if(length(tau)!=1) stop("The parameter 'tau' must be a single number between 0 and 1")
  if(!is.numeric(tau)|(tau>=1)|(tau<=0)) stop("The parameter 'tau' must be a single number between 0 and 1")
  
  if (!is.numeric(x)){stop("'x' must be numeric")}
  if (!is.numeric(y)){stop("'y' must be numeric")}
  if(!length(x)>1 | !length(y)>1){stop("'x' and 'y' must be numeric vectors with the same length")}
  if(length(x)!=length(y)){stop("'x' and 'y' must be numeric vectors with the same length")}
  
  iaux=complete.cases(x,y)
  if(sum(!iaux)!=0){ii=which(iaux==FALSE);x=x[-ii];y=y[-ii];warning("Missing values have been removed from 'x' and 'y'")}

	n=length(x)
	N<-function(x,mu,sigma){exp(-(x-mu)^2/(2*sigma^2))/sqrt(2*sigma^2*pi)}
	RK<-integrate(function(x,mu,sigma){N(x,mu,sigma)^2},lower=-Inf,upper=Inf,mu=0,sigma=1)$value
	mu2=integrate(function(x){(x^2)*exp(-x^2/2)/sqrt(2*pi)},lower=-Inf,upper=Inf)$value

	ind=sort(x,ind=TRUE)$ix;x_ord=x[ind];y_ord=y[ind];dm=cbind(x_ord,x_ord^2,x_ord^3,x_ord^4)
	
	N_star=5; Nmax=max(min(n%/%20,N_star),1); RSS=numeric(Nmax)
	
	if(Nmax==1){Nb=1}else{
		RSSb=numeric(Nmax);ent=n%/%Nmax
		
		for(b in 1:(Nmax-1)){
			modelb<-try(rq(y_ord[(1+(b-1)*ent):(b*ent)]~dm[(1+(b-1)*ent):(b*ent),],tau=tau),TRUE)
			if(class(modelb)=="rq"){
				RSSb[b]=sum(tau*residuals(modelb)*(residuals(modelb)>=0)+(tau-1)*residuals(modelb)*(residuals(modelb)<0))
			}else{RSSb[b]=Inf}				
		}

		modelB=try(rq(y_ord[((Nmax-1)*ent+1):n]~dm[((Nmax-1)*ent+1):n,],tau=tau),TRUE)
		if(class(modelB)=="rq"){
				RSSb[Nmax]=sum(tau*residuals(modelB)*(residuals(modelB)>=0)+(tau-1)*residuals(modelB)*(residuals(modelB)<0))
		}else{RSSb[Nmax]=Inf}	
		
		RSS[Nmax]=sum(RSSb)
		
		modelb=rq(y_ord~dm,tau=tau)
		RSS[1]=sum(tau*residuals(modelb)*(residuals(modelb)>=0)+(tau-1)*residuals(modelb)*(residuals(modelb)<0))

		for(i in 2:(Nmax-1)){

			RSSb=numeric(i)
			for(b in 1:(i-1)){
				ent=n%/%i
				modelb<-try(rq(y_ord[(1+(b-1)*ent):(b*ent)]~dm[(1+(b-1)*ent):(b*ent),],tau=tau),TRUE)
				if(class(modelb)=="rq"){
					RSSb[b]=sum(tau*residuals(modelb)*(residuals(modelb)>=0)+(tau-1)*residuals(modelb)*(residuals(modelb)<0))
				}else{RSSb[b]=Inf}	
			}
			
			modelb=try(rq(y_ord[((i-1)*ent+1):n]~dm[((i-1)*ent+1):n,],tau=tau),TRUE)
			if(class(modelb)=="rq"){
				RSSb[i]=sum(tau*residuals(modelb)*(residuals(modelb)>=0)+(tau-1)*residuals(modelb)*(residuals(modelb)<0))
			}else{RSSb[i]=Inf}	
			
			RSS[i]=sum(RSSb)
		}
		
		Nmax=max((RSS!=Inf)*1:Nmax);Cp=numeric(Nmax)
		Cp[Nmax]=RSS[Nmax]/(RSS[Nmax]/(n-5*Nmax))-(n-10*Nmax)
		for(b in 1:(Nmax-1)){Cp[b]=RSS[b]/(RSS[Nmax]/(n-5*Nmax))-(n-10*b)}
		Nb=which.min(Cp)
	}

	ent=n%/%Nb;curvb1=numeric(Nb);sparsityb=numeric(Nb);lb=numeric(Nb);nb=numeric(Nb)

	for(b in 1:Nb){
		if(b==Nb){
			xb=dm[((Nb-1)*ent+1):n,];yb=y_ord[((Nb-1)*ent+1):n]
		}else{
			xb=dm[(1+(b-1)*ent):(b*ent),];yb=y_ord[(1+(b-1)*ent):(b*ent)]
		}
		
		nb[b]=dim(xb)[1];lb[b]=max(xb[,1])-min(xb[,1])
		modelb=rq(yb~xb,tau=tau); coefb=coef(modelb)		
		curvb1[b]=sum((2*coefb[3]+6*coefb[4]*xb[,1]+12*coefb[5]*xb[,2])^2)

		residuosb=as.numeric(sort(residuals(modelb)))
		hB=bwB(residuosb,tau)
		tau1=min(tau+hB,1);j1=ceiling(tau1*nb[b])
		tau2=max(tau-hB,0);j2=max(ceiling(tau2*nb[b]),1)
		sparsityb[b]=max(0.005,(residuosb[j1]-residuosb[j2]))/(2*hB)	
	}

	curv1=sum(curvb1)/n
	sparsity=sum(sparsityb^2*lb)
	h_pil=((RK*tau*(1-tau)*sparsity)/(n*(mu2^2)*curv1))^(1/5)
	
	return(h_pil)
}
