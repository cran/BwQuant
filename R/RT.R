RT <- function(x,y,tau){

  if(!is.numeric(tau)|!is.vector(tau)|any(!is.finite(tau))) stop("The parameter 'tau' must be a single number between 0 and 1")
  if(length(tau)!=1) stop("The parameter 'tau' must be a single number between 0 and 1")
  if(!is.numeric(tau)|(tau>=1)|(tau<=0)) stop("The parameter 'tau' must be a single number between 0 and 1")
  
  if (!is.numeric(x)){stop("'x' must be numeric")}
  if (!is.numeric(y)){stop("'y' must be numeric")}
  if(!length(x)>1 | !length(y)>1){stop("'x' and 'y' must be numeric vectors with the same length")}
  if(length(x)!=length(y)){stop("'x' and 'y' must be numeric vectors with the same length")}
  
  iaux=complete.cases(x,y)
  if(sum(!iaux)!=0){ii=which(iaux==FALSE);x=x[-ii];y=y[-ii];warning("Missing values have been removed from 'x' and 'y'")}

  n=length(x);ran=max(x)-min(x)
  
  RK<-1/(2*sqrt(pi))
  Rd1K<-integrate(function(x){(-x*exp(-x^2/2)/sqrt(2*pi))^2},lower=-Inf,upper=Inf)$value	
  Rd2K<-integrate(function(x){((x^2-1)*exp(-x^2/2)/sqrt(2*pi))^2},lower=-Inf,upper=Inf)$value
  mu2=1 ; mu4=3 ; mu6=15
  RKc=integrate(function(x){(exp(-x^2/4)/sqrt(4*pi))^2},lower=-Inf,upper=Inf)$value
  
  IN2x4<-integrate(function(x){x^4*(exp(-x^2/2)/sqrt(2*pi))^2},lower=-Inf,upper=Inf)$value
  IN2x2<-integrate(function(x){x^2*(exp(-x^2/2)/sqrt(2*pi))^2},lower=-Inf,upper=Inf)$value
  
  d1N<-function (x, mu, sigma){-(exp(-((x - mu)^2/(2 * sigma^2))) * (x - mu)/(sigma^2 * sqrt(2 * (pi * sigma^2))))}
  d2N<-function (x, mu, sigma) {-((1 - (x - mu)^2/sigma^2) * exp(-((x - mu)^2/(2 * sigma^2)))/(sigma^2* sqrt(2 * (pi *sigma^2))))}
  d3N<-function(x,mu,sigma){2 * (((1 - (x - mu)^2/sigma^2)/2 + 1) * exp(-((x - mu)^2/(2 *sigma^2))) * (x - mu)/(sigma^4 * sqrt(2 * (pi *sigma^2))))}
  d4N<-function (x, mu, sigma){
    .e1 <- sigma^2 ; .e2 <- (x - mu)^2 ; .e4 <- (1 - .e2/.e1)/2 + 1
    2 * ((.e4 - 2 * ((.e4/2 + 0.5) * .e2/.e1)) * exp(-(.e2/(2 *.e1)))/(sigma^4 * sqrt(2 * (pi * .e1))))
  }
  
	D=mu2*mu4*mu6-mu4^3-mu2^3*mu6+mu2^2*mu4^2
	a31=-mu2^2*mu6+mu2*mu4^2;alfa31=a31/D
	a33=mu2*mu6-mu4^2;alfa33=a33/D
	delta1=(alfa31*mu4+alfa33*mu6)/6; delta2=4*(alfa31^2*RK+alfa33^2*IN2x4+2*alfa31*alfa33*IN2x2)
	C2I=((5*delta2)/(2*delta1))^(1/7);C2II=(delta2/delta1)^(1/7)


	ind=sort(x,ind=TRUE)$ix;x_ord=x[ind];y_ord=y[ind];dm=cbind(x_ord,x_ord^2,x_ord^3,x_ord^4)
	
	N_star=5;Nmax=max(min(n%/%20,N_star),1)
	RSS=numeric(Nmax)
	
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

	
	ent=n%/%Nb
	curvb1=numeric(Nb);curvb2=numeric(Nb);sparsityb=numeric(Nb)
	d2sb=numeric(Nb); Idtauq2=numeric(Nb);Iaux1=numeric(Nb);Iaux2=numeric(Nb);Iaux3=numeric(Nb)
	Ig=numeric(Nb);lb=numeric(Nb);nb=numeric(Nb)

	for(b in 1:Nb){
		if(b==Nb){
			xb=dm[((Nb-1)*ent+1):n,];yb=y_ord[((Nb-1)*ent+1):n]
		}else{
			xb=dm[(1+(b-1)*ent):(b*ent),];yb=y_ord[(1+(b-1)*ent):(b*ent)]
		}
		
		nb[b]=dim(xb)[1];lb[b]=max(xb[,1])-min(xb[,1])
		xseq=seq(min(xb[,1]),max(xb[,1]),by=0.01)
		modelb=rq(yb~xb,tau=tau); coefb=coef(modelb)
		
		curvb1[b]=sum((2*coefb[3]+6*coefb[4]*xb[,1]+12*coefb[5]*xb[,2])^2)
		curvb2[b]=sum((2*coefb[3]+6*coefb[4]*xb[,1]+12*coefb[5]*xb[,2])*24*coefb[5])
		q2<-function(x){2*coefb[3]+6*coefb[4]*x+12*coefb[5]*x^2}
		q2_2<-function(x){(2*coefb[3]+6*coefb[4]*x+12*coefb[5]*x^2)^2}

		residuosb=as.numeric(sort(residuals(modelb)))
		mu_est=mean(residuosb);sigma_est=sd(residuosb)

		hB=bwB(residuosb,tau)
		tau1=min(tau+hB,1);j1=ceiling(tau1*nb[b])
		tau2=max(tau-hB,0);j2=max(ceiling(tau2*nb[b]),1)
		sparsityb[b]=max(0.005,(residuosb[j1]-residuosb[j2]))/(2*hB)
		
		m=(0.25*nb[b]^(8/9))%/%1
		ntau=(nb[b]*tau)%/%1
		d2sb[b]=0.5*(nb[b]/m)^3*(residuosb[min(ntau+2*m,nb[b])]-2*residuosb[min(ntau+m,nb[b])]+2*residuosb[max(ntau-m,1)]-residuosb[max(ntau-2*m,1)])
	
		hBn=nb[b]^(-0.2) * ((4.5 * dnorm(qnorm(tau))^4)/(2 * qnorm(tau)^2 + 1)^2)^0.2
		tauplus=min(tau+hBn,0.99);tauminus=max(tau-hBn,0.01)
		modelbplus=rq(yb~xb,tau=tauplus);coefbplus=coef(modelbplus)
		modelbminus=rq(yb~xb,tau=tauminus);coefbminus=coef(modelbminus)
		dtauq2=function(x){(2*coefbplus[3]+6*coefbplus[4]*x+12*coefbplus[5]*x^2-2*coefbminus[3]-6*coefbminus[4]*x-12*coefbminus[5]*x^2)/(2*hBn)}
		Idtauq2[b]=integrate(dtauq2,min(xb[,1]),max(xb[,1]))$value

		
		Rf3<- integrate(function(x,mu,sigma){d3N(x,mu,sigma)*d3N(x,mu,sigma)},lower=min(xb),upper=max(xb),mu=mu_est,sigma=sigma_est)$value;Rf3
		Rf4<-integrate(function(x,mu,sigma){d4N(x,mu,sigma)*d4N(x,mu,sigma)},lower=min(xb),upper=max(xb),mu=mu_est,sigma=sigma_est)$value;Rf4
	
		h_pil_d1<-(3*Rd1K/(mu2^2*Rf3))^(1/7)*nb[b]^(-1/7); h_pil_d2<-(5*Rd2K/(mu2^2*Rf4))^(1/7)*nb[b]^(-1/9)

		prediction=predict(modelb)
		fd1<-function(x,h){sum(d1N((x-prediction)/h,mu=0,sigma=1))/(nb[b]*h^2)}
		fd2<-function(x,h){sum(d2N((x-prediction)/h,mu=0,sigma=1))/(nb[b]*h^3)}
	
		pred_new=as.numeric(cbind(rep(1,length(xseq)),xseq,xseq^2,xseq^3,xseq^4)%*%coefb)
		Iaux1[b]=simpson(sapply(pred_new,fd1,h=h_pil_d1)^2*q2_2(xseq),min(xb[,1]),max(xb[,1]))
		Iaux2[b]=simpson(sapply(pred_new,fd2,h=h_pil_d2)*q2_2(xseq),min(xb[,1]),max(xb[,1]))
		Iaux3[b]=simpson(sapply(pred_new,fd1,h=h_pil_d1)*q2(xseq)*dtauq2(xseq),min(xb[,1]),max(xb[,1]))
	
		Ig[b]=n*lb[b]/nb[b]
		
	}

	
	curv1=sum(curvb1)/n
	curv2=sum(curvb2)/n
	sparsity=sum(sparsityb^2*lb)
	
	IA=RK*sum(sparsityb^2*Ig*lb)/2
	IB=sum(sparsityb*d2sb*lb)/3
	IC=sum(sparsityb*Idtauq2)*mu2
	IC2=(sum(sparsityb^3*Iaux1)-sum(sparsityb^2*Iaux2)-sum(2*sparsityb*Iaux3))*mu4/4
	ID=sum(sparsityb^4*Ig*lb)*2
	IE=(0.5*RKc-RK)*sum(sparsityb^4*Ig^2*lb)

	if(curv2>0){g1=C2I*(((tau*(1-tau)*sparsity)/(curv2*n))^(1/7))}else{g1=C2II*(((tau*(1-tau)*sparsity)/(abs(curv2)*n))^(1/7))}
	h_pil=((RK*tau*(1-tau)*sparsity)/(n*(mu2^2)*curv1))^(1/5)
	
	return(list("h"=h_pil,"g"=g1,"IA"=IA,"IB"=IB,"IC"=IC,"IC2"=IC2,"ID"=ID,"IE"=IE))
}
