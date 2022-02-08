bwPI <- function(x,y,tau){

  if(!is.numeric(tau)|!is.vector(tau)|any(!is.finite(tau))) stop("The parameter 'tau' must be a single number between 0 and 1")
  if(length(tau)!=1) stop("The parameter 'tau' must be a single number between 0 and 1")
  if(!is.numeric(tau)|(tau>=1)|(tau<=0)) stop("The parameter 'tau' must be a single number between 0 and 1")
  
  if (!is.numeric(x)){stop("'x' must be numeric")}
  if (!is.numeric(y)){stop("'y' must be numeric")}
  if(!length(x)>1 | !length(y)>1){stop("'x' and 'y' must be numeric vectors with the same length")}
  if(length(x)!=length(y)){stop("'x' and 'y' must be numeric vectors with the same length")}
  
  iaux=complete.cases(x,y)
  if(sum(!iaux)!=0){ii=which(iaux==FALSE);x=x[-ii];y=y[-ii];warning("Missing values have been removed from 'x' and 'y'")}


	n=length(x);a=min(x);b=max(x);Rk=1/(2*sqrt(pi));mu=1
	point=seq(min(x),max(x),length=100);lp=length(point)
	x0=numeric(2);sparsity_p=numeric(length(point));q_second=numeric(n)
	
	alfa=0.05;c1=(1-alfa)*a+alfa*b;c2=alfa*a+(1-alfa)*b

	RuleThumb=RT(x,y,tau); g1=RuleThumb$g

	for(i in 1:n){
		if((c1<x[i])&(c2>x[i])){
			z1<-x-x[i];z2=z1^2;z3=z1^3
			aux<-dnorm(z1/g1);weight=aux/sum(aux)
			model=rq(y~z1+z2+z3,weights=weight,tau=tau);q_second[i]=2*coef(model)[3]
		}else{q_second[i]=0}

	}
	curvature=mean(q_second^2)

	x0[1]=max(0.0051,RuleThumb$h)
	x0[2]=n^(-0.2) * ((4.5 * dnorm(qnorm(tau))^4)/(2 * qnorm(tau)^2 + 1)^2)^0.2
	if(x0[2]>(1-tau)){x0[2]=0.99-tau}else if(tau<x0[2]){x0[2]=tau-0.01}
	IA=RuleThumb$IA;	IB=RuleThumb$IB;	IC=RuleThumb$IC;	IC2=RuleThumb$IC2; ID=RuleThumb$ID; IE=RuleThumb$IE

	Fun2 <- function(p){(IA/(n*p[1]*p[2])+IB*p[2]^2+IC*p[1]^2+IC2*p[1]^4)^2+ID/(n*p[2])+IE/(n^2*p[2]^2*p[1])}
	Amat <- matrix(c(0,0,0,1,-1,-1,1,0), nrow=4, ncol=2)
	bvec <- c((-tau+0.01), (tau-0.99),0.01, 0.005)
	band.pil<-constrOptim(x0, Fun2, NULL, ui = Amat, ci = bvec,control = list(reltol = 1e-100))$par

	y_estimate_sup<-numeric(lp);y_estimate_inf<-numeric(lp)
	beta0=1;beta1=0;z=0;weight=0
	storage.mode(n)="integer";storage.mode(z)="double";storage.mode(y)="double";storage.mode(weight)="double"
	storage.mode(beta0)="double";storage.mode(beta1)="double"
	tau_sup=min((tau+band.pil[2]),0.99);storage.mode(tau_sup)="double"
	tau_inf=max((tau-band.pil[2]),0.01);storage.mode(tau_inf)="double"
	for (i in 1:lp){
		z<-x-point[i]
		aux<-dnorm(z/band.pil[1])
		weight=aux/sum(aux)
		if((max(weight)>=0.9999)){
			y_estimate_sup[i]<-y[which.max(weight)]
			y_estimate_inf[i]<-y[which.max(weight)]
		}else{
			modelo_sup<-.Fortran("barro", n, z, y, weight, beta0, beta1,tau_sup)
			y_estimate_sup[i]<-modelo_sup[[5]]
			
			modelo_inf<-.Fortran("barro", n, z, y, weight, beta0, beta1,tau_inf)
			y_estimate_inf[i]<-modelo_inf[[5]]
		}
	}
		
	sparsity_p=(y_estimate_sup-y_estimate_inf)/(2*band.pil[2]); fp=(sparsity_p)^2; int_p=simpson(fp,a,b)
	h_PI=((Rk*tau*(1-tau)*int_p)/(n*mu^2*curvature))^(1/5)
	return(h_PI)
}
