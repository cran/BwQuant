bwCV <- function(x,y,hseq,tau){

  if(!is.numeric(tau)|!is.vector(tau)|any(!is.finite(tau))) stop("The parameter 'tau' must be a single number between 0 and 1")
  if(length(tau)!=1) stop("The parameter 'tau' must be a single number between 0 and 1")
  if(!is.numeric(tau)|(tau>=1)|(tau<=0)) stop("The parameter 'tau' must be a single number between 0 and 1")
  
  if (!is.numeric(x)){stop("'x' must be numeric")}
  if (!is.numeric(y)){stop("'y' must be numeric")}
  if(!length(x)>1 | !length(y)>1){stop("'x' and 'y' must be numeric vectors with the same length")}
  if(length(x)!=length(y)){stop("'x' and 'y' must be numeric vectors with the same length")}
  
  if(!is.numeric(hseq)|!is.vector(hseq)|any(!is.finite(hseq))) stop("'hseq' must be a sequence of values where the cross-validation function will be evaluated")
  if(!length(hseq)>1){stop("'hseq' must be a a sequence of values where the cross-validation function will be evaluated")}
  if(is.null(hseq)){ hseq=seq(sort(abs(outer(x,x,"-"))[outer(x,x,"-")!=0])[2],diff(range(x))/2,length=20)}
  
  iaux=complete.cases(x,y)
  if(sum(!iaux)!=0){ii=which(iaux==FALSE);x=x[-ii];y=y[-ii];warning("Missing values have been removed from 'x' and 'y'")}

	nh=length(hseq); n=length(y)-1;result=numeric(nh);beta0=1;beta1=0;zz=0;yy=0;weight=0
	storage.mode(n)="integer";storage.mode(zz)="double";storage.mode(yy)="double";storage.mode(weight)="double"
	storage.mode(beta0)="double";storage.mode(beta1)="double";storage.mode(tau)="double"

	fun_CV<-function(h){
		result=0
		for (i in 1:length(y)){ 
			zz<-x[-i]-x[i]
			aux<-dnorm(zz/h)
			if((sum(aux)!=0)&(max(aux/sum(aux))<0.9999)){
				weight=aux/sum(aux)					
				yy=y[-i]
				yi<-.Fortran("barro", n, zz, yy, weight, beta0, beta1,tau)[[5]]
			}else{yi<-Inf}
			result=result+(y[i]-yi)*(tau-((y[i]-yi)<0))
		}
		return(result)
	}
	return(hseq[which.min(lapply(hseq,fun_CV))])
}
