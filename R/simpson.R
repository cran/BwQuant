simpson <- function(fxs,a,b){
	np<-length(fxs);h=(b-a)/(np-1)
	int<-3*(fxs[1]+fxs[np])/8+7*(fxs[2]+fxs[np-1])/6+23*(fxs[3]+fxs[np-2])/24+sum(fxs[4:(np-3)])
	return(int*h)
}
