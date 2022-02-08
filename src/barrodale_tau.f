        subroutine barro(n,x,y,w,b1,b2,tau)
        implicit none
        integer n
*        parameter(nf=100)
        double precision x(n),y(n),w(n),wordy(n),yord(n),xordy(n)
        double precision c,mtot,med,xmed,b1,b2,tau
        integer indy(n),indf(n)
        integer i,corte,indc(2),auxint

        double precision v(n),z1(n),z2(n),coc(n),ww(n)
        double precision cv,cz1,cz2,pivot,prv,prz1,s
        double precision b1v,b1z1,b1z2,b2v,b2z1,b2z2
        integer indz2p(n),indcoc(n),auxi(n)
        integer r,rr,fu,pendiente
        
        double precision auxs,auxv(n),crit1,crit2
        integer criterio,iter, maxiter


*       Primera etapa del Algoritmo de Barrodale and Roberts

*       Se introduce la ordenada en el origen en la base

        do i=1,n
           indy(i)=i
        end do
        call quicksort(n,y,indy)
        do i=1,n
           yord(i)=y(indy(i))
           xordy(i)=x(indy(i))
        end do

        wordy=w(indy)

        mtot=sum(wordy)*tau
        c=0
        i=0
        do while ((i<n).and.(c<=mtot))
           i=i+1
           c=c+wordy(i)
        end do
        corte=i-1

        med=yord(corte+1)
        xmed=xordy(corte+1)
        cv=yord(corte+1)*wordy(corte+1)*tau
        cz2=xordy(corte+1)*wordy(corte+1)*tau
        do i=1,corte
           indf(i)=i
           v(i)=med-yord(i)
           cv=cv-yord(i)*wordy(i)*(1-tau)
           z1(i)=1
           z2(i)=xmed-xordy(i)
           cz2=cz2-xordy(i)*wordy(i)*(1-tau)
        end do
        do i=(corte+2),n
           auxint=i-1
           indf(auxint)=i
           v(auxint)=yord(i)-med
           cv=cv+yord(i)*wordy(i)*tau
           z1(auxint)=-1
           z2(auxint)=xordy(i)-xmed
           cz2=cz2+xordy(i)*wordy(i)*tau
        end do
        cz1=sum(wordy((corte+1):n))*tau-sum(wordy(1:corte))*(1-tau)

        indc(1)=corte+1
        b1v=med
        b1z1=1
        b1z2=xmed

        cv=cv-yord(corte+1)*cz1
        cz2=cz2-xordy(corte+1)*cz1
        cz1=-cz1


*       Se introduce la pendiente en la base

        if (.not.(cz2==0)) then
        
           if (cz2<0) then
              z2(1:(n-1))=-z2(1:(n-1))
              b1z2=-b1z2
              cz2=-cz2
              pendiente=-1
           else
              pendiente=1
           end if

           r=0
           do i=1,(n-1)
              if (z2(i)>0) then
                 r=r+1
                 indz2p(r)=i
                 coc(r)=v(i)/z2(i)
                 indcoc(r)=r
              end if
           end do
           rr=r
           call quicksort(rr,coc,indcoc)

           s=0
           r=0
           do while ((r<rr).and.(s<cz2))
              r=r+1
              auxi(r)=indz2p(indcoc(r))
              ww(r)=wordy(indf(auxi(r)))
              s=s+z2(auxi(r))*ww(r)
           end do
           fu=auxi(r)
           rr=r-1


           do r=1,rr
              auxint=auxi(r)
              cv=cv-v(auxint)*ww(r)
              cz1=cz1-z1(auxint)*ww(r)
              cz2=cz2-z2(auxint)*ww(r)

              v(auxint)=-v(auxint)
              z1(auxint)=-z1(auxint)
              z2(auxint)=-z2(auxint)
           end do

           pivot=z2(fu)
           prv=v(fu)/pivot
           prz1=z1(fu)/pivot

           indc(2)=indf(fu)
           do i=(fu+1),n
              auxint=i-1
              v(auxint)=v(i)
              z1(auxint)=z1(i)
              z2(auxint)=z2(i)
              indf(auxint)=indf(i)
           end do

           do i=1,(n-2)
              v(i)=v(i)-prv*z2(i)
              z1(i)=z1(i)-prz1*z2(i)
              z2(i)=-z2(i)/pivot
           end do

           b1v=b1v-prv*b1z2
           b1z1=b1z1-prz1*b1z2
           b1z2=-b1z2/pivot
           
           b2v=prv
           b2z1=prz1
           b2z2=1/pivot
           
           cv=cv-prv*cz2
           cz1=cz1-prz1*cz2
           cz2=-cz2/pivot

        end if


!        write(*,*) 'b1v, b2v:',b1v,b2v
!        write(*,*) 'cz1, cz2:',cz1,cz2

*       Segunda etapa del Algoritmo de Barrodale and Roberts

*Bucle donde se itera hasta que los costes marginales est‚n entre cotas

        crit1=max(cz1,-wordy(indc(1))-cz1)
        crit2=max(cz2,-wordy(indc(2))-cz2)

        if ((crit1>0).or.(crit2>0)) then
           criterio=1
        else
           criterio=0
        end if

        maxiter=10  ! N£mero m ximo de iteraciones

        iter=0
        do while (criterio==1)
           iter=iter+1
        ! Se escoge la columna que va a entrar en la base
        ! Si es la columna 1, se cambia con la 2
        
           if (crit1>crit2) then
                auxs=cz1
                cz1=cz2
                cz2=auxs

                auxv(1:(n-2))=z1(1:(n-2))
                z1(1:(n-2))=z2(1:(n-2))
                z2(1:(n-2))=auxv(1:(n-2))

                auxs=b1z1
                b1z1=b1z2
                b1z2=auxs

                auxs=b2z1
                b2z1=b2z2
                b2z2=auxs

                auxs=indc(2)
                indc(2)=indc(1)
                indc(1)=auxs
           end if

        ! Se mira si es la columna o su opuesta, la que se debe cambiar

           if (-wordy(indc(2))-cz2>cz2) then
                z2(1:(n-2))=-z2(1:(n-2))
                b1z2=-b1z2
                b2z2=-b2z2
                cz2=-wordy(indc(2))-cz2
           end if

           r=0
           do i=1,(n-2)
              if (z2(i)>0) then
                 r=r+1
                 indz2p(r)=i
                 coc(r)=v(i)/z2(i)
                 indcoc(r)=r
              end if
           end do
           rr=r
           call quicksort(rr,coc,indcoc)


           s=0
           r=0
           do while ((r<rr).and.(s<cz2))
              r=r+1
              auxi(r)=indz2p(indcoc(r))
              ww(r)=wordy(indf(auxi(r)))
              s=s+z2(auxi(r))*ww(r)
           end do
           fu=auxi(r)
           rr=r-1

           do r=1,rr
              auxint=auxi(r)
              cv=cv-v(auxint)*ww(r)
              cz1=cz1-z1(auxint)*ww(r)
              cz2=cz2-z2(auxint)*ww(r)

              v(auxint)=-v(auxint)
              z1(auxint)=-z1(auxint)
              z2(auxint)=-z2(auxint)
           end do

 !       write(*,*) 'iter,cv,cz1,cz2:',iter,cv,cz1,cz2


           pivot=z2(fu)
           prv=v(fu)/pivot
           prz1=z1(fu)/pivot

           b1v=b1v-prv*b1z2    ! Fila de la ordenada en el origen
           b1z1=b1z1-prz1*b1z2
           b1z2=-b1z2/pivot

           b2v=b2v-prv*b2z2   ! Fila de la pendiente
           b2z1=b2z1-prz1*b2z2
           b2z2=-b2z2/pivot

           v(fu)=prv   ! Fila que se cambia en la base
           z1(fu)=prz1
           z2(fu)=1/pivot

           do i=1,(fu-1)  ! Filas de los puntos muestrales que siguen en la base
              v(i)=v(i)-prv*z2(i)
              z1(i)=z1(i)-prz1*z2(i)
              z2(i)=-z2(i)/pivot
           end do
           do i=(fu+1),(n-2)  ! Filas de los puntos muestrales que siguen en la base
              v(i)=v(i)-prv*z2(i)
              z1(i)=z1(i)-prz1*z2(i)
              z2(i)=-z2(i)/pivot
           end do

           cv=cv-prv*cz2   ! Fila de los costes marginales
           cz1=cz1-prz1*cz2
           cz2=-cz2/pivot

           auxint=indf(fu)
           indf(fu)=indc(2)
           indc(2)=auxint

           crit1=max(cz1,-wordy(indc(1))-cz1)
           crit2=max(cz2,-wordy(indc(2))-cz2)

           if ((max(crit1,crit2)<=0).or.(iter>maxiter)) then
              criterio=0
           end if

!           write(*,*) 'iter, cz1, cz2:',iter,cz1,cz2
!           write(*,*) 'b1v, b2v:',b1v,b2v
!           write(*,*) 'criterio:',criterio
!           write(*,*) 'crit1, crit2:',crit1,crit2

        end do

        b1=b1v
        if (pendiente==1) then
            b2=b2v
        else
            b2=-b2v
        end if

        end subroutine barro




        subroutine quicksort(n,x,ind)
        implicit none
        integer n,nf
!        parameter(nf=100)
        double precision x(n)
        integer i,ic,nc,newnc,ns,nl,iit,fft,e,iter,crit
        integer ind(n),it(n),ft(n),newit(n),newft(n),inf(n),sup(n)
        double precision chosen

*        n=size(x)
        crit=0
        iter=0
        nc=1
        it(1)=1
        ft(1)=n
*        do i=1,n
*           ind(i)=i
*        end do
        do while (crit==0)
           iter=iter+1
           newnc=0
*           write(*,*) iter,nc
           do ic=1,nc
                 iit=it(ic)
                 fft=ft(ic)
                 e=ind(iit)
                 chosen=x(ind(iit))
                 ns=0
                 nl=0

*                B£squeda de ¡ndices por debajo y por arriba
                 do i=(iit+1),fft
                    if (x(ind(i))<chosen) then
                       ns=ns+1
                       inf(ns)=ind(i)
                    else
                       nl=nl+1
                       sup(nl)=ind(i)
                    end if
                 end do

*                Recolocaci¢n de los ¡ndices obtenidos
                 if (ns>0) then
                    ind(iit:iit+ns-1)=inf(1:ns)
                    if (ns>1) then
                       newnc=newnc+1
                       newit(newnc)=iit
                       newft(newnc)=iit+ns-1
                    end if
                 end if

                 ind(iit+ns)=e

                 if (nl>0) then
                    ind(iit+ns+1:fft)=sup(1:nl)
                    if (nl>1) then
                       newnc=newnc+1
                       newit(newnc)=iit+ns+1
                       newft(newnc)=fft
                    end if
                 end if

           end do
           if ((newnc==0).or.(iter>n)) then
              crit=1
           else
              it(1:newnc)=newit(1:newnc)
              ft(1:newnc)=newft(1:newnc)
              nc=newnc
           end if


        end do
           
        end subroutine quicksort

