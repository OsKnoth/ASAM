      do 451 k=1,l
         do 451 j=1,mp                           
            do 451 i=1,np
               Rs=0.4
               rdc=sqrt((REAL(i-1)/REAL(n-1)-0.5)**2
     .              +(REAL(k-1)/REAL(l-1)-0.5)**2)/Rs

               ue(i,j,k)=u00
               ve(i,j,k)=v00
               the(i,j,k)=1./(rh00)

               if(rdc<=1) then
                  thetapol=atan2((REAL(k-1)/REAL(l-1)-0.5),
     .                 (REAL(i-1)/REAL(n-1)-0.5))

                  if(thetapol<0) thetapol=thetapol+2*pi


                  ue(i,j,k)=u00-1024.*sin(thetapol)*
     .                 (1.-rdc)**6.*rdc**6.
                  ve(i,j,k)=v00+1024.*cos(thetapol)*
     .                 (1.-rdc)**6.*rdc**6.

                  the(i,j,k)=1./(rh00+0.5*(1.-rdc**2.)**6.)

c     first term is offset from discretization, result of gcrk()
                  p(i,j,k)=p(i,j,k)+1024.**2.*(rdc**(36.)/72. - (6.*rdc**(35.))/35.
     .                 + (15.*rdc**(34.))/17. - (74.*rdc**(33.))/33.
     .                 + (57.*rdc**(32.))/32. + (174.*rdc**(31.))/31.
     .                 - (269.*rdc**(30.))/15. + (450.*rdc**(29.))/29.
     .                 + (153.*rdc**(28.))/8. - (1564.*rdc**(27.))/27.
     .                 + (510.*rdc**(26.))/13. + (204.*rdc**(25.))/5.
     .                 - (1473.*rdc**(24.))/16. + (1014.*rdc**(23.))/23.
     .                 + (1053.*rdc**(22.))/22. - (558.*rdc**(21.))/7.
     .                 + (783.*rdc**(20.))/20. + (54.*rdc**(19.))/19.
     .                 - (38.*rdc**(18.))/9. - (222.*rdc**(17.))/17.
     .                 + (609.*rdc**(16.))/32. - (184.*rdc**(15.))/15.
     .                 + (9.*rdc**(14.))/2. - (12.*rdc**(13.))/13.
     .                 + rdc**(12.)/12. - 34373./1805044411170.)*0.25
c     last term is value of integral at r=1

               endif
 451        continue
