function [elem1]=tAbB_td_wtc(distv,cz)

   a    =   2.4795         ;    
   a_cc = a/sqrt(3.);
   c0 = 3.261 ;
   
   gamma0 = 2.7 ;
   gamma1 = 0.6602; %ABp  (tA&bB)

   q_sigma = c0 * 2.434857  ;
   q_pi = a_cc * 2.434857 ;


%dz = 0             % For in-plane hopping terms.

dx = distv(1)*a ;  % in A^o
dy = distv(2)*a ;
dz = cz*a ; 

   dis = sqrt( dx.^2 + dy.^2 + dz.^2 )  ; % in A^o
   
   WTC = true;
   
   if (WTC)
   v_pi = -gamma0 * exp(q_pi * (1.0 - dis / a_cc)) ;
   v_sigma = gamma1 * exp(q_sigma * (1.0 - dis / c0)) ;

   end

%     n =  cos(dz./ dis);  % rad
    n =  (dz./ dis);  

    elem1 =    ( v_pi.*(1 - n.^2  ) +   v_sigma.*n.^2 )  ;
end

