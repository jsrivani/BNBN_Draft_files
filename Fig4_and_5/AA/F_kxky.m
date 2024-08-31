function [elem]=F_kxky(kx,ky,tf)

a = 2.4795000553; % lattice constant in Ao                                     


f1 = exp(1j.*a.*ky./sqrt(3))+2.*exp(-1j.*a.*ky./(sqrt(3).*2)).*cos(a.*kx./2);
f3 =  exp(-1j.*ky.*2.*a./sqrt(3))+2.*exp(1j.*ky.*a./sqrt(3)).*cos(kx.*a);

f4 = 2*exp(1j*5*ky*a/(2*sqrt(3)))*cos(kx*a/2)+ 2*exp(-1j*ky*a/(2*sqrt(3)))*cos(3*kx*a/2)...
                                         + 2*exp(-1j*2*ky*a/sqrt(3))*cos(kx*a) ; 
                                     
f7 =  2*exp(1j*ky*a/sqrt(3))*cos(2*kx*a) +2*exp(1j*ky*5*a/(2*sqrt(3)))*cos(3*kx*a/2) + ...
         2*exp(-1j*ky*7*a/(2*sqrt(3)))*cos(kx*a/2); 
     
f8 =  exp(1j*ky*4*a/sqrt(3)) + 2*cos(2*kx*a)*exp(-1j*ky*2*a/(sqrt(3))); 

f9 = 2*exp(1j*4*ky*a/sqrt(3))*cos(kx*a) +2*exp(-1j*ky*a/(2*sqrt(3)))*cos(5*kx*a/2) + ...
                                 2*exp(-1j*ky*7*a/(2*sqrt(3)))*cos(3*kx*a/2);
                             
f11 =  exp(-1j*5*ky*a/sqrt(3)) + 2*exp(1j*5*ky*a/(2*sqrt(3)))*cos(5*kx*a/2) ;

f13 = 2*cos(2*kx*a)*exp(1j*ky*4*a/sqrt(3)) + ...
        2*cos(3*kx*a)*exp(1j*ky*a/sqrt(3)) + 2*cos(kx*a)*exp(-1j*ky*5*a/sqrt(3)) ; 
    
f14 = 2*cos(kx*a/2)*exp(1j*ky*11*a/(2*sqrt(3))) +2*cos(3*kx*a)*exp(-1j*ky*2*a/sqrt(3)) + ...
        2*cos(5*kx*a/2)*exp(-1j*7*ky*a/(2*sqrt(3))) ;
    
f16 =  2*cos(kx*7*a/2)*exp(-1j*ky*a/(2*sqrt(3))) +2*cos(kx*3*a/2)*exp(1j*ky*a*11/(2*sqrt(3))) + ...
        2*cos(kx*a*2)*exp(-1j*5*ky*a/sqrt(3)) ;
    
    Fkxky = zeros(10,1);
    Fkxky(1) = f1;
    Fkxky(2) = f3;
    Fkxky(3) = f4;
    Fkxky(4) = f7;
    Fkxky(5) = f8;
    Fkxky(6) = f9;
    Fkxky(7) = f11;
    Fkxky(8) = f13;
    Fkxky(9) = f14;
    Fkxky(10) = f16;
    

elem = tf'*Fkxky  ;