function [elem]=G_kxky(kx,ky,tg)
a = 2.4795000553; % lattice constant in Ao                                     


g0 =  1;
g2 =  2.*cos(kx.*a)+4.*cos(ky.*sqrt(3).*a./2).*cos(kx.*a./2);
g5 =  exp(1j.*ky.*sqrt(3).*a)+2.*exp(-1j.*ky.*sqrt(3).*a./2).*cos(kx.*3.*a./2);

g5p = g5';

g6 = 2*cos(2*kx*a) + 4*cos(kx*a)*cos(ky*sqrt(3)*a) ;

g10 =  (2*cos(kx*2*a)*exp(1j*ky*sqrt(3)*a))+(2*cos(kx*a./2)*exp(-1j*ky*3*sqrt(3)*a./2))...
                      +(cos(kx*5*a./2)*exp(-1j*ky*sqrt(3)*a./2));
g10p = g10';
                      
g12 =  2*cos(3*kx*a) + 4*cos(3*kx*a/2)*cos(3*ky*sqrt(3)*a/2) ;

g15 =  exp(1j*ky*6*a./sqrt(3))+(2*cos(kx*3*a)*exp(-1j*ky*sqrt(3)*a)); 
g15p = g15';
g17 =  (2*cos(kx*7*a./2).*exp(1j*ky*sqrt(3)*a./2))+(2*cos(kx*a)*exp(-1j*ky*2*sqrt(3)*a))...
        +(2*cos(kx*5*a./2)*exp(-1j*ky*3*sqrt(3)*a./2));
g17p = g17'; 
   
     Gkxky = zeros(10,1);
    Gkxky(1) = g0;
    Gkxky(2) = g2;
    Gkxky(3) = g5;
    Gkxky(4) = g5p;
    Gkxky(5) = g6;
    Gkxky(6) = g10;
    Gkxky(7) = g10p;
    Gkxky(8) = g12;
    Gkxky(9) = g15;
    Gkxky(10) = g15p;
    Gkxky(11) = g17;
    Gkxky(12) = g17p;
    
elem = tg'*Gkxky  ;










