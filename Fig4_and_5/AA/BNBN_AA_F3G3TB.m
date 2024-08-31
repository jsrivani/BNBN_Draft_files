%                1     2
%            (A) C-----C (B)
%                '     '
%                '     '             -- AA-stack
%                '     '
%           (A') C-----C (B')
%                3     4      
%% ------------------------------------------------------------------------
clear; 
clc ;

nlayer = 2;
a = 2.4795000553; % lattice constant in Ao                                     
%% G functions
load BNBN_AA_pi1pi1.dat ;
load BNBN_AA_pi2pi2.dat ;
load BNBN_AA_pi3pi3.dat ;
load BNBN_AA_pi4pi4.dat ;
GAA = BNBN_AA_pi1pi1;
GBB = BNBN_AA_pi2pi2;
GApAp = BNBN_AA_pi3pi3;
GBpBp = BNBN_AA_pi4pi4;

load BNBN_AA_pi1pi3.dat;
load BNBN_AA_pi2pi4.dat;
GAAp = BNBN_AA_pi1pi3 ;
GBBp = BNBN_AA_pi2pi4 ;

%% F functions
load BNBN_AA_pi1pi2.dat;   
load BNBN_AA_pi3pi4.dat;
FAB = BNBN_AA_pi1pi2;  
FApBp = BNBN_AA_pi3pi4;

load BNBN_AA_pi1pi4.dat;
load BNBN_AA_pi2pi3.dat;
FABp = BNBN_AA_pi1pi4;
FApB = BNBN_AA_pi2pi3;
%%
ndim=4; % dimension of the hamiltonian                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
N =1000; 
%% Hopping energies

for i =1
tgaa = zeros(12,1); % AA
tgaa(1,1) = GAA(1,3);
tgaa(2,1) = GAA(2,3) ; 
tgaa(3,1) = GAA(8,3); 
tgaa(4,1) = GAA(8,3);
tgaa(5,1) = GAA(14,3);
tgaa(6,1) = GAA(20,3);     
tgaa(7,1) = GAA(20,3);      
tgaa(8,1) = GAA(32,3);       
tgaa(9,1) = GAA(38,3);      
tgaa(10,1) = GAA(38,3);      
tgaa(11,1) = GAA(44,3);       
tgaa(12,1) = GAA(44,3);  

tgbb = zeros(12,1); % BB
tgbb(1,1) = GBB(1,3);
tgbb(2,1) = GBB(2,3) ; 
tgbb(3,1) = GBB(8,3); 
tgbb(4,1) = GBB(8,3);
tgbb(5,1) = GBB(14,3);
tgbb(6,1) = GBB(20,3);     
tgbb(7,1) = GBB(20,3);      
tgbb(8,1) = GBB(32,3);       
tgbb(9,1) = GBB(38,3);      
tgbb(10,1) = GBB(38,3);      
tgbb(11,1) = GBB(44,3);       
tgbb(12,1) = GBB(44,3); 

tfab = zeros(10,1);  % AB
tfab(1,1) = FAB(1,3);                       
tfab(2,1) = FAB(4,3); 
tfab(3,1) = FAB(7,3);
tfab(4,1) = FAB(13,3);
tfab(5,1) = FAB(19,3); 
tfab(6,1) = FAB(22,3);        
tfab(7,1) = FAB(28,3);        
tfab(8,1) = FAB(31,3);       
tfab(9,1) = FAB(37,3);       
tfab(10,1) = FAB(43,3);

tgapap = zeros(12,1);  % ApAp
tgapap(1,1) = GApAp(1,3);
tgapap(2,1) = GApAp(2,3) ; 
tgapap(3,1) = GApAp(8,3); 
tgapap(4,1) = GApAp(8,3);
tgapap(5,1) = GApAp(14,3);
tgapap(6,1) = GApAp(20,3);     
tgapap(7,1) = GApAp(20,3);      
tgapap(8,1) = GApAp(32,3);       
tgapap(9,1) = GApAp(38,3);      
tgapap(10,1) = GApAp(38,3);      
tgapap(11,1) = GApAp(44,3);       
tgapap(12,1) = GApAp(44,3);  

tgbpbp = zeros(12,1);  % BpBp
tgbpbp(1,1) = GBpBp(1,3);       
tgbpbp(2,1) =  GBpBp(2,3);                
tgbpbp(3,1) = GBpBp(8,3) ;              
tgbpbp(4,1) = GBpBp(8,3);
tgbpbp(5,1) = GBpBp(14,3); 
tgbpbp(6,1) = GBpBp(20,3);       
tgbpbp(7,1) = GBpBp(20,3);      
tgbpbp(8,1) = GBpBp(32,3);        
tgbpbp(9,1) = GBpBp(38,3);       
tgbpbp(10,1) = GBpBp(38,3);      
tgbpbp(11,1) = GBpBp(44,3);       
tgbpbp(12,1) = GBpBp(44,3);      

tfbpap = zeros(10,1);   % ApBp
tfbpap(1,1) = FApBp(1,3);                       
tfbpap(2,1) = FApBp(4,3); 
tfbpap(3,1) = FApBp(7,3);
tfbpap(4,1) = FApBp(13,3);
tfbpap(5,1) = FApBp(19,3); 
tfbpap(6,1) = FApBp(22,3);        
tfbpap(7,1) = FApBp(28,3);        
tfbpap(8,1) = FApBp(31,3);       
tfbpap(9,1) = FApBp(37,3);       
tfbpap(10,1) = FApBp(43,3);


tgaap = zeros(12,1);   % AAp
tgaap(1,1) = GAAp(1,3);  
tgaap(2,1) = GAAp(2,3);  
tgaap(3,1) = GAAp(8,3) ;                
tgaap(4,1) = GAAp(11,3);
tgaap(5,1) = GAAp(14,3); 
tgaap(6,1) = GAAp(20,3);
tgaap(7,1) = GAAp(26,3);
tgaap(8,1) = GAAp(32,3); 
tgaap(9,1) = GAAp(38,3);
tgaap(10,1) = GAAp(41,3);
tgaap(11,1) = GAAp(44,3);
tgaap(12,1) = GAAp(50,4);

tgbbp = zeros(12,1);   % BBp
tgbbp(1,1) = GBBp(1,3);  
tgbbp(2,1) = GBBp(2,3);  
tgbbp(3,1) = GBBp(8,3) ;                
tgbbp(4,1) = GBBp(11,3);
tgbbp(5,1) = GBBp(14,3); 
tgbbp(6,1) = GBBp(20,3);
tgbbp(7,1) = GBBp(26,3);
tgbbp(8,1) = GBBp(32,3); 
tgbbp(9,1) = GBBp(38,3);
tgbbp(10,1) = GBBp(41,3);
tgbbp(11,1) = GBBp(44,3);
tgbbp(12,1) = GBBp(50,4);

tfbpa = zeros(10,1);  % ABp
tfbpa(1,1) = FABp(1,3);
tfbpa(2,1)=  FABp(4,3) ;
tfbpa(3,1)= FABp(7,3);  
tfbpa(4,1)= FABp(13,3);
tfbpa(5,1)= FABp(19,3) ;
tfbpa(6,1)= FABp(22,3);
tfbpa(7,1)= FABp(28,3); 
tfbpa(8,1)= FABp(31,3);
tfbpa(9,1)= FABp(37,3);
tfbpa(10,1)= FABp(43,3); 

tfapb = zeros(10,1);  % ApB
tfapb(1,1) = FApB(1,3);
tfapb(2,1)=  FApB(4,3) ;
tfapb(3,1)= FApB(7,3);  
tfapb(4,1)= FApB(13,3);
tfapb(5,1)= FApB(19,3) ;
tfapb(6,1)= FApB(22,3);
tfapb(7,1)= FApB(28,3); 
tfapb(8,1)= FApB(31,3);
tfapb(9,1)= FApB(37,3);
tfapb(10,1)= FApB(43,3); 
end

fdata(:,1) = tfab; %apbp
fdata(:,2) = tfbpa; %bap

dlmwrite('fdata_AA_bBN_jsv.dat',fdata)

gdata(:,1) = tgaa; %apap
gdata(:,2) = tgbb; %bpbp
gdata(:,3) = tgaap;
gdata(:,4) = tgbbp;

dlmwrite('gdata_AA_bBN_jsv.dat',gdata)

%% FULL_TB

hbb =@(kx,ky) G_kxky(kx,ky,tgbb);
hab =@(kx,ky) F_kxky(kx,ky,tfab);
haa =@(kx,ky) G_kxky(kx,ky,tgaa);

hbpbp =@(kx,ky) G_kxky(kx,ky,tgbpbp);
hapbp =@(kx,ky) F_kxky(kx,ky,tfbpap);
hapap =@(kx,ky) G_kxky(kx,ky,tgapap);  

haap =@(kx,ky) G_kxky(kx,ky,tgaap);
hbbp =@(kx,ky) G_kxky(kx,ky,tgbbp);
habp =@(kx,ky) F_kxky(kx,ky,tfbpa); 
hbap =@(kx,ky) conj(F_kxky(kx,ky,tfapb));

% Hamiltonian 
 H1 =@(kx,ky) [real(haa(kx,ky)) hab(kx,ky); hab(kx,ky)' real(hbb(kx,ky))];
 H2 =@(kx,ky) [ haap(kx,ky) (habp(kx,ky)); (hbap(kx,ky)) hbbp(kx,ky)];
 H3 =@(kx,ky) H2(kx,ky)';
 H4 =@(kx,ky) [real(hapap(kx,ky)) hapbp(kx,ky);hapbp(kx,ky)' real(hbpbp(kx,ky))];
 
 H =@(kx,ky) [H1(kx,ky) H2(kx,ky); H3(kx,ky) H4(kx,ky)];
 
%% F2G2_TB
%{
 %g-func
CAA0 = (GAA(1,3))-(3.*GAA(2,3))+ (6.*GAA(8,3))-(3.*GAA(14,3))-(6.*GAA(20,3))+(6.*GAA(32,3))+(6.*GAA(38,3))-(6.*GAA(44,3)); 
tgaa0 = GAA(1,3); 
tgaa2 = GAA(2,3); 
tgaa5n = (1/6).*(CAA0-tgaa0+(3*tgaa2));

CBB0 = (GBB(1,3))-(3.*GBB(2,3))+ (6.*GBB(8,3))-(3.*GBB(14,3))-(6.*GBB(20,3))+(6.*GBB(32,3))+(6.*GBB(38,3))-(6.*GBB(44,3)); 
tgbb0 = GBB(1,3); 
tgbb2 = GBB(2,3); 
tgbb5n = (1/6).*(CBB0-tgbb0+(3*tgbb2));

CApAp0 = (GApAp(1,3))-(3.*GApAp(2,3))+ (6.*GApAp(8,3))-(3.*GApAp(14,3))-(6.*GApAp(20,3))+(6.*GApAp(32,3))+(6.*GApAp(38,3))-(6.*GApAp(44,3)); 
tgapap0 = GApAp(1,3); 
tgapap2 = GApAp(2,3); 
tgapap5n = (1/6).*(CApAp0-tgapap0+(3*tgapap2));

CBpBp0 = (GBpBp(1,3))-(3.*GBpBp(2,3))+ (6.*GBpBp(8,3))-(3.*GBpBp(14,3))-(6.*GBpBp(20,3))+(6.*GBpBp(32,3))+(6.*GBpBp(38,3))-(6.*GBpBp(44,3)); 
tgbpbp0 = GBpBp(1,3); 
tgbpbp2 = GBpBp(2,3); 
tgbpbp5n = (1/6).*(CBpBp0-tgbpbp0+(3*tgbpbp2));

CAAp0 = (GAAp(1,3))-(3.*GAAp(2,3))+ (6.*GAAp(8,3))-(3.*GAAp(14,3))-(6.*GAAp(20,3))+(6.*GAAp(32,3))+(6.*GAAp(38,3))-(6.*GAAp(44,3)); 
tgaap0 = GAAp(1,3); 
tgaap2 = GAAp(2,3); 
tgaap5n = (1/6).*(CAAp0-tgaap0+(3*tgaap2));

CBBp0 = (GBBp(1,3))-(3.*GBBp(2,3))+ (6.*GBBp(8,3))-(3.*GBBp(14,3))-(6.*GBBp(20,3))+(6.*GBBp(32,3))+(6.*GBBp(38,3))-(6.*GBBp(44,3)); 
tgbbp0 = GBBp(1,3); 
tgbbp2 = GBBp(2,3); 
tgbbp5n = (1/6).*(CBBp0-tgbbp0+(3*tgbbp2));

%f-func

CAB1= (sqrt(3)*a./2).* (-FAB(1,3)+(2.*FAB(4,3)) +FAB(7,3)-(5.*FAB(13,3)) -(4.*FAB(19,3)) +(7.*FAB(22,3)) +(5.*FAB(28,3)) + (2.*FAB(31,3))-(4.*FAB(37,3)) + (11.*FAB(43,3)));
tfab1 = FAB(1,3);%-2*CAB1/(sqrt(3)*a);%
tfab3n = (CAB1./(sqrt(3)*a))+(tfab1./2);

CApBp1= (sqrt(3)*a./2).* (-FApBp(1,3)+(2.*FApBp(4,3)) +FApBp(7,3)-(5.*FApBp(13,3)) -(4.*FApBp(19,3)) +(7.*FApBp(22,3)) +(5.*FApBp(28,3)) + (2.*FApBp(31,3))-(4.*FApBp(37,3)) + (11.*FApBp(43,3)));
tfapbp1 = FApBp(1,3);%-2*CApBp1/(sqrt(3)*a);%
tfapbp3n = (CApBp1./(sqrt(3)*a))+(tfapbp1./2);

CABp1= (sqrt(3)*a./2).* (-FABp(1,3)+(2.*FABp(4,3)) +FABp(7,3)-(5.*FABp(13,3)) -(4.*FABp(19,3)) +(7.*FABp(22,3)) +(5.*FABp(28,3)) + (2.*FABp(31,3))-(4.*FABp(37,3)) + (11.*FABp(43,3)));
tfabp1 = FABp(1,3);%-2*CABp1/(sqrt(3)*a);%
tfabp3n = (CABp1./(sqrt(3)*a))+(tfabp1./2);

CApB1= (sqrt(3)*a./2).* (-FApB(1,3)+(2.*FApB(4,3)) +FApB(7,3)-(5.*FApB(13,3)) -(4.*FApB(19,3)) +(7.*FApB(22,3)) +(5.*FApB(28,3)) + (2.*FApB(31,3))-(4.*FApB(37,3)) + (11.*FApB(43,3)));
tfapb1 = FApB(1,3);%-2*CApB1/(sqrt(3)*a);%
tfapb3n = (CApB1./(sqrt(3)*a))+(tfapb1./2);

for i =1
tgaan = zeros(12,1); % AAx
tgaan(1,1) = tgaa0;
tgaan(2,1) = tgaa2 ; 
tgaan(3,1) = tgaa5n; 
tgaan(4,1) = tgaa5n;
 
tgbbn = zeros(12,1); % BB
tgbbn(1,1) = tgbb0;
tgbbn(2,1) = tgbb2 ; 
tgbbn(3,1) = tgbb5n; 
tgbbn(4,1) = tgbb5n;

tfabn = zeros(10,1);  % AB
tfabn(1,1) = tfab1;                       
tfabn(2,1) = tfab3n; 

tgapapn = zeros(12,1);  % ApAp
tgapapn(1,1) = tgapap0;
tgapapn(2,1) = tgapap2 ; 
tgapapn(3,1) = tgapap5n; 
tgapapn(4,1) = tgapap5n;

tgbpbpn = zeros(12,1);  % BpBp
tgbpbpn(1,1) = tgbpbp0;       
tgbpbpn(2,1) = tgbpbp2;                
tgbpbpn(3,1) = tgbpbp5n ;              
tgbpbpn(4,1) = tgbpbp5n;
     
tfbpapn = zeros(10,1);   % ApBp
tfbpapn(1,1) = tfapbp1;                       
tfbpapn(2,1) = tfapbp3n; 

tgaapn = zeros(12,1);   % AAp
tgaapn(1,1) = tgaap0;  
tgaapn(2,1) = tgaap2;  
tgaapn(3,1) = tgaap5n ;                
tgaapn(4,1) = tgaap5n;

tgbbpn = zeros(12,1);   % BBp
tgbbpn(1,1) = tgbbp0;  
tgbbpn(2,1) = tgbbp2;  
tgbbpn(3,1) = tgbbp5n ;                
tgbbpn(4,1) = tgbbp5n;
   
tfapbn = zeros(10,1);   % BAp
tfapbn(1,1) = tfapb1;
tfapbn(2,1) = tfapb3n;
   
tfbpan = zeros(10,1);  % ABp
tfbpan(1,1) = tfabp1;
tfbpan(2,1)=  tfabp3n;

end
%}
%% F3G3_TB
%{
 %g-func
CAA0 = (GAA(1,3))-(3.*GAA(2,3))+ (6.*GAA(8,3))-(3.*GAA(14,3))-(6.*GAA(20,3))+(6.*GAA(32,3))+(6.*GAA(38,3))-(6.*GAA(44,3)); 
tgaa0 = GAA(1,3); 
tgaa2 = GAA(2,3); 
tgaa5 = GAA(8,3) ;                
tgaa5p = GAA(11,3);
% tgaa5n = (1/6).*(CAA0-tgaa0+(3*tgaa2));
tgaa6n = 1./3.*(-CAA0+tgaa0-3.*tgaa2+6.*tgaa5);

CBB0 = (GBB(1,3))-(3.*GBB(2,3))+ (6.*GBB(8,3))-(3.*GBB(14,3))-(6.*GBB(20,3))+(6.*GBB(32,3))+(6.*GBB(38,3))-(6.*GBB(44,3)); 
tgbb0 = GBB(1,3); 
tgbb2 = GBB(2,3); 
tgbb5 = GBB(8,3) ;                
tgbb5p = GBB(11,3);
tgbb6n = 1./3.*(-CBB0+tgbb0-3.*tgbb2+6.*tgbb5);

CApAp0 = (GApAp(1,3))-(3.*GApAp(2,3))+ (6.*GApAp(8,3))-(3.*GApAp(14,3))-(6.*GApAp(20,3))+(6.*GApAp(32,3))+(6.*GApAp(38,3))-(6.*GApAp(44,3)); 
tgapap0 = GApAp(1,3); 
tgapap2 = GApAp(2,3); 
tgapap5 = GApAp(8,3) ;                
tgapap5p = GApAp(11,3);
tgapap6n = 1./3.*(-CApAp0+tgapap0-3.*tgapap2+6.*tgapap5);

CBpBp0 = (GBpBp(1,3))-(3.*GBpBp(2,3))+ (6.*GBpBp(8,3))-(3.*GBpBp(14,3))-(6.*GBpBp(20,3))+(6.*GBpBp(32,3))+(6.*GBpBp(38,3))-(6.*GBpBp(44,3)); 
tgbpbp0 = GBpBp(1,3); 
tgbpbp2 = GBpBp(2,3); 
tgbpbp5 = GBpBp(8,3) ;                
tgbpbp5p = GBpBp(11,3);
tgbpbp6n = 1./3.*(-CBpBp0+tgbpbp0-3.*tgbpbp2+6.*tgbpbp5);

CAAp0 = (GAAp(1,3))-(3.*GAAp(2,3))+ (6.*GAAp(8,3))-(3.*GAAp(14,3))-(6.*GAAp(20,3))+(6.*GAAp(32,3))+(6.*GAAp(38,3))-(6.*GAAp(44,3)); 
tgaap0 = GAAp(1,3); 
tgaap2 = GAAp(2,3); 
tgaap5 = GAAp(8,3);
tgaap5p = GAAp(11,3);
tgaap6n = 1./3.*(-CAAp0+tgaap0-3.*tgaap2+6.*tgaap5);

CBBp0 = (GBBp(1,3))-(3.*GBBp(2,3))+ (6.*GBBp(8,3))-(3.*GBBp(14,3))-(6.*GBBp(20,3))+(6.*GBBp(32,3))+(6.*GBBp(38,3))-(6.*GBBp(44,3)); 
tgbbp0 = GBBp(1,3); 
tgbbp2 = GBBp(2,3); 
tgbbp5 = GBBp(8,3); 
tgbbp5p = GBBp(11,3); 
tgbbp6n = 1./3.*(-CBBp0+tgbbp0-3.*tgbbp2+6.*tgbbp5);

%f-func
CAB1= (sqrt(3)*a./2).* (-FAB(1,3)+(2.*FAB(4,3)) +FAB(7,3)-(5.*FAB(13,3)) -(4.*FAB(19,3)) +(7.*FAB(22,3)) +(5.*FAB(28,3)) + (2.*FAB(31,3))-(4.*FAB(37,3)) + (11.*FAB(43,3)));
tfab1 = FAB(1,3);
tfab3 = FAB(4,3);
tfab4n = (2./(sqrt(3).*a)).*CAB1+tfab1-2.*tfab3;

CApBp1= (sqrt(3)*a./2).* (-FApBp(1,3)+(2.*FApBp(4,3)) +FApBp(7,3)-(5.*FApBp(13,3)) -(4.*FApBp(19,3)) +(7.*FApBp(22,3)) +(5.*FApBp(28,3)) + (2.*FApBp(31,3))-(4.*FApBp(37,3)) + (11.*FApBp(43,3)));
tfapbp1 = FApBp(1,3);
tfapbp3 = FApBp(4,3);
tfapbp4n = (2./(sqrt(3).*a)).*CApBp1+tfapbp1-2.*tfapbp3;

CABp1= (sqrt(3)*a./2).* (-FABp(1,3)+(2.*FABp(4,3)) +FABp(7,3)-(5.*FABp(13,3)) -(4.*FABp(19,3)) +(7.*FABp(22,3)) +(5.*FABp(28,3)) + (2.*FABp(31,3))-(4.*FABp(37,3)) + (11.*FABp(43,3)));
tfabp1 = FABp(1,3);
tfabp3 = FABp(4,3);
tfabp4n = (2./(sqrt(3).*a)).*CABp1+tfabp1-2.*tfabp3;

CApB1= (sqrt(3)*a./2).* (-FApB(1,3)+(2.*FApB(4,3)) +FApB(7,3)-(5.*FApB(13,3)) -(4.*FApB(19,3)) +(7.*FApB(22,3)) +(5.*FApB(28,3)) + (2.*FApB(31,3))-(4.*FApB(37,3)) + (11.*FApB(43,3)));
tfapb1 = FApB(1,3);
tfapb3 = FApB(4,3);
tfapb4n = (2./(sqrt(3).*a)).*CApB1+tfapb1-2.*tfapb3;

for i =1
tgaan = zeros(12,1); % AA
tgaan(1,1) = tgaa0;
tgaan(2,1) = tgaa2 ; 
tgaan(3,1) = tgaa5; 
tgaan(4,1) = tgaa5p;
tgaan(5,1) = tgaa6n;
 
tgbbn = zeros(12,1); % BB
tgbbn(1,1) = tgbb0;
tgbbn(2,1) = tgbb2 ; 
tgbbn(3,1) = tgbb5; 
tgbbn(4,1) = tgbb5p;
tgbbn(5,1) = tgbb6n;

tfabn = zeros(10,1);  % AB
tfabn(1,1) = tfab1;                       
tfabn(2,1) = tfab3; 
tfabn(3,1) = tfab4n; 

tgapapn = zeros(12,1);  % ApAp
tgapapn(1,1) = tgapap0;
tgapapn(2,1) = tgapap2 ; 
tgapapn(3,1) = tgapap5; 
tgapapn(4,1) = tgapap5p;
tgapapn(5,1) = tgapap6n;

tgbpbpn = zeros(12,1);  % BpBp
tgbpbpn(1,1) = tgbpbp0;       
tgbpbpn(2,1) = tgbpbp2;                
tgbpbpn(3,1) = tgbpbp5 ;              
tgbpbpn(4,1) = tgbpbp5p;
tgbpbpn(5,1) = tgbpbp6n;
     
tfbpapn = zeros(10,1);   % ApBp
tfbpapn(1,1) = tfapbp1;                       
tfbpapn(2,1) = tfapbp3; 
tfbpapn(3,1) = tfapbp4n; 

tgaapn = zeros(12,1);   % AAp
tgaapn(1,1) = tgaap0;  
tgaapn(2,1) = tgaap2;  
tgaapn(3,1) = tgaap5 ;                
tgaapn(4,1) = tgaap5p;
tgaapn(5,1) = tgaap6n;

tgbbpn = zeros(12,1);   % BBp
tgbbpn(1,1) = tgbbp0;  
tgbbpn(2,1) = tgbbp2;  
tgbbpn(3,1) = tgbbp5 ;                
tgbbpn(4,1) = tgbbp5p;
tgbbpn(5,1) = tgbbp6n;
   
tfapbn = zeros(10,1);   % BAp
tfapbn(1,1) = tfapb1;
tfapbn(2,1) = tfapb3;
tfapbn(3,1) = tfapb4n;
   
tfbpan = zeros(10,1);  % ABp
tfbpan(1,1) = tfabp1;
tfbpan(2,1)=  tfabp3;
tfbpan(3,1)=  tfabp4n;
end
%}
%% F4G4_TB
% {
 %g-func
CAA0 = (GAA(1,3))-(3.*GAA(2,3))+ (6.*GAA(8,3))-(3.*GAA(14,3))-(6.*GAA(20,3))+(6.*GAA(32,3))+(6.*GAA(38,3))-(6.*GAA(44,3)); 
tgaa0 = GAA(1,3); 
tgaa2 = GAA(2,3); 
tgaa5 = GAA(8,3) ;                
tgaa5p = GAA(11,3);
% tgaa5n = (1/6).*(CAA0-tgaa0+(3*tgaa2));
% tgaa6n = 1./3.*(-CAA0+tgaa0-3.*tgaa2+6.*tgaa5);
tgaa6 = GAA(14,3);
tgaa10n = (1/6).*(-CAA0+tgaa0-3.*tgaa2+6.*tgaa5-3*tgaa6);

CBB0 = (GBB(1,3))-(3.*GBB(2,3))+ (6.*GBB(8,3))-(3.*GBB(14,3))-(6.*GBB(20,3))+(6.*GBB(32,3))+(6.*GBB(38,3))-(6.*GBB(44,3)); 
tgbb0 = GBB(1,3); 
tgbb2 = GBB(2,3); 
tgbb5 = GBB(8,3) ;                
tgbb5p = GBB(11,3);
tgbb6 = GBB(14,3);
tgbb10n = (1/6).*(-CBB0+tgbb0-3.*tgbb2+6.*tgbb5-3*tgbb6);

CApAp0 = (GApAp(1,3))-(3.*GApAp(2,3))+ (6.*GApAp(8,3))-(3.*GApAp(14,3))-(6.*GApAp(20,3))+(6.*GApAp(32,3))+(6.*GApAp(38,3))-(6.*GApAp(44,3)); 
tgapap0 = GApAp(1,3); 
tgapap2 = GApAp(2,3); 
tgapap5 = GApAp(8,3) ;                
tgapap5p = GApAp(11,3);
tgapap6 = GApAp(14,3);
tgapap10n = (1/6).*(-CApAp0+tgapap0-3.*tgapap2+6.*tgapap5-3*tgapap6);

CBpBp0 = (GBpBp(1,3))-(3.*GBpBp(2,3))+ (6.*GBpBp(8,3))-(3.*GBpBp(14,3))-(6.*GBpBp(20,3))+(6.*GBpBp(32,3))+(6.*GBpBp(38,3))-(6.*GBpBp(44,3)); 
tgbpbp0 = GBpBp(1,3); 
tgbpbp2 = GBpBp(2,3); 
tgbpbp5 = GBpBp(8,3) ;                
tgbpbp5p = GBpBp(11,3);
tgbpbp6 = GBpBp(14,3);
tgbpbp10n = (1/6).*(-CBpBp0+tgbpbp0-3.*tgbpbp2+6.*tgbpbp5-3*tgbpbp6);

CAAp0 = (GAAp(1,3))-(3.*GAAp(2,3))+ (6.*GAAp(8,3))-(3.*GAAp(14,3))-(6.*GAAp(20,3))+(6.*GAAp(32,3))+(6.*GAAp(38,3))-(6.*GAAp(44,3)); 
tgaap0 = GAAp(1,3); 
tgaap2 = GAAp(2,3); 
tgaap5 = GAAp(8,3);
tgaap5p = GAAp(11,3);
tgaap6 = GAAp(14,3);
tgaap10n = (1/6).*(-CAAp0+tgaap0-3.*tgaap2+6.*tgaap5-3*tgaap6);

CBBp0 = (GBBp(1,3))-(3.*GBBp(2,3))+ (6.*GBBp(8,3))-(3.*GBBp(14,3))-(6.*GBBp(20,3))+(6.*GBBp(32,3))+(6.*GBBp(38,3))-(6.*GBBp(44,3)); 
tgbbp0 = GBBp(1,3); 
tgbbp2 = GBBp(2,3); 
tgbbp5 = GBBp(8,3); 
tgbbp5p = GBBp(11,3); 
tgbbp6 = GBBp(14,3);
tgbbp10n = (1/6).*(-CBBp0+tgbbp0-3.*tgbbp2+6.*tgbbp5-3*tgbbp6);

%f-func
CAB1= (sqrt(3)*a./2).* (-FAB(1,3)+(2.*FAB(4,3)) +FAB(7,3)-(5.*FAB(13,3)) -(4.*FAB(19,3)) +(7.*FAB(22,3)) +(5.*FAB(28,3)) + (2.*FAB(31,3))-(4.*FAB(37,3)) + (11.*FAB(43,3)));
tfab1 = FAB(1,3);
tfab3 = FAB(4,3);
% tfab3n = (CAB1./(sqrt(3)*a))+(tfab1./2);
% tfab4n = (2./(sqrt(3).*a)).*CAB1+tfab1-2.*tfab3;
tfab4 = FAB(7,3);
tfab7n = (-1/5).*((2*CAB1/(sqrt(3).*a))+tfab1-2.*tfab3-tfab4 );

CApBp1= (sqrt(3)*a./2).* (-FApBp(1,3)+(2.*FApBp(4,3)) +FApBp(7,3)-(5.*FApBp(13,3)) -(4.*FApBp(19,3)) +(7.*FApBp(22,3)) +(5.*FApBp(28,3)) + (2.*FApBp(31,3))-(4.*FApBp(37,3)) + (11.*FApBp(43,3)));
tfapbp1 = FApBp(1,3);
tfapbp3 = FApBp(4,3);
tfapbp4 = FApBp(7,3);
tfapbp7n = (-1/5).*((2*CApBp1/(sqrt(3).*a))+tfapbp1-2.*tfapbp3-tfapbp4 );

CABp1= (sqrt(3)*a./2).* (-FABp(1,3)+(2.*FABp(4,3)) +FABp(7,3)-(5.*FABp(13,3)) -(4.*FABp(19,3)) +(7.*FABp(22,3)) +(5.*FABp(28,3)) + (2.*FABp(31,3))-(4.*FABp(37,3)) + (11.*FABp(43,3)));
tfabp1 = FABp(1,3);
tfabp3 = FABp(4,3);
tfabp4 = FABp(7,3);
tfabp7n = (-1/5).*((2*CABp1/(sqrt(3).*a))+tfabp1-2.*tfabp3-tfabp4 );

CApB1= (sqrt(3)*a./2).* (-FApB(1,3)+(2.*FApB(4,3)) +FApB(7,3)-(5.*FApB(13,3)) -(4.*FApB(19,3)) +(7.*FApB(22,3)) +(5.*FApB(28,3)) + (2.*FApB(31,3))-(4.*FApB(37,3)) + (11.*FApB(43,3)));
tfapb1 = FApB(1,3);
tfapb3 = FApB(4,3);
tfapb4 = FApB(7,3);
tfapb7n = (-1/5).*((2*CApB1/(sqrt(3).*a))+tfapb1-2.*tfapb3-tfapb4 );

for i =1
tgaan = zeros(12,1); % AA
tgaan(1,1) = tgaa0;
tgaan(2,1) = tgaa2 ; 
tgaan(3,1) = tgaa5; 
tgaan(4,1) = tgaa5p;
tgaan(5,1) = tgaa6;
tgaan(6,1) = tgaa10n;
tgaan(7,1) = tgaa10n;

tgbbn = zeros(12,1); % BB
tgbbn(1,1) = tgbb0;
tgbbn(2,1) = tgbb2 ; 
tgbbn(3,1) = tgbb5; 
tgbbn(4,1) = tgbb5p;
tgbbn(5,1) = tgbb6;
tgbbn(6,1) = tgbb10n;
tgbbn(7,1) = tgbb10n;

tfabn = zeros(10,1);  % AB
tfabn(1,1) = tfab1;                       
tfabn(2,1) = tfab3; 
tfabn(3,1) = tfab4; 
tfabn(4,1) = tfab7n; 

tgapapn = zeros(12,1);  % ApAp
tgapapn(1,1) = tgapap0;
tgapapn(2,1) = tgapap2 ; 
tgapapn(3,1) = tgapap5; 
tgapapn(4,1) = tgapap5p;
tgapapn(5,1) = tgapap6;
tgapapn(6,1) = tgapap10n;
tgapapn(7,1) = tgapap10n;

tgbpbpn = zeros(12,1);  % BpBp
tgbpbpn(1,1) = tgbpbp0;       
tgbpbpn(2,1) = tgbpbp2;                
tgbpbpn(3,1) = tgbpbp5 ;              
tgbpbpn(4,1) = tgbpbp5p;
tgbpbpn(5,1) = tgbpbp6;
tgbpbpn(6,1) = tgbpbp10n;
tgbpbpn(7,1) = tgbpbp10n;

tfbpapn = zeros(10,1);   % ApBp
tfbpapn(1,1) = tfapbp1;                       
tfbpapn(2,1) = tfapbp3; 
tfbpapn(3,1) = tfapbp4; 
tfbpapn(4,1) = tfapbp7n;

tgaapn = zeros(12,1);   % AAp
tgaapn(1,1) = tgaap0;  
tgaapn(2,1) = tgaap2;  
tgaapn(3,1) = tgaap5 ;                
tgaapn(4,1) = tgaap5p;
tgaapn(5,1) = tgaap6;
tgaapn(6,1) = tgaap10n;
tgaapn(7,1) = tgaap10n;

tgbbpn = zeros(12,1);   % BBp
tgbbpn(1,1) = tgbbp0;  
tgbbpn(2,1) = tgbbp2;  
tgbbpn(3,1) = tgbbp5 ;                
tgbbpn(4,1) = tgbbp5p;
tgbbpn(5,1) = tgbbp6;
tgbbpn(6,1) = tgbbp10n;
tgbbpn(7,1) = tgbbp10n;
   
tfapbn = zeros(10,1);   % BAp
tfapbn(1,1) = tfapb1;
tfapbn(2,1) = tfapb3;
tfapbn(3,1) = tfapb4;
tfapbn(4,1) = tfapb7n;
   
tfbpan = zeros(10,1);  % ABp
tfbpan(1,1) = tfabp1;
tfbpan(2,1)=  tfabp3;
tfbpan(3,1)=  tfabp4;
tfbpan(4,1) = tfabp7n;

end

fdatan(:,1) = tfabn; 
fdatan(:,2) = tfbpapn;
dlmwrite('nn_fdata_AA_bBN_jsv.dat',fdatan)

gdatan(:,1) = tgaan; 
gdatan(:,2) = tgbbn; 
gdatan(:,3) = tgapapn;
gdatan(:,4) = tgbpbpn;
dlmwrite('nn_gdata_AA_bBN_jsv.dat',gdatan)

%}


hbbn =@(kx,ky) G_kxky(kx,ky,tgbbn);
habn =@(kx,ky) F_kxky(kx,ky,tfabn);
haan =@(kx,ky) G_kxky(kx,ky,tgaan);

hbpbpn =@(kx,ky) G_kxky(kx,ky,tgbpbpn);
hapbpn =@(kx,ky) F_kxky(kx,ky,tfbpapn);
hapapn =@(kx,ky) G_kxky(kx,ky,tgapapn);  

haapn =@(kx,ky) G_kxky(kx,ky,tgaapn);
hbbpn =@(kx,ky) G_kxky(kx,ky,tgbbpn);
habpn =@(kx,ky) F_kxky(kx,ky,tfbpan); 
hbapn =@(kx,ky) conj(F_kxky(kx,ky,tfapbn));

% Hamiltonian 
 H1n =@(kx,ky) [real(haan(kx,ky)) habn(kx,ky); habn(kx,ky)' real(hbbn(kx,ky))];
 H2n =@(kx,ky) [ haapn(kx,ky) (habpn(kx,ky)); (hbapn(kx,ky)) hbbpn(kx,ky)];
 H3n =@(kx,ky) H2n(kx,ky)';
 H4n =@(kx,ky) [real(hapapn(kx,ky)) hapbpn(kx,ky);hapbpn(kx,ky)' real(hbpbpn(kx,ky))];
 
 Hn =@(kx,ky) [H1n(kx,ky) H2n(kx,ky); H3n(kx,ky) H4n(kx,ky)];
 %}
%% defining the K-space

%Define special points of the BZ

Spoints = ((4.*pi)./(3.*a)).*[0 (3./4) 1 0;
                              0 (sqrt(3)./4) 0 0];
Kpoint = Spoints(1,3) ; 


% G -> M
  Gstart = Spoints(1,1);
  Gend = Kpoint.*(sqrt(3)./2);
  kxG2M = linspace(Spoints(1,1),Spoints(1,2),N);
  kyG2M = linspace(Spoints(2,1),Spoints(2,2),N);
  rG2M = linspace(Gstart,Gend,N);
  eG2M= zeros(N,ndim);
 
  for i = 1:N
     
      kx = kxG2M(i);
      ky = kyG2M(i);
                 
      [VC_G2M,engC_G2M] = eig(H(kx,ky));
      eG2M(i,:) = diag(engC_G2M);

      [VCn_G2M,engCn_G2M] = eig(Hn(kx,ky));
      eG2Mn(i,:) = diag(engCn_G2M);
  end 

% M -> K
 Mstart = Gend;
 Mend = (Kpoint.*(sqrt(3)+1))./2;
 kxM2K = linspace(Spoints(1,2),Spoints(1,3),N);
 kyM2K = linspace(Spoints(2,2),Spoints(2,3),N);
 rM2K = linspace(Mstart,Mend,N);
 eM2K= zeros(N,ndim);
  for i = 1:N
   
      kx = kxM2K(i);
      ky = kyM2K(i);
            
      [VC_M2K,engC_M2K] = eig(H(kx,ky));
      eM2K(i,:) = diag(engC_M2K);

      [VCn_M2K,engCn_M2K] = eig(Hn(kx,ky));
      eM2Kn(i,:) = diag(engCn_M2K);
         
  end
  
  % K -> G
  Kstart = Mend;
  Kend = Mend+Kpoint;
  kxK2G = linspace(Spoints(1,3),Spoints(1,4),N);
  kyK2G = linspace(Spoints(2,3),Spoints(2,4),N);
  rK2G = linspace(Kstart,Kend,N);
  eK2G = zeros(N,ndim);
  for i = 1:N
      
      kx = kxK2G(i);
      ky = kyK2G(i);
                  
      [VC_K2G,engC_K2G] = eig(H(kx,ky));
      eK2G(i,:) = diag(engC_K2G);
      
      [VCn_K2G,engCn_K2G] = eig(Hn(kx,ky));
      eK2Gn(i,:) = diag(engCn_K2G);
  end
   
  
%       [VC_G2M,engC_G2M] = eig(H(kxG2M(N),kyG2M(N)))
%       [VC_M2K,engC_M2K] = eig(H(kxM2K(N),kyM2K(N)))
%       [VC_K2G,engC_K2G] = eig(H(kxK2G(N),kyK2G(N)))
 
      [VCn_G2M,engCn_G2M] = eig(Hn(kxG2M(N),kyG2M(N)))    
      [VCn_M2K,engCn_M2K] = eig(Hn(kxM2K(N),kyM2K(N)))
      [VCn_K2G,engCn_K2G] = eig(Hn(kxK2G(N),kyK2G(N)))   


  %%
Ef = -0.0290;  
Efn = -0.0260;

figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
% DFT
load DFT_bands.dat
DFT = DFT_bands;

h1 = plot((DFT(:,1)*6.2827),DFT(:,2),'.k','Markersize',8);

%Wannier 90
load 'AA_band.dat'
Y = AA_band;

p = Y(:,1);
q = Y(:,2);
pp = p;
qq = q ;

h2 = plot(pp,qq,'.r','Markersize',5);
hold on

line([0 ; 4],[0 ; 0],'Linestyle',':','color','r','LineWidth',1) 
% Full-TB
x = 0.0;
 h3 = plot(rG2M+x,eG2M(:,1)-Ef,'b--','LineWidth',3);
  hold on
  plot(rG2M+x,eG2M(:,2)-Ef,'b--','LineWidth',3)
   hold on
  plot(rG2M+x,eG2M(:,3)-Ef,'b--','LineWidth',3)
   hold on
   plot(rG2M+x,eG2M(:,4)-Ef,'b--','LineWidth',3)
  hold on
  plot(rM2K+x,eM2K(:,1)-Ef,'b--','LineWidth',3)
   hold on
  plot(rM2K+x,eM2K(:,2)-Ef,'b--','LineWidth',3)
  plot(rM2K+x,eM2K(:,3)-Ef,'b--','LineWidth',3)
  plot(rM2K+x,eM2K(:,4)-Ef,'b--','LineWidth',3)
  plot(rK2G+x,eK2G(:,1)-Ef,'b--','LineWidth',3)
  plot(rK2G+x,eK2G(:,2)-Ef,'b--','LineWidth',3)
  plot(rK2G+x,eK2G(:,3)-Ef,'b--','LineWidth',3)
  plot(rK2G+x,eK2G(:,4)-Ef,'b--','LineWidth',3) 
  
 h4 = plot(rG2M+x,eG2Mn(:,1)-Efn,':','Color',[0.9290 0.6940 0.1250],'LineWidth',3);
  hold on
  plot(rG2M+x,eG2Mn(:,2)-Efn,':','Color',[0.9290 0.6940 0.1250],'LineWidth',3)
   hold on
  plot(rG2M+x,eG2Mn(:,3)-Efn,':','Color',[0.9290 0.6940 0.1250],'LineWidth',3)
   hold on
   plot(rG2M+x,eG2Mn(:,4)-Efn,':','Color',[0.9290 0.6940 0.1250],'LineWidth',3)
  hold on
  plot(rM2K+x,eM2Kn(:,1)-Efn,':','Color',[0.9290 0.6940 0.1250],'LineWidth',3)
   hold on
  plot(rM2K+x,eM2Kn(:,2)-Efn,':','Color',[0.9290 0.6940 0.1250],'LineWidth',3)
  plot(rM2K+x,eM2Kn(:,3)-Efn,':','Color',[0.9290 0.6940 0.1250],'LineWidth',3)
  plot(rM2K+x,eM2Kn(:,4)-Efn,':','Color',[0.9290 0.6940 0.1250],'LineWidth',3)
  plot(rK2G+x,eK2Gn(:,1)-Efn,':','Color',[0.9290 0.6940 0.1250],'LineWidth',3)
  plot(rK2G+x,eK2Gn(:,2)-Efn,':','Color',[0.9290 0.6940 0.1250],'LineWidth',3) 
  plot(rK2G+x,eK2Gn(:,3)-Efn,':','Color',[0.9290 0.6940 0.1250],'LineWidth',3) 
  plot(rK2G+x,eK2Gn(:,4)-Efn,':','Color',[0.9290 0.6940 0.1250],'LineWidth',3) 
  
legend([h1,h2,h3,h4],'DFT','W90','Full-TB','F4G4-TB')
ylabel('E (eV)');
xlabel('\Gamma-M-K-\Gamma');
title('AA-stacked BN/BN');
xlim(axes1,[0 4.029]);
box(axes1,'on');
set(axes1,'FontSize',24,'FontWeight','bold');
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.455137444489359 0.760942760942761 0.166425470332851 0.142255892255892]);
ylim([-10 10])
F1=550;F2=600;F3=550;F4=600;
figure1.Position=[F1 F2 F3 F4];

xlim([1.054  2.896])
ylim([-4.467 3.808])
%%
%{
% Ef = -0.1260;  
% Efn = (min(eM2Kn(:,3))+max(eM2Kn(:,2)))/2;  
Efn = -0.0260;


figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
% DFT
load DFT_bands.dat
DFT = DFT_bands;

h1 = plot((DFT(:,1)*6.2827),DFT(:,2),'.k','Markersize',8);

line([0 ; 4],[0 ; 0],'Linestyle',':','color','r','LineWidth',1) 

x = 0.0;
  
FG2 = [0.3010 0.7450 0.9330];
FG3 = [0.8500 0.3250 0.0980];
FG4 = [0.9290 0.6940 0.1250];


n_TB = FG2;

h4 = plot(rG2M+x,eG2Mn(:,1)-Efn,'Color',n_TB,'LineWidth',2);
  plot(rG2M+x,eG2Mn(:,2)-Efn,'Color',n_TB,'LineWidth',2)
  plot(rG2M+x,eG2Mn(:,3)-Efn,'Color',n_TB,'LineWidth',2)
  plot(rG2M+x,eG2Mn(:,4)-Efn,'Color',n_TB,'LineWidth',2)

  plot(rM2K+x,eM2Kn(:,1)-Efn,'Color',n_TB,'LineWidth',2)
  plot(rM2K+x,eM2Kn(:,2)-Efn,'Color',n_TB,'LineWidth',2)
  plot(rM2K+x,eM2Kn(:,3)-Efn,'Color',n_TB,'LineWidth',2)
  plot(rM2K+x,eM2Kn(:,4)-Efn,'Color',n_TB,'LineWidth',2)

  plot(rK2G+x,eK2Gn(:,1)-Efn,'Color',n_TB,'LineWidth',2)
  plot(rK2G+x,eK2Gn(:,2)-Efn,'Color',n_TB,'LineWidth',2) 
  plot(rK2G+x,eK2Gn(:,3)-Efn,'Color',n_TB,'LineWidth',2) 
  plot(rK2G+x,eK2Gn(:,4)-Efn,'Color',n_TB,'LineWidth',2) 
  
axes1.TickLabelInterpreter = 'latex';
axes1.FontWeight = 'bold';
set(gca,'XTick',[0 1.463 2.308 3.997])
set(gca,'XTickLabel',({'$\bf \Gamma$' , '$\bf M$' , '$\bf K$','$\bf \Gamma$','interpreter', 'latex'}))

axes1 = gca;
axes1.YAxis.TickLabelInterpreter = 'latex';
axes1.YAxis.TickLabelFormat      = '\\textbf{%g}';

leg1 = legend([h1,h4],'\bf DFT','\bf F2G2-TB');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',20);

ylabel('\bf E (eV)','Interpreter','latex');
title('\bf AA-stacked BN/BN','Interpreter','latex');
box(axes1,'on');
set(axes1,'FontSize',24,'FontWeight','bold');
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.455137444489359 0.760942760942761 0.166425470332851 0.142255892255892]);

F1=550;F2=600;F3=550;F4=600;
figure1.Position=[F1 F2 F3 F4];

xlim([0 3.997]);
ylim([-8.58 10.28])

 
xlim([1.054  2.896])
ylim([-4.467 3.808])
%}

%%
figure2 = figure;
axes2 = axes('Parent',figure2);
hold(axes2,'on');
% h5 =  plot(rG2M+x,(eG2M(:,2)-eG2Mn(:,2)),'-b','LineWidth',2)
h5 =  plot(rG2M+x,(eG2M(:,3)-Ef)-(eG2Mn(:,3)-Efn),'-r','LineWidth',2)
%   plot(rM2K+x,(eM2K(:,2)-eM2Kn(:,2)),'-b','LineWidth',2)
  plot(rM2K+x,(eM2K(:,3)-Ef)-(eM2Kn(:,3)-Efn),'-r','LineWidth',2)
%   plot(rK2G+x,(eK2G(:,2)-eK2Gn(:,2)),'-b','LineWidth',2) 
  plot(rK2G+x,(eK2G(:,3)-Ef)-(eK2Gn(:,3)-Efn),'-r','LineWidth',2) 
  
  line([0 ; 4],[0 ; 0],'Linestyle',':','color','r','LineWidth',1) 

axes2.TickLabelInterpreter = 'latex';
axes2.FontWeight = 'bold';
set(gca,'XTick',[0 1.463 2.308 3.997])
set(gca,'XTickLabel',({'$\bf \Gamma$' , '$\bf M$' , '$\bf K$','$\bf \Gamma$','interpreter', 'latex'}))

axes2 = gca;
axes2.YAxis.TickLabelInterpreter = 'latex';
axes2.YAxis.TickLabelFormat      = '\\textbf{%g}';

leg2 = legend(h5,'\bf (FTB-F2G2)');
set(leg2,'Interpreter','latex');
set(leg2,'FontSize',20);

ylabel('\bf E (eV)','Interpreter','latex');
% title('\bf AA-stacked BN/BN','Interpreter','latex');
box(axes2,'on');
set(axes2,'FontSize',24,'FontWeight','bold');

F1=550;F2=600;F3=550;F4=600;
figure2.Position=[F1 F2 F3 F4];

 
%%
%{

a  = 2.46;


CBBp0 = (0.4295)-(3.* 0.22349)+ (6.*0.03692)-(3.*-0.00253)-(6.*0.00076)+(6.* 0.00327)+(6.*-0.00085)-(6.*-0.00031 )
tgbbp0 = 0.4295
tgbbp2 =  0.22349
tgbbp5n = (1/6).*(CBBp0-tgbbp0+(3*tgbbp2))

%f-func
t1 = -3.010;
t2 = -0.2387;
t3 = 0.01900;
t4 = -0.01165;
t5 = -0.01167;
t6 = -0.00824;
t7 = 0.00386;
t8 = 0.00250;
t9 = 0.00224;
t10 = -0.00012;
CAB1= (sqrt(3)*a./2).* (-t1+(2.*t2) +t3-(5.*t4) -(4.*t5) +(7.*t6) +(5.*t7) + (2.*t8)-(4.*t9) + (11.*t10))
tfab1 = t1%-2*CAB1/(sqrt(3)*a)%FAB(1,3);%
% tfab3n = (CAB1./(sqrt(3)*a))+(tfab1./2)
tfab3n = (CAB1./(sqrt(3)*a))+(tfab1./2)


%%
5.567/(sqrt(3)*a)-3.010/2

CBBp0 = (GBBp(1,3))-(3.*GBBp(2,3))+ (6.*GBBp(8,3))-(3.*GBBp(14,3))-(6.*GBBp(20,3))+(6.*GBBp(32,3))+(6.*GBBp(38,3))-(6.*GBBp(44,3)); 
tgbbp0 = GBBp(1,3); 
tgbbp2 = GBBp(2,3); 
tgbbp5n = (1/6).*(CBBp0-tgbbp0+(3*tgbbp2));

%f-func

CAB1= (sqrt(3)*a./2).* (-FAB(1,3)+(2.*FAB(4,3)) +FAB(7,3)-(5.*FAB(13,3)) -(4.*FAB(19,3)) +(7.*FAB(22,3)) +(5.*FAB(28,3)) + (2.*FAB(31,3))-(4.*FAB(37,3)) + (11.*FAB(43,3)));
tfab1 = FAB(1,3);
tfab3n = (CAB1./(sqrt(3)*a))+(tfab1./2);
%}
% min(emod(:,ndim/2+1)-EF)-max(emod(:,ndim/2)-EF)
% 
% min(emod(:,ndim/2+1)-EF),max(emod(:,ndim/2)-EF)



%% 
clc

tgaan
tgbbn
tfabn

tgapapn
tgbpbpn
tfbpapn


tgaapn
tgbbpn
tfbpan
tfapbn


%{
tgaan =

    1.7666
    0.0053
    0.0223
    0.0223
   -0.0483
   -0.0007
   -0.0007
         0
         0
         0
         0
         0


tgbbn =

   -2.1843
    0.1923
    0.0195
    0.0195
   -0.0373
    0.0011
    0.0011
         0
         0
         0
         0
         0


tfabn =

   -2.7001
   -0.2102
    0.0797
   -0.0240
         0
         0
         0
         0
         0
         0


tgapapn =

    1.7667
    0.0052
    0.0223
    0.0223
   -0.0483
   -0.0007
   -0.0007
         0
         0
         0
         0
         0


tgbpbpn =

   -2.1840
    0.1923
    0.0195
    0.0195
   -0.0373
    0.0011
    0.0011
         0
         0
         0
         0
         0


tfbpapn =

   -2.7000
   -0.2102
    0.0797
   -0.0240
         0
         0
         0
         0
         0
         0


tgaapn =

    0.7270
    0.0498
    0.0089
    0.0089
    0.0139
   -0.0000
   -0.0000
         0
         0
         0
         0
         0


tgbbpn =

    0.2705
   -0.0185
   -0.0025
   -0.0025
    0.0070
    0.0001
    0.0001
         0
         0
         0
         0
         0


tfbpan =

    0.0265
    0.0082
   -0.0176
    0.0028
         0
         0
         0
         0
         0
         0


tfapbn =

    0.0265
    0.0082
   -0.0176
    0.0028
         0
         0
         0
         0
         0
         0

%}

%{
tgaap

tgaap =

    0.7270
    0.0498
    0.0089
    0.0089
    0.0139
    0.0017
    0.0017
    0.0012
    0.0003
    0.0003
   -0.0002
    0.0002

tgbbp

tgbbp =

    0.2705
   -0.0185
   -0.0025
   -0.0025
    0.0070
    0.0019
    0.0019
    0.0018
    0.0004
    0.0004
    0.0005
    0.0005

tfapb

tfapb =

    0.0265
    0.0082
   -0.0176
   -0.0045
   -0.0005
   -0.0046
   -0.0007
   -0.0009
   -0.0008
   -0.0003

tfbpa

tfbpa =

    0.0265
    0.0082
   -0.0176
   -0.0045
   -0.0005
   -0.0046
   -0.0007
   -0.0009
   -0.0008
   -0.0003
   %}













