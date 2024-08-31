%                       
%                1     2
%            (A) B-----N (B)
%                '     
%                '                  -- BA-stack
%                '     
%     (A') B-----N (B')
%          3     4         
%% ------------------------------------------------------------------------
clear; 
clc ;

nlayer = 2;
a = 2.4795000553; % lattice constant in Ao                                     
%% G functions
load BNNB_BA_pi1pi1.dat ;
load BNNB_BA_pi2pi2.dat ;
load BNNB_BA_pi3pi3.dat ;
load BNNB_BA_pi4pi4.dat ;
GAA = BNNB_BA_pi1pi1;
GBB = BNNB_BA_pi2pi2;
GApAp = BNNB_BA_pi3pi3;
GBpBp = BNNB_BA_pi4pi4;

load BNNB_BA_pi1pi3.dat;
load BNNB_BA_pi2pi4.dat;
FAAp = BNNB_BA_pi1pi3 ;
FBBp = BNNB_BA_pi2pi4 ;

%% F functions
load BNNB_BA_pi1pi2.dat;   
load BNNB_BA_pi3pi4.dat;
load BNNB_BA_pi2pi3.dat;
FAB = BNNB_BA_pi1pi2;  
FApBp = BNNB_BA_pi3pi4;
FApB = BNNB_BA_pi2pi3;

load BNNB_BA_pi1pi4.dat;
GABp = BNNB_BA_pi1pi4;
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


tgabp = zeros(12,1);   % ABp
tgabp(1,1) = GABp(1,3);  
tgabp(2,1) = GABp(2,3);  
tgabp(3,1) = GABp(8,3) ;                
tgabp(4,1) = GABp(11,3);
tgabp(5,1) = GABp(14,3); 
tgabp(6,1) = GABp(20,3);
tgabp(7,1) = GABp(26,3);
tgabp(8,1) = GABp(32,3); 
tgabp(9,1) = GABp(38,3);
tgabp(10,1) = GABp(41,3);
tgabp(11,1) = GABp(44,3);
tgabp(12,1) = GABp(50,4);

tfbpb = zeros(10,1);   % BBp
tfbpb(1,1) = FBBp(1,3);
tfbpb(2,1) = FBBp(4,3);
tfbpb(3,1) = FBBp(7,3); 
tfbpb(4,1) = FBBp(13,3); 
tfbpb(5,1) = FBBp(19,3);  
tfbpb(6,1) = FBBp(22,3);          
tfbpb(7,1) = FBBp(28,3);         
tfbpb(8,1) = FBBp(31,3);         
tfbpb(9,1) = FBBp(37,3);         
tfbpb(10,1) = FBBp(43,3);    

tfapa = zeros(10,1);   % AAp
tfapa(1,1) = FAAp(1,3);
tfapa(2,1) = FAAp(4,3);
tfapa(3,1) = FAAp(7,3); 
tfapa(4,1) = FAAp(13,3); 
tfapa(5,1) = FAAp(19,3);  
tfapa(6,1) = FAAp(22,3);          
tfapa(7,1) = FAAp(28,3);         
tfapa(8,1) = FAAp(31,3);         
tfapa(9,1) = FAAp(37,3);         
tfapa(10,1) = FAAp(43,3);   

tfbap = zeros(10,1);  % ApB
tfbap(1,1) = FApB(1,3);
tfbap(2,1)=  FApB(4,3) ;
tfbap(3,1)= FApB(7,3);  
tfbap(4,1)= FApB(13,3);
tfbap(5,1)= FApB(19,3) ;
tfbap(6,1)= FApB(22,3);
tfbap(7,1)= FApB(28,3); 
tfbap(8,1)= FApB(31,3);
tfbap(9,1)= FApB(37,3);
tfbap(10,1)= FApB(43,3); 
end

fdata(:,1) = tfab; 
fdata(:,2) = tfbpap;
dlmwrite('fdata_BAp_bBN_jsv.dat',fdata)

gdata(:,1) = tgaa; 
gdata(:,2) = tgbb; 
gdata(:,3) = tgapap;
gdata(:,4) = tgbpbp;
dlmwrite('gdata_BAp_bBN_jsv.dat',gdata)

%% FULL_TB

hbb =@(kx,ky) G_kxky(kx,ky,tgbb);
hab =@(kx,ky) F_kxky(kx,ky,tfab);
haa =@(kx,ky) G_kxky(kx,ky,tgaa);

hbpbp =@(kx,ky) G_kxky(kx,ky,tgbpbp);
hapbp =@(kx,ky) F_kxky(kx,ky,tfbpap);
hapap =@(kx,ky) G_kxky(kx,ky,tgapap);  

hbbp =@(kx,ky) conj(F_kxky(kx,ky,tfbpb));
haap =@(kx,ky) conj(F_kxky(kx,ky,tfapa));
hbap =@(kx,ky) (F_kxky(kx,ky,tfbap));
habp =@(kx,ky) G_kxky(kx,ky,tgabp); 

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

CABp0 = (GABp(1,3))-(3.*GABp(2,3))+ (6.*GABp(8,3))-(3.*GABp(14,3))-(6.*GABp(20,3))+(6.*GABp(32,3))+(6.*GABp(38,3))-(6.*GABp(44,3)); 
tgabp0 = GABp(1,3); 
tgabp2 = GABp(2,3); 
tgabp5n = (1/6).*(CABp0-tgabp0+(3*tgabp2));

%f-func

CAB1= (sqrt(3)*a./2).* (-FAB(1,3)+(2.*FAB(4,3)) +FAB(7,3)-(5.*FAB(13,3)) -(4.*FAB(19,3)) +(7.*FAB(22,3)) +(5.*FAB(28,3)) + (2.*FAB(31,3))-(4.*FAB(37,3)) + (11.*FAB(43,3)));
tfab1 = FAB(1,3);
tfab3n = (CAB1./(sqrt(3)*a))+(tfab1./2);

CApBp1= (sqrt(3)*a./2).* (-FApBp(1,3)+(2.*FApBp(4,3)) +FApBp(7,3)-(5.*FApBp(13,3)) -(4.*FApBp(19,3)) +(7.*FApBp(22,3)) +(5.*FApBp(28,3)) + (2.*FApBp(31,3))-(4.*FApBp(37,3)) + (11.*FApBp(43,3)));
tfapbp1 = FApBp(1,3);
tfapbp3n = (CApBp1./(sqrt(3)*a))+(tfapbp1./2);

CApB1= (sqrt(3)*a./2).* (-FApB(1,3)+(2.*FApB(4,3)) +FApB(7,3)-(5.*FApB(13,3)) -(4.*FApB(19,3)) +(7.*FApB(22,3)) +(5.*FApB(28,3)) + (2.*FApB(31,3))-(4.*FApB(37,3)) + (11.*FApB(43,3)));
tfapb1 = FApB(1,3);
tfapb3n = (CApB1./(sqrt(3)*a))+(tfapb1./2);

CAAp1= (sqrt(3)*a./2).* (-FAAp(1,3)+(2.*FAAp(4,3)) +FAAp(7,3)-(5.*FAAp(13,3)) -(4.*FAAp(19,3)) +(7.*FAAp(22,3)) +(5.*FAAp(28,3)) + (2.*FAAp(31,3))-(4.*FAAp(37,3)) + (11.*FAAp(43,3)));
tfaap1 = FAAp(1,3);
tfaap3n = (CAAp1./(sqrt(3)*a))+(tfaap1./2);


CBBp1= (sqrt(3)*a./2).* (-FBBp(1,3)+(2.*FBBp(4,3)) +FBBp(7,3)-(5.*FBBp(13,3)) -(4.*FBBp(19,3)) +(7.*FBBp(22,3)) +(5.*FBBp(28,3)) + (2.*FBBp(31,3))-(4.*FBBp(37,3)) + (11.*FBBp(43,3)));
tfbbp1 = FBBp(1,3);
tfbbp3n = (CBBp1./(sqrt(3)*a))+(tfbbp1./2);

for i =1
tgaan = zeros(12,1); % AA
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

tgabpn = zeros(12,1);   % ABp
tgabpn(1,1) = tgabp0;  
tgabpn(2,1) = tgabp2;  
tgabpn(3,1) = tgabp5n ;                
tgabpn(4,1) = tgabp5n;

tfbpbn = zeros(10,1);   % BBp
tfbpbn(1,1) = tfbbp1;
tfbpbn(2,1) = tfbbp3n;
   
tfapan = zeros(10,1);   % AAp
tfapan(1,1) = tfaap1;
tfapan(2,1) = tfaap3n;
   
tfapbn = zeros(10,1);  % ApB
tfapbn(1,1) = tfapb1;
tfapbn(2,1)=  tfapb3n;

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

CABp0 = (GABp(1,3))-(3.*GABp(2,3))+ (6.*GABp(8,3))-(3.*GABp(14,3))-(6.*GABp(20,3))+(6.*GABp(32,3))+(6.*GABp(38,3))-(6.*GABp(44,3)); 
tgabp0 = GABp(1,3); 
tgabp2 = GABp(2,3); 
tgabp5 = GABp(8,3) ;                
tgabp5p = GABp(11,3);
tgabp6n = 1./3.*(-CABp0+tgabp0-3.*tgabp2+6.*tgabp5);

%f-func
CAB1= (sqrt(3)*a./2).* (-FAB(1,3)+(2.*FAB(4,3)) +FAB(7,3)-(5.*FAB(13,3)) -(4.*FAB(19,3)) +(7.*FAB(22,3)) +(5.*FAB(28,3)) + (2.*FAB(31,3))-(4.*FAB(37,3)) + (11.*FAB(43,3)));
tfab1 = FAB(1,3);
% tfab3n = (CAB1./(sqrt(3)*a))+(tfab1./2);
tfab3 = FAB(4,3);
tfab4n = (2./(sqrt(3).*a)).*CAB1+tfab1-2.*tfab3;

CApBp1= (sqrt(3)*a./2).* (-FApBp(1,3)+(2.*FApBp(4,3)) +FApBp(7,3)-(5.*FApBp(13,3)) -(4.*FApBp(19,3)) +(7.*FApBp(22,3)) +(5.*FApBp(28,3)) + (2.*FApBp(31,3))-(4.*FApBp(37,3)) + (11.*FApBp(43,3)));
tfapbp1 = FApBp(1,3);
tfapbp3 = FApBp(4,3);
tfapbp4n = (2./(sqrt(3).*a)).*CApBp1+tfapbp1-2.*tfapbp3;

CApB1= (sqrt(3)*a./2).* (-FApB(1,3)+(2.*FApB(4,3)) +FApB(7,3)-(5.*FApB(13,3)) -(4.*FApB(19,3)) +(7.*FApB(22,3)) +(5.*FApB(28,3)) + (2.*FApB(31,3))-(4.*FApB(37,3)) + (11.*FApB(43,3)));
tfapb1 = FApB(1,3);
tfapb3 = FApB(4,3);
tfapb4n = (2./(sqrt(3).*a)).*CApB1+tfapb1-2.*tfapb3;

CAAp1= (sqrt(3)*a./2).* (-FAAp(1,3)+(2.*FAAp(4,3)) +FAAp(7,3)-(5.*FAAp(13,3)) -(4.*FAAp(19,3)) +(7.*FAAp(22,3)) +(5.*FAAp(28,3)) + (2.*FAAp(31,3))-(4.*FAAp(37,3)) + (11.*FAAp(43,3)));
tfaap1 = FAAp(1,3);
tfaap3 = FAAp(4,3);
tfaap4n = (2./(sqrt(3).*a)).*CAAp1+tfaap1-2.*tfaap3;


CBBp1= (sqrt(3)*a./2).* (-FBBp(1,3)+(2.*FBBp(4,3)) +FBBp(7,3)-(5.*FBBp(13,3)) -(4.*FBBp(19,3)) +(7.*FBBp(22,3)) +(5.*FBBp(28,3)) + (2.*FBBp(31,3))-(4.*FBBp(37,3)) + (11.*FBBp(43,3)));
tfbbp1 = FBBp(1,3);
tfbbp3 = FBBp(4,3);
tfbbp4n = (2./(sqrt(3).*a)).*CBBp1+tfbbp1-2.*tfbbp3;

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

tgabpn = zeros(12,1);   % ABp
tgabpn(1,1) = tgabp0;  
tgabpn(2,1) = tgabp2;  
tgabpn(3,1) = tgabp5 ;                
tgabpn(4,1) = tgabp5p;
tgabpn(5,1) = tgabp6n;

tfbpbn = zeros(10,1);   % BBp
tfbpbn(1,1) = tfbbp1;
tfbpbn(2,1) = tfbbp3;
tfbpbn(3,1) = tfbbp4n;
   
tfapan = zeros(10,1);   % AAp
tfapan(1,1) = tfaap1;
tfapan(2,1) = tfaap3;
tfapan(3,1) = tfaap4n;
      
tfapbn = zeros(10,1);  % ApB
tfapbn(1,1) = tfapb1;
tfapbn(2,1)=  tfapb3;
tfapbn(3,1)=  tfapb4n;

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

CABp0 = (GABp(1,3))-(3.*GABp(2,3))+ (6.*GABp(8,3))-(3.*GABp(14,3))-(6.*GABp(20,3))+(6.*GABp(32,3))+(6.*GABp(38,3))-(6.*GABp(44,3)); 
tgabp0 = GABp(1,3); 
tgabp2 = GABp(2,3); 
tgabp5 = GABp(8,3) ;                
tgabp5p = GABp(11,3);
tgabp6 = GABp(14,3);
tgabp10n = (1/6).*(-CABp0+tgabp0-3.*tgabp2+6.*tgabp5-3*tgabp6);

%f-func
CAB1= (sqrt(3)*a./2).* (-FAB(1,3)+(2.*FAB(4,3)) +FAB(7,3)-(5.*FAB(13,3)) -(4.*FAB(19,3)) +(7.*FAB(22,3)) +(5.*FAB(28,3)) + (2.*FAB(31,3))-(4.*FAB(37,3)) + (11.*FAB(43,3)));
tfab1 = FAB(1,3);
tfab3 = FAB(4,3);
tfab4 = FAB(7,3);
tfab7n = (-1/5).*((2*CAB1/(sqrt(3).*a))+tfab1-2.*tfab3-tfab4 );

CApBp1= (sqrt(3)*a./2).* (-FApBp(1,3)+(2.*FApBp(4,3)) +FApBp(7,3)-(5.*FApBp(13,3)) -(4.*FApBp(19,3)) +(7.*FApBp(22,3)) +(5.*FApBp(28,3)) + (2.*FApBp(31,3))-(4.*FApBp(37,3)) + (11.*FApBp(43,3)));
tfapbp1 = FApBp(1,3);
tfapbp3 = FApBp(4,3);
tfapbp4 = FApBp(7,3);
tfapbp7n = (-1/5).*((2*CApBp1/(sqrt(3).*a))+tfapbp1-2.*tfapbp3-tfapbp4 );

CApB1= (sqrt(3)*a./2).* (-FApB(1,3)+(2.*FApB(4,3)) +FApB(7,3)-(5.*FApB(13,3)) -(4.*FApB(19,3)) +(7.*FApB(22,3)) +(5.*FApB(28,3)) + (2.*FApB(31,3))-(4.*FApB(37,3)) + (11.*FApB(43,3)));
tfapb1 = FApB(1,3);
tfapb3 = FApB(4,3);
tfapb4 = FApB(7,3);
tfapb7n = (-1/5).*((2*CApB1/(sqrt(3).*a))+tfapb1-2.*tfapb3-tfapb4 );

CAAp1= (sqrt(3)*a./2).* (-FAAp(1,3)+(2.*FAAp(4,3)) +FAAp(7,3)-(5.*FAAp(13,3)) -(4.*FAAp(19,3)) +(7.*FAAp(22,3)) +(5.*FAAp(28,3)) + (2.*FAAp(31,3))-(4.*FAAp(37,3)) + (11.*FAAp(43,3)));
tfaap1 = FAAp(1,3);
tfaap3 = FAAp(4,3);
tfaap4 = FAAp(7,3);
tfaap7n = (-1/5).*((2*CAAp1/(sqrt(3).*a))+tfaap1-2.*tfaap3-tfaap4 );

CBBp1= (sqrt(3)*a./2).* (-FBBp(1,3)+(2.*FBBp(4,3)) +FBBp(7,3)-(5.*FBBp(13,3)) -(4.*FBBp(19,3)) +(7.*FBBp(22,3)) +(5.*FBBp(28,3)) + (2.*FBBp(31,3))-(4.*FBBp(37,3)) + (11.*FBBp(43,3)));
tfbbp1 = FBBp(1,3);
tfbbp3 = FBBp(4,3);
tfbbp4 = FBBp(7,3);
tfbbp7n = (-1/5).*((2*CBBp1/(sqrt(3).*a))+tfbbp1-2.*tfbbp3-tfbbp4 );

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

tgabpn = zeros(12,1);   % ABp
tgabpn(1,1) = tgabp0;  
tgabpn(2,1) = tgabp2;  
tgabpn(3,1) = tgabp5 ;                
tgabpn(4,1) = tgabp5p;
tgabpn(5,1) = tgabp6;
tgabpn(6,1) = tgabp10n;
tgabpn(7,1) = tgabp10n;

tfbpbn = zeros(10,1);   % BBp
tfbpbn(1,1) = tfbbp1;
tfbpbn(2,1) = tfbbp3;
tfbpbn(3,1) = tfbbp4;
tfbpbn(4,1) = tfbbp7n;  

tfapan = zeros(10,1);   % AAp
tfapan(1,1) = tfaap1;
tfapan(2,1) = tfaap3;
tfapan(3,1) = tfaap4;
tfapan(4,1) = tfaap7n;

tfapbn = zeros(10,1);  % ApB
tfapbn(1,1) = tfapb1;
tfapbn(2,1)=  tfapb3;
tfapbn(3,1)=  tfapb4;
tfapbn(4,1)=  tfapb7n;

end

fdatan(:,1) = tfabn; 
fdatan(:,2) = tfbpapn;
dlmwrite('nn_fdata_BAp_bBN_jsv.dat',fdatan)

gdatan(:,1) = tgaan; 
gdatan(:,2) = tgbbn; 
gdatan(:,3) = tgapapn;
gdatan(:,4) = tgbpbpn;
dlmwrite('nn_gdata_BAp_bBN_jsv.dat',gdatan)

%}

hbbn =@(kx,ky) G_kxky(kx,ky,tgbbn);
habn =@(kx,ky) F_kxky(kx,ky,tfabn);
haan =@(kx,ky) G_kxky(kx,ky,tgaan);

hbpbpn =@(kx,ky) G_kxky(kx,ky,tgbpbpn);
hapbpn =@(kx,ky) F_kxky(kx,ky,tfbpapn);
hapapn =@(kx,ky) G_kxky(kx,ky,tgapapn);  

hbbpn =@(kx,ky) conj(F_kxky(kx,ky,tfbpbn));
haapn =@(kx,ky) conj(F_kxky(kx,ky,tfapan));
hbapn =@(kx,ky) F_kxky(kx,ky,tfapbn);
habpn =@(kx,ky) G_kxky(kx,ky,tgabpn);

% Hamiltonian 
 H1n =@(kx,ky) [real(haan(kx,ky)) habn(kx,ky); habn(kx,ky)' real(hbbn(kx,ky))];
 H2n =@(kx,ky) [ haapn(kx,ky) (habpn(kx,ky)); (hbapn(kx,ky)) hbbpn(kx,ky)];
 H3n =@(kx,ky) H2n(kx,ky)';
 H4n =@(kx,ky) [real(hapapn(kx,ky)) hapbpn(kx,ky);hapbpn(kx,ky)' real(hbpbpn(kx,ky))];
 
 Hn =@(kx,ky) [H1n(kx,ky) H2n(kx,ky); H3n(kx,ky) H4n(kx,ky)];

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
Ef = 0;  
Efn = 0;

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
  
legend([h1,h2,h3,h4],'DFT','W90','Full-TB','F2G2-TB')
ylabel('E (eV)');
xlabel('\Gamma-M-K-\Gamma');
title('BA-stacked BN/BN');
xlim(axes1,[0 4.029]);
box(axes1,'on');
set(axes1,'FontSize',24,'FontWeight','bold');
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.455137444489359 0.760942760942761 0.166425470332851 0.142255892255892]);
ylim([-10 10])
F1=550;F2=600;F3=550;F4=600;
figure1.Position=[F1 F2 F3 F4];


%%

% Ef = -0.1260;  
% Efn = (min(eM2Kn(:,3))+max(eM2Kn(:,2)))/2;  
% Efn = -0.1250;


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
title('\bf BA''-stacked BN/BN','Interpreter','latex');
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


%%
clc

tgaan
tgbbn
tfabn

tgapapn
tgbpbpn
tfbpapn

tfapan
tfbpbn
tfapbn
tgabpn
