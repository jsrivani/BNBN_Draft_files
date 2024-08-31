%                1     2
%            (A) C-----C (B)
%                      '
%                      '             -- AB-stack
%                      '
%                 (A') C-----C (B')
%                      3     4      
%% ------------------------------------------------------------------------
clear; 
clc ;

nlayer = 2;
a = 2.4795000553; % lattice constant in Ao                                     
%% G functions
load BN_pi1pi1.dat ;
load BN_pi2pi2.dat ;

GAA = BN_pi1pi1;
GBB = BN_pi2pi2;

%% F functions
load BN_pi1pi2.dat;   

FAB = BN_pi1pi2;   

%%
ndim=2; % dimension of the hamiltonian                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
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

end
%% FULL_TB

hbb =@(kx,ky) G_kxky(kx,ky,tgbb);
hab =@(kx,ky) F_kxky(kx,ky,tfab);
haa =@(kx,ky) G_kxky(kx,ky,tgaa);


% Hamiltonian 
 H =@(kx,ky) [real(haa(kx,ky)) hab(kx,ky); hab(kx,ky)' real(hbb(kx,ky))];

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


%f-func

CAB1= (sqrt(3)*a./2).* (-FAB(1,3)+(2.*FAB(4,3)) +FAB(7,3)-(5.*FAB(13,3)) -(4.*FAB(19,3)) +(7.*FAB(22,3)) +(5.*FAB(28,3)) + (2.*FAB(31,3))-(4.*FAB(37,3)) + (11.*FAB(43,3)));
tfab1 = FAB(1,3);
tfab3n = (CAB1./(sqrt(3)*a))+(tfab1./2);

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

%f-func
CAB1= (sqrt(3)*a./2).* (-FAB(1,3)+(2.*FAB(4,3)) +FAB(7,3)-(5.*FAB(13,3)) -(4.*FAB(19,3)) +(7.*FAB(22,3)) +(5.*FAB(28,3)) + (2.*FAB(31,3))-(4.*FAB(37,3)) + (11.*FAB(43,3)));
tfab1 = FAB(1,3);
% tfab3n = (CAB1./(sqrt(3)*a))+(tfab1./2);
tfab3 = FAB(4,3);
tfab4n = (2./(sqrt(3).*a)).*CAB1+tfab1-2.*tfab3;

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

%f-func
CAB1= (sqrt(3)*a./2).* (-FAB(1,3)+(2.*FAB(4,3)) +FAB(7,3)-(5.*FAB(13,3)) -(4.*FAB(19,3)) +(7.*FAB(22,3)) +(5.*FAB(28,3)) + (2.*FAB(31,3))-(4.*FAB(37,3)) + (11.*FAB(43,3)));
tfab1 = FAB(1,3);
tfab3 = FAB(4,3);
% tfab3n = (CAB1./(sqrt(3)*a))+(tfab1./2);
% tfab4n = (2./(sqrt(3).*a)).*CAB1+tfab1-2.*tfab3;
tfab4 = FAB(7,3);
tfab7n = (-1/5).*((2*CAB1/(sqrt(3).*a))+tfab1-2.*tfab3-tfab4 );

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

end
%}
%^^^^^^^^^^^^^^^^
hbbn =@(kx,ky) G_kxky(kx,ky,tgbbn);
habn =@(kx,ky) F_kxky(kx,ky,tfabn);
haan =@(kx,ky) G_kxky(kx,ky,tgaan);

% Hamiltonian 
 Hn =@(kx,ky) [real(haan(kx,ky)) habn(kx,ky); habn(kx,ky)' real(hbbn(kx,ky))];
 
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
  hold on
  plot(rM2K+x,eM2K(:,1)-Ef,'b--','LineWidth',3)
   hold on
  plot(rM2K+x,eM2K(:,2)-Ef,'b--','LineWidth',3)
  plot(rK2G+x,eK2G(:,1)-Ef,'b--','LineWidth',3)
  plot(rK2G+x,eK2G(:,2)-Ef,'b--','LineWidth',3) 

  
 h4 = plot(rG2M+x,eG2Mn(:,1)-Efn,':','Color',[0.9290 0.6940 0.1250],'LineWidth',3);
  hold on
  plot(rG2M+x,eG2Mn(:,2)-Efn,':','Color',[0.9290 0.6940 0.1250],'LineWidth',3)
  hold on
  plot(rM2K+x,eM2Kn(:,1)-Efn,':','Color',[0.9290 0.6940 0.1250],'LineWidth',3)
   hold on
  plot(rM2K+x,eM2Kn(:,2)-Efn,':','Color',[0.9290 0.6940 0.1250],'LineWidth',3)
  plot(rK2G+x,eK2Gn(:,1)-Efn,':','Color',[0.9290 0.6940 0.1250],'LineWidth',3)
  plot(rK2G+x,eK2Gn(:,2)-Efn,':','Color',[0.9290 0.6940 0.1250],'LineWidth',3) 
 
legend([h1,h2,h3,h4],'DFT','W90','Full-TB','F4G4-TB')
% legend([h1,h2,h3],'DFT','W90','Full-TB')
% legend([h3,h4],'Full-TB','F2G2-TB')

ylabel('E (eV)');
xlabel('\Gamma-M-K-\Gamma');
title('Monolayer BN');
xlim(axes1,[0 4.029]);
box(axes1,'on');
set(axes1,'FontSize',24,'FontWeight','bold');
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.455137444489359 0.760942760942761 0.166425470332851 0.142255892255892]);
ylim([-20 12])
F1=550;F2=600;F3=550;F4=600;
figure1.Position=[F1 F2 F3 F4];

% Eg = 4.6200 eV

%%

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
  hold on
  plot(rG2M+x,eG2Mn(:,2)-Efn,'Color',n_TB,'LineWidth',2)
  hold on
  plot(rM2K+x,eM2Kn(:,1)-Efn,'Color',n_TB,'LineWidth',2)
   hold on
  plot(rM2K+x,eM2Kn(:,2)-Efn,'Color',n_TB,'LineWidth',2)
  plot(rK2G+x,eK2Gn(:,1)-Efn,'Color',n_TB,'LineWidth',2)
  plot(rK2G+x,eK2Gn(:,2)-Efn,'Color',n_TB,'LineWidth',2) 

  
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
title('\bf Monolayer BN','Interpreter','latex');
box(axes1,'on');
set(axes1,'FontSize',24,'FontWeight','bold');
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.455137444489359 0.760942760942761 0.166425470332851 0.142255892255892]);

F1=550;F2=600;F3=550;F4=600;
figure1.Position=[F1 F2 F3 F4];

xlim([0.0454  3.997])
ylim([-9.763 8.063])


xlim([1.093  2.854])
ylim([-6.194 2.028])

























