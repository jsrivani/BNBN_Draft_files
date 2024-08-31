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

ndim=4; % dimension of the hamiltonian                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
N =1000; 
%% F4G4_TB

tgaan = zeros(12,1);
tgbbn = tgaan; tgapapn = tgaan; tgbpbpn = tgaan; 
tgapbn= tgaan;

tfabn = zeros(10,1);
tfbpapn = tfabn;
tfapan = tfabn;tfbpbn = tfabn;
tfbpan= tfabn;

dis = 3.1:0.1:3.5;  % interlayer distance

for i = 1:length(dis) %[3.1,3.2,3.261,3.3,3.4,3.5]
    
    dd = dis(i);
    
tgaan(1:6) = tg_AB_d4TB(11,dd);
tgbbn(1:6) = tg_AB_d4TB(22,dd);
tgapapn(1:6) = tg_AB_d4TB(33,dd);
tgbpbpn(1:6) = tg_AB_d4TB(44,dd);

tgapbn(1:6)= tg_AB_d4TB(23,dd);


tfabn(1:4) = tf_AB_d4TB(12,dd);
tfbpapn(1:4) = tf_AB_d4TB(34,dd);
tfapan(1:4) = tf_AB_d4TB(13,dd);
tfbpbn(1:4) = tf_AB_d4TB(24,dd);

tfbpan(1:4)= tf_AB_d4TB(14,dd);

% {
hbbn =@(kx,ky) G_kxky(kx,ky,tgbbn);
habn =@(kx,ky) F_kxky(kx,ky,tfabn);
haan =@(kx,ky) G_kxky(kx,ky,tgaan);

hbpbpn =@(kx,ky) G_kxky(kx,ky,tgbpbpn);
hapbpn =@(kx,ky) F_kxky(kx,ky,tfbpapn);
hapapn =@(kx,ky) G_kxky(kx,ky,tgapapn);  

hbbpn =@(kx,ky) F_kxky(kx,ky,tfbpbn);
habpn =@(kx,ky) conj(F_kxky(kx,ky,tfbpan));
hbapn =@(kx,ky) G_kxky(kx,ky,tgapbn); 
haapn =@(kx,ky) F_kxky(kx,ky,tfapan);

% Hamiltonian 
 H1n =@(kx,ky) [real(haan(kx,ky)) habn(kx,ky); habn(kx,ky)' real(hbbn(kx,ky))];
 H2n =@(kx,ky) [ haapn(kx,ky) (habpn(kx,ky)); (hbapn(kx,ky)) hbbpn(kx,ky)];
 H3n =@(kx,ky) H2n(kx,ky)';
 H4n =@(kx,ky) [real(hapapn(kx,ky)) hapbpn(kx,ky);hapbpn(kx,ky)' real(hbpbpn(kx,ky))];
 
 Hn =@(kx,ky) [H1n(kx,ky) H2n(kx,ky); H3n(kx,ky) H4n(kx,ky)];
 % }
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
Efn = 0;

x=0;

figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
% DFT

load AB_DFT_bands_31.dat
load AB_DFT_bands_32.dat
load DFT_bands3261.dat
load AB_DFT_bands_33.dat
load AB_DFT_bands_34.dat
load AB_DFT_bands_35.dat

if dd == 3.1
DFT = AB_DFT_bands_31;
Efn = -0.0670;
h1 = plot((DFT(:,1)*6.2827),DFT(:,2)+Efn,'.k','Markersize',8);
elseif dd == 3.2
DFT = AB_DFT_bands_32;  
Efn =  -0.0860;
h1 = plot((DFT(:,1)*6.2827),DFT(:,2)+Efn,'.k','Markersize',8);
elseif dd == 3.261
DFT = DFT_bands3261;  
Efn = -0.0960;
h1 = plot((DFT(:,1)*6.2827),DFT(:,2)+Efn,'.k','Markersize',8);
elseif dd == 3.3   
DFT = AB_DFT_bands_33;  
  Efn =  -0.1000;
  h1 = plot((DFT(:,1)*6.2827),DFT(:,2)+Efn,'.k','Markersize',8);
elseif dd == 3.4
DFT = AB_DFT_bands_34;
Efn = -0.1080;
h1 = plot((DFT(:,1)*6.2827),DFT(:,2)+Efn,'.k','Markersize',8);
elseif dd == 3.5    
DFT = AB_DFT_bands_35;
Efn = -0.1060;
h1 = plot((DFT(:,1)*6.2827),DFT(:,2)+Efn,'.k','Markersize',8);
end

% h1 = plot((DFT(:,1)*6.2827),DFT(:,2),'.k','Markersize',8);

line([0 ; 4],[0 ; 0],'Linestyle',':','color','r','LineWidth',1) 
  
 h4 = plot(rG2M+x,eG2Mn(:,1) ,':','Color',[0.9290 0.6940 0.1250],'LineWidth',3);
  hold on
  plot(rG2M+x,eG2Mn(:,2) ,':','Color',[0.9290 0.6940 0.1250],'LineWidth',3)
   hold on
  plot(rG2M+x,eG2Mn(:,3) ,':','Color',[0.9290 0.6940 0.1250],'LineWidth',3)
   hold on
   plot(rG2M+x,eG2Mn(:,4) ,':','Color',[0.9290 0.6940 0.1250],'LineWidth',3)
  hold on
  plot(rM2K+x,eM2Kn(:,1) ,':','Color',[0.9290 0.6940 0.1250],'LineWidth',3)
   hold on
  plot(rM2K+x,eM2Kn(:,2) ,':','Color',[0.9290 0.6940 0.1250],'LineWidth',3)
  plot(rM2K+x,eM2Kn(:,3) ,':','Color',[0.9290 0.6940 0.1250],'LineWidth',3)
  plot(rM2K+x,eM2Kn(:,4) ,':','Color',[0.9290 0.6940 0.1250],'LineWidth',3)
  plot(rK2G+x,eK2Gn(:,1) ,':','Color',[0.9290 0.6940 0.1250],'LineWidth',3)
  plot(rK2G+x,eK2Gn(:,2) ,':','Color',[0.9290 0.6940 0.1250],'LineWidth',3) 
  plot(rK2G+x,eK2Gn(:,3),':','Color',[0.9290 0.6940 0.1250],'LineWidth',3) 
  plot(rK2G+x,eK2Gn(:,4) ,':','Color',[0.9290 0.6940 0.1250],'LineWidth',3) 
  
% legend([h1,h4],'DFT','F4G4-TB')

ylabel('E (eV)');
xlabel('\Gamma-M-K-\Gamma');
title(['AB-stacked BN/BN (d = ',num2str(dd),' A^o)']);
xlim(axes1,[0 4.029]);
box(axes1,'on');
set(axes1,'FontSize',24,'FontWeight','bold');
% legend1 = legend(axes1,'show');
% set(legend1,...
%     'Position',[0.455137444489359 0.760942760942761 0.166425470332851 0.142255892255892]);
ylim([-10 10])
F1=550;F2=600;F3=550;F4=600;
figure1.Position=[F1 F2 F3 F4];

end
%%
% tgaan
% tgbbn
% tfabn 
% tgaapn
% tgbbpn
% tfapbn



















