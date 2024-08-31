%                       
%                1     2
%            (A) B-----N (B)
%                '     '
%                '     '             -- AA-stack
%                '     '
%           (A') B-----N (B')
%                3     4       

%% ------------------------------------------------------------------------
clear; 
clc ;

nlayer = 2;
a = 2.4795000553; % lattice constant in Ao                                     
%% G functions
load BNNB_AA_pi1pi1.dat ;
load BNNB_AA_pi2pi2.dat ;
load BNNB_AA_pi3pi3.dat ;
load BNNB_AA_pi4pi4.dat ;
GAA = BNNB_AA_pi1pi1;
GBB = BNNB_AA_pi2pi2;
GApAp = BNNB_AA_pi3pi3;
GBpBp = BNNB_AA_pi4pi4;

load BNNB_AA_pi1pi3.dat;
load BNNB_AA_pi2pi4.dat;
GAAp = BNNB_AA_pi1pi3 ;
GBBp = BNNB_AA_pi2pi4 ;

%% F functions
load BNNB_AA_pi1pi2.dat;   
load BNNB_AA_pi3pi4.dat;
FAB = BNNB_AA_pi1pi2;  
FApBp = BNNB_AA_pi3pi4;

load BNNB_AA_pi1pi4.dat;
load BNNB_AA_pi2pi3.dat;
FABp = BNNB_AA_pi1pi4;
FApB = BNNB_AA_pi2pi3;
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

%% FULL_TB

hbb =@(kx,ky) a1G_kxky(kx,ky,tgbb);
hab =@(kx,ky) a1F_kxky(kx,ky,tfab);
haa =@(kx,ky) a1G_kxky(kx,ky,tgaa);

hbpbp =@(kx,ky) a1G_kxky(kx,ky,tgbpbp);
hapbp =@(kx,ky) a1F_kxky(kx,ky,tfbpap);
hapap =@(kx,ky) a1G_kxky(kx,ky,tgapap);  

haap =@(kx,ky) a1G_kxky(kx,ky,tgaap);
hbbp =@(kx,ky) a1G_kxky(kx,ky,tgbbp);
habp =@(kx,ky) a1F_kxky(kx,ky,tfbpa); 
hbap =@(kx,ky) conj(a1F_kxky(kx,ky,tfapb));

% Hamiltonian 
 H1 =@(kx,ky) [real(haa(kx,ky)) hab(kx,ky); hab(kx,ky)' real(hbb(kx,ky))];
 H2 =@(kx,ky) [ haap(kx,ky) (habp(kx,ky)); (hbap(kx,ky)) hbbp(kx,ky)];
 H3 =@(kx,ky) H2(kx,ky)';
 H4 =@(kx,ky) [real(hapap(kx,ky)) hapbp(kx,ky);hapbp(kx,ky)' real(hbpbp(kx,ky))];
 
 H =@(kx,ky) [H1(kx,ky) H2(kx,ky); H3(kx,ky) H4(kx,ky)];

%%
a =1;
b=4.*pi./(sqrt(3).*a);%Reciprocal lattice vector

%Define the FBZ K-space
tic

nn = 10; % Increase the K-density here
pp = 4*nn;

nnx = round(pp*2/sqrt(3));
ddnx = (2*b./sqrt(3))./nnx;

nny = pp;
ddny = b./nny;

kx1 = linspace((-b./sqrt(3))-ddnx./2,(b./sqrt(3))+ddnx./2,nnx+1);
ky1 = linspace(-b-ddny./2,b+ddny./2,nny+1);

% %%%% Conditions to confine triangular part of BZ %%%
bb=b;

cond1 = @(kx,ky)   ((sqrt(3).*kx-b)<= (ky) & (ky<=(-sqrt(3).*kx+b)));
cond2 = @(kx,ky)   ((sqrt(3).*kx+b))>=(ky) & (ky>=-sqrt(3).*kx-b);
cond3 = @(kx,ky)   ((-b./2)<=ky) & (ky<=(b./2));

% cond4 =@(kx,ky) ((ky== (sqrt(3).*kx)+b) & (ky == (-sqrt(3).*kx)-b) & (ky == -b./2));

% %%% triangular BZ
% cond4=@(kx,ky) ((-kx./sqrt(3))<= ky) & (ky < ((-kx./sqrt(3))+b));  
% cond5=@(kx,ky) (0<=kx) & (kx < (b.*sqrt(3)./2));%3./2);
% cond6 =@(kx,ky) (ky <= ((-kx./sqrt(3))+b)) & (0<=kx) &  (kx./sqrt(3)<=ky) ; % at one K point 

kk1 = 0;

temp1 = zeros((nnx.*nny).*(3/4),2);

for ii = 1:nnx
    for jj=1:nny
        kx = kx1(ii);
        ky = ky1(jj);
        if cond1(kx,ky) && cond2(kx,ky) && cond3(kx,ky)
%           if cond4(kx,ky) %&& cond5(kx,ky) %&& cond6(kx,ky)
            kk1 =kk1+1;
            temp1(kk1,1) = kx;
            temp1(kk1,2) = ky;
      
        end
    end
end

temp1((kk1+1):((nnx*nny).*(3/4)),:)=[];

kx=temp1(:,1);
ky=temp1(:,2);
scatter(kx,ky)

lkx = length(kx);
lky = length(ky);
 

band = zeros(lkx,lky,4);
hk = zeros(4,4);
emod = zeros(lkx,4);

for ii=1: lkx
   for jj = 1:lky
        kxx = kx(ii);
        kyy = ky(ii);
       
       hk = H(kxx,kyy);
      
% diagonalize Hamiltonian
   [V,D] = eig(hk) ;
        EE = diag(D);
        band(ii,jj,:) = EE(:);
       emod(ii,:) =   eig(hk) ;  
   end
end

CBM=min(emod(:,2));
VBM=max(emod(:,1));
 
 %plot bands
%%
marksize  = 15;

figure1 = figure('Color',[1 1 1]);
colormap(jet);
axes1 = axes('Parent',figure1);
hold(axes1,'on');
scatter3(kx,ky,emod(:,1),marksize.*ones(lkx,1),emod(:,1),'MarkerFaceColor','flat','MarkerEdgeColor','none');
view(axes1,[0 -90]);
box(axes1,'on');
axis(axes1,'tight');
set(axes1,'FontSize',24,'XTick',zeros(1,0),'YTick',zeros(1,0));
colorbar('peer',axes1,'southoutside');
title('VB_1')
F1=400;F2=400;F3=400;F4=400;
figure1.Position=[F1 F2 F3 F4];
savefig('AAp_bz_VB1.fig')
%%
figure2 = figure('Color',[1 1 1]);
colormap(jet);
axes2 = axes('Parent',figure2);
hold(axes2,'on');
scatter3(kx,ky,emod(:,2),marksize.*ones(lkx,1),emod(:,2),'MarkerFaceColor','flat','MarkerEdgeColor','none');
view(axes2,[0 -90]);
box(axes2,'on');
axis(axes2,'tight');
set(axes2,'FontSize',24,'XTick',zeros(1,0),'YTick',zeros(1,0));
colorbar('peer',axes2,'southoutside');
title('VB_2')
F1=400;F2=400;F3=400;F4=400;
figure2.Position=[F1 F2 F3 F4];
savefig('AAp_bz_VB2.fig')
%%
figure3 = figure('Color',[1 1 1]);
colormap(jet);
axes3 = axes('Parent',figure3);
hold(axes3,'on');
scatter3(kx,ky,emod(:,3),marksize.*ones(lkx,1),emod(:,3),'MarkerFaceColor','flat','MarkerEdgeColor','none');
view(axes3,[0 -90]);
box(axes3,'on');
axis(axes3,'tight');
set(axes3,'FontSize',24,'XTick',zeros(1,0),'YTick',zeros(1,0));
colorbar('peer',axes3,'southoutside');
title('CB_1')
F1=400;F2=400;F3=400;F4=400;
figure3.Position=[F1 F2 F3 F4];
savefig('AAp_bz_CB1.fig')
%%
figure4 = figure('Color',[1 1 1]);
colormap(jet);
axes4 = axes('Parent',figure4);
hold(axes4,'on');
scatter3(kx,ky,emod(:,4),marksize.*ones(lkx,1),emod(:,4),'MarkerFaceColor','flat','MarkerEdgeColor','none');
view(axes4,[0 -90]);
box(axes4,'on');
axis(axes4,'tight');
set(axes4,'FontSize',24,'XTick',zeros(1,0),'YTick',zeros(1,0));
colorbar('peer',axes4,'southoutside');
title('CB_2')
F1=400;F2=400;F3=400;F4=400;
figure4.Position=[F1 F2 F3 F4];
savefig('AAp_bz_CB2.fig')





%{
nlayer = 2;
a = 2.4795000553; % lattice constant in Ao                                     
%% G functions
load BNNB_AA_pi1pi1.dat ;
load BNNB_AA_pi2pi2.dat ;
load BNNB_AA_pi3pi3.dat ;
load BNNB_AA_pi4pi4.dat ;

GAA = BNNB_AA_pi1pi1;
GBB = BNNB_AA_pi2pi2;
GApAp = BNNB_AA_pi3pi3;
GBpBp = BNNB_AA_pi4pi4;

load BNNB_AA_pi1pi3.dat;
load BNNB_AA_pi2pi4.dat;

GAAp = BNNB_AA_pi1pi3 ;
GBBp = BNNB_AA_pi2pi4 ;

%% F functions
load BNNB_AA_pi1pi2.dat;   
load BNNB_AA_pi3pi4.dat;
load BNNB_AA_pi2pi3.dat;
load BNNB_AA_pi1pi4.dat;

FAB = BNNB_AA_pi1pi2;  %% ######################################################################
FApBp = BNNB_AA_pi3pi4;
FABp = BNNB_AA_pi1pi4;
FApB = BNNB_AA_pi2pi3;

%%
ndim=4; % dimension of the hamiltonian                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
N =21; 

a = 1; %2.4795000553; % lattice constant in Ao
bb = 4.*pi./((sqrt(3).*a));  % reciprocal lattice constant

kmin = -(1.5)*pi./a;
kmax = (1.5)*pi./a;

kx = linspace(kmin,kmax,N);
ky = linspace(kmin,kmax,N);

dkx = kx(2) -kx(1);
dky = ky(2) -ky(1);

band = zeros(N,N,ndim);

% %%%% Conditions to confine triangular primitive cell %%%
% cond1=@(kx,ky) ((-kx./sqrt(3))<= ky) & (ky < ((-kx./sqrt(3))+bb)); 
% cond2=@(kx,ky)  (0<=kx) & (kx < (bb.*sqrt(3)./2));%3./2); 

%%%% Conditions to confine 1st Brillouin zone %%%
% ^^^^^^^^^^^^^^^^^^^^   
cond1 = @(kx,ky)   ((sqrt(3).*kx-bb)<= (ky) & (ky<(-sqrt(3).*kx+bb)));
cond2 = @(kx,ky)   ((sqrt(3).*kx+bb))>=(ky) & (ky>-sqrt(3).*kx-bb);
cond3 = @(kx,ky)   ((-bb./2)<=ky) & (ky<(bb./2));

% ^^^^^^^^^^^^^^^^^^^^

for ii = 1:N
    for jj = 1:N
        
    kxx = kx(ii);
    kyy = ky(jj);
            
                 kv = [kxx ; kyy] ;
                 
haa = 0; hbb = 0; hab = 0; hapap = 0; hbpbp = 0; hapbp = 0;
haap = 0; hbbp = 0; habp = 0; hbap = 0;
                 
for ij = 1:300
    
     haa = haa + GAA(ij,3)*exp(1j*(dot(kv,[GAA(ij,5);GAA(ij,6)])));
     hbb = hbb + GBB(ij,3)*exp(1j*(dot(kv,[GBB(ij,5);GBB(ij,6)])));
     hab = hab + FAB(ij,3)*exp(1j*(dot(kv,[FAB(ij,5);FAB(ij,6)])));
    
     hapap = hapap + GApAp(ij,3)*exp(1j*(dot(kv,[GApAp(ij,5);GApAp(ij,6)])));
     hbpbp = hbpbp + GBpBp(ij,3)*exp(1j*(dot(kv,[GBpBp(ij,5);GBpBp(ij,6)])));
     hapbp = hapbp + FApBp(ij,3)*exp(1j*(dot(kv,[FApBp(ij,5);FApBp(ij,6)])));
     
     haap = haap + GAAp(ij,3)*exp(1j*(dot(kv,[GAAp(ij,5);GAAp(ij,6)])));
     hbbp = hbbp + GBBp(ij,3)*exp(1j*(dot(kv,[GBBp(ij,5);GBBp(ij,6)])));
    
     habp = habp + FABp(ij,3)*exp(1j*(dot(kv,[FABp(ij,5);FABp(ij,6)])));
     hbap = hbap + FApB(ij,3)*exp(1j*(dot(kv,[FApB(ij,5);FApB(ij,6)])));

end

 
H = [    real(haa)          hab         haap      habp         ;   ...
           conj(hab)    real(hbb)         hbap      hbbp         ; ...
           conj(haap)  conj(hbap)   real(hapap)     hapbp        ; ...
           conj(habp)  conj(hbbp)  conj(hapbp)  real(hbpbp)    ] ;

                         
     [V,eng]=eig(H);                             
     Vtot = V(:,:);
    E = diag(eng);
   band(ii,jj,:) = E(:);    

    end
end
%%
figure1 = figure;
colormap(jet);
axes1 = axes('Parent',figure1);
hold(axes1,'on');
% surf(kx,ky,band(:,:,1))
% surf(kx,ky,band(:,:,2))
surf(kx,ky,band(:,:,3))
surf(kx,ky,band(:,:,4))
title('AA''-stacked BN/BN :CB')
shading interp
box(axes1,'on');
grid(axes1,'on');
% view(axes1,[-37.5 30]);
view(axes1,[0 90]);
set(axes1,'FontSize',20,'FontWeight','bold');
colorbar
ylabel('k_y')
xlabel('k_x')
zlabel('E (eV)')
F1=700;F2=600;F3=700;F4=600;
figure1.Position=[F1 F2 F3 F4];
savefig(figure1,'AAp-sCB_1.fig','compact')

figure1 = figure;
colormap(jet);
axes1 = axes('Parent',figure1);
hold(axes1,'on');
surf(kx,ky,band(:,:,1))
surf(kx,ky,band(:,:,2))
% surf(kx,ky,band(:,:,3))
% surf(kx,ky,band(:,:,4))
title('AA''-stacked BN/BN:VB')
shading interp
box(axes1,'on');
grid(axes1,'on');
% view(axes1,[-37.5 30]);
view(axes1,[0 90]);
set(axes1,'FontSize',20,'FontWeight','bold');
colorbar
ylabel('k_y')
xlabel('k_x')
zlabel('E (eV)')
F1=700;F2=600;F3=700;F4=600;
figure1.Position=[F1 F2 F3 F4];
savefig(figure1,'AAp-sVB_1.fig','compact')
%% 


Ef = -max(band(:,:,2))+min(band(:,:,3))
Eg = min(Ef)  %
%%
figure2 = figure;
axes2 = axes('Parent',figure2);
hold(axes2,'on');
contour(kx,ky,band(:,:,2),27,'Linewidth',2)
box(axes2,'on');
grid(axes2,'on');
set(axes2,'FontSize',20,'FontWeight','bold');
colorbar
ylabel('k_y')
xlabel('k_x')
title('VB_1 (eV)')
F1=700;F2=600;F3=700;F4=600;
figure2.Position=[F1 F2 F3 F4];
savefig(figure2,'AAp-VB_1.fig','compact')

%%
figure3 = figure;
axes3 = axes('Parent',figure3);
hold(axes3,'on');
contour(kx,ky,band(:,:,3),33,'Linewidth',2)
box(axes3,'on');
grid(axes3,'on');
set(axes3,'FontSize',20,'FontWeight','bold');
colorbar
ylabel('k_y')
xlabel('k_x')
title('CB_1 (eV)')
F1=700;F2=600;F3=700;F4=600;
figure3.Position=[F1 F2 F3 F4];
savefig(figure3,'AAp-CB_1.fig','compact')
%}


% https://scholar.google.co.jp/citations?user=fjhc75QAAAAJ&hl=cs

% https://scholar.google.co.jp/citations?user=fjhc75QAAAAJ&hl=cs











