clear all
clc



num    = 3  ;   % Number of k-points.
ndim   = 4  ;     % We have four dimensions in a bilayer BN.


% --- --- --- --- --- INPUT HERE --- ---  --- --- ---

orient = 0; % Set orientation to 0 degrees (use 180 for rotated stacks)
% Stack configurations:
%   - For orient = 0: stack = 1 (AA), 2 (AB), 3 (BA)
%   - For orient = 180: stack = 1 (AA'), 2 (AB'), 3 (BA')
stack = 2; % BA stacking

% --- --- --- --- --- INPUT HERE ☝🏼 --- ---  --- --- ---

[stack,orient]

if stack == 1 && orient == 0 
dvec   =  [ 0 ; 0 ]             ;  % In units of a -AA-stack
elseif stack == 2 && orient == 0
dvec   =  [ 0 ; (1/sqrt(3)) ]   ;  % In units of a -AB-stack
elseif stack == 3 && orient == 0
dvec   =  [ 0 ; -(1/sqrt(3)) ]  ;  % In units of a -BA-stack
elseif stack == 1 && orient == 180
dvec   =  [ 0 ; 0 ]             ;  % In units of a -AA'-stack
elseif stack == 2 && orient == 180
dvec   =  [ 0 ; (1/sqrt(3)) ]   ;  % In units of a -AB'-stack
elseif stack == 3 && orient == 180
dvec   =  [ 0 ; -(1/sqrt(3)) ]  ;  % In units of a -BA'-stack
end

a = 2.4795000553;
cz     =  3.261 / a  ; 
%% K-path selection
%
n3 = 20000;
n1 = round( n3*sqrt(3) / 2 );
n2 = 10000;

n3 = 8000;
n1 = round( n3*sqrt(3) / 2 );
n2 = 4000;

n3 = 1000;
n1 = round( n3*sqrt(3) / 2 );
n2 = 500;

nsum = n1 + n2 + n3 + 1 ;

K = [  2/3    0 ];
G = [   0     0 ];
M = [  1/2   1/(2*sqrt(3)) ];
KM = K - M;

% Define kvecs  ----------------------------------------------------
kvec = zeros(3,nsum);

for jk = 1:n1
kvec(1,jk+1) = M(1) * (jk/ n1) ;
kvec(2,jk+1) = M(2) * (jk/ n1) ;
kvec(3,jk+1) = kvec(3,jk) + sqrt( (M(1) / n1)^2 + (M(2)/n1)^2 );
end


for jk = 1:n2
kvec(1,jk + n1+1) = M(1) + KM(1) * (jk/ n2) ;
kvec(2,jk + n1+1) = M(2) + KM(2) * (jk/ n2) ;
kvec(3,jk + n1+1) = kvec(3,jk + n1) + sqrt( (KM(1) / n2)^2 + (KM(2)/n2)^2 );
end


for jk = 1:n3
kvec(1,jk + n1 + n2+1) = K(1) - K(1) * (jk/ n3) ;
kvec(2,jk + n1 + n2+1) = K(2) - K(2) * (jk/ n3) ;
kvec(3,jk + n1 + n2+1) = kvec(3,jk + n1 + n2) + sqrt( (K(1) / n3)^2 + (K(2)/n3)^2 ) ;
end

kvec = kvec*2*pi; % Kpath

%%%% do not modify the following%%%%%%%%%%%%%%%%%%
num1 = num;      % dimension of array
kd = 4*pi/3 ;    % Dirac point, subtract


%% ----------------------------------------------
N    = 5 ;
L    = 2*(2*N+1)^2 ;
z    =  ones(L,1);
% Define the list of lattice sites to be used in our calculations

% Reset data
mn   = zeros(2, L)  ;
vecs = zeros(2, L)  ;
inShell = zeros(1, L) ;
shell = 3 ; % Cutoff in the disk range

a1   = [ 1 ; 0 ] ;
a2   = [ 0.5 ; sqrt(3)./2 ] ;
di   = [0; 0] ;
df   = [0; 0] ;
dfia = [0; 0] ;
dfib = [0; 0] ;

shiftv =   [ 0 ;  1./sqrt(3) ];
nco = 0  ;

for j = 1:2:L-1
    
    nco = nco + 1  ;
    m   = number2coordinate(j,N) ;
    av = a1.*m(1) + a2.*m(2)     ;
    avec(:,j)    = av(:)         ;
    avec(:,j+1)  = av(:)  + shiftv(:) ; % single layer lattice in 'a' units
    
    subA(:,j) = avec(:,j);
    subB(:,j) = avec(:,j+1);
    
    norma1 = norm(av)  ;
    norma2 = norm(av)  ;
    normvec1(nco) = norma1  ;
    normvec2(nco) = norma2  ;

% We may use this cutoff disk
    if  norma1 < shell
        inShell(j)   =  1  ;
    end
    
    if  norma2 < shell
        inShell(j+1) =  1  ;
    end 
    
    
end
indInShell = find(inShell==1);
avecc = [avec(1,indInShell);avec(2,indInShell)];
zz = z(indInShell);

%
% figure
% h1 = scatter(subA(1,:),subA(2,:),'bo','Displayname','A');
% hold on
% h2 = scatter(subB(1,:),subB(2,:),'ro','Displayname','B');
% h3 = scatter(avecc(1,:),avecc(2,:),'kx','Displayname','< cutoff');
% legend([h1,h2,h3])
% circle(0,0,shell)
% title('Monolayer graphene')
%% Intra-layer hopping energies

if stack ==1 && orient == 0
    
load fdata_AA_bBN_jsv.dat
load gdata_AA_bBN_jsv.dat   % FTB 

load nn_fdata_AA_bBN_jsv.dat
load nn_gdata_AA_bBN_jsv.dat  % F4G4

gaa(1:12,1) = nn_gdata_AA_bBN_jsv(1:12,1);
gbb(1:12,1) = nn_gdata_AA_bBN_jsv(1:12,2);
fab(1:10,1) = nn_fdata_AA_bBN_jsv(1:10,1);

fapbp = fab;
gapap = gaa; %0^o
gbpbp = gbb;
fabp(1:10,1) = nn_fdata_AA_bBN_jsv(1:10,2);
fbap = fabp;
gaap = nn_gdata_AA_bBN_jsv(1:12,3);
gbbp = nn_gdata_AA_bBN_jsv(1:12,4);

% if orient == 0
% fapbp = fab;
% gapap = gaa; %0^o
% gbpbp = gbb;
% elseif orient == 180
% fapbp = fab;
% gapap = gbb; %180^o
% gbpbp = gaa;
% end

elseif stack ==2 && orient == 0
    
load fdata_AB_bBN_jsv.dat
load gdata_AB_bBN_jsv.dat   % FTB 

load nn_fdata_AB_bBN_jsv.dat
load nn_gdata_AB_bBN_jsv.dat  % F4G4

gaa(1:12,1) = nn_gdata_AB_bBN_jsv(1:12,1);
gbb(1:12,1) = nn_gdata_AB_bBN_jsv(1:12,2);
gapap(1:12,1) = nn_gdata_AB_bBN_jsv(1:12,3);
gbpbp(1:12,1) = nn_gdata_AB_bBN_jsv(1:12,4);

fab(1:10,1) = nn_fdata_AB_bBN_jsv(1:10,1);
fapbp(1:10,1) = nn_fdata_AB_bBN_jsv(1:10,2);

elseif stack ==3 && orient == 0
    
load fdata_BA_bBN_jsv.dat
load gdata_BA_bBN_jsv.dat   % FTB 

load nn_fdata_BA_bBN_jsv.dat
load nn_gdata_BA_bBN_jsv.dat  % F4G4

gaa(1:12,1) = gdata_BA_bBN_jsv(1:12,1);
gbb(1:12,1) = gdata_BA_bBN_jsv(1:12,2);
gapap(1:12,1) = gdata_BA_bBN_jsv(1:12,3);
gbpbp(1:12,1) = gdata_BA_bBN_jsv(1:12,4);

fab(1:10,1) = fdata_BA_bBN_jsv(1:10,1);
fapbp(1:10,1) = fdata_BA_bBN_jsv(1:10,2);

elseif stack ==1 && orient == 180
    
load fdata_AAp_bBN_jsv.dat
load gdata_AAp_bBN_jsv.dat   % FTB 

load nn_fdata_AAp_bBN_jsv.dat
load nn_gdata_AAp_bBN_jsv.dat  % F4G4

gaa(1:12,1) = gdata_AAp_bBN_jsv(1:12,1);
gbb(1:12,1) = gdata_AAp_bBN_jsv(1:12,2);
gapap(1:12,1) = gdata_AAp_bBN_jsv(1:12,3);
gbpbp(1:12,1) = gdata_AAp_bBN_jsv(1:12,4);

fab(1:10,1) = fdata_AAp_bBN_jsv(1:10,1);
fapbp(1:10,1) = fdata_AAp_bBN_jsv(1:10,2);

elseif stack ==2 && orient == 180
    
load fdata_ABp_bBN_jsv.dat
load gdata_ABp_bBN_jsv.dat   % FTB 

load nn_fdata_ABp_bBN_jsv.dat
load nn_gdata_ABp_bBN_jsv.dat  % F4G4

gaa(1:12,1) = nn_gdata_ABp_bBN_jsv(1:12,1);
gbb(1:12,1) = nn_gdata_ABp_bBN_jsv(1:12,2);
gapap(1:12,1) = nn_gdata_ABp_bBN_jsv(1:12,3);
gbpbp(1:12,1) = nn_gdata_ABp_bBN_jsv(1:12,4);

fab(1:10,1) = nn_fdata_ABp_bBN_jsv(1:10,1);
fapbp(1:10,1) = nn_fdata_ABp_bBN_jsv(1:10,2);

elseif stack ==3 && orient == 180
    
load fdata_BAp_bBN_jsv.dat
load gdata_BAp_bBN_jsv.dat   % FTB 

load nn_fdata_BAp_bBN_jsv.dat
load nn_gdata_BAp_bBN_jsv.dat  % F4G4

gaa(1:12,1) = nn_gdata_BAp_bBN_jsv(1:12,1);
gbb(1:12,1) = nn_gdata_BAp_bBN_jsv(1:12,2);
gapap(1:12,1) = nn_gdata_BAp_bBN_jsv(1:12,3);
gbpbp(1:12,1) = nn_gdata_BAp_bBN_jsv(1:12,4);

fab(1:10,1) = nn_fdata_BAp_bBN_jsv(1:10,1);
fapbp(1:10,1) = nn_fdata_BAp_bBN_jsv(1:10,2);

end


%From Fengping
%{
gaa = zeros(12,1); % AA
gaa(1,1) = 1.753484;   % on-site 
gaa(2,1) = 0.007210 ;
gaa(3,1) = 0.022764 ; 
gaa(4,1) = -0.048173 ; 
gaa(5,1) = -0.002824 ;
gaa(6,1) = -0.003307 ;
gaa(7,1) = 0.000136 ;     
 
gbb = zeros(12,1); % BB
gbb(1,1) = -2.293102;  % on-site
gbb(2,1) = 0.193699 ;
gbb(3,1) = 0.019136 ; 
gbb(4,1) = -0.037515 ; 
gbb(5,1) = -0.002884 ;
gbb(6,1) = -0.004746 ;
gbb(7,1) =  0.000033;     

fab = zeros(10,1);  % AB
fab(1,1) = -2.720361 ;                       
fab(2,1) = -0.211012 ; 
fab(3,1) =  0.079806 ;
fab(4,1) =  0.005005;
fab(5,1) =  -0.008063; 
fab(6,1) =  0.015574;   

if orient == 0
fapbp = fab;
gapap = gaa; %0^o
gbpbp = gbb;
elseif orient == 180
fapbp = fab;
gapap = gbb; %180^o
gbpbp = gaa;
end
%}

%%
emod = zeros(nsum,ndim); % eigen vals 


% Uncomment this for 1d representation and comment loop
    for ii = 1:nsum   
            kx = kvec(1,ii)  ;  
            ky = kvec(2,ii)  ;
            
                 kv = [kx ; ky] ; %kvec(1:2,ii)  ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Full model  %AA
[haa]   = gfunc(kx,ky,gaa);[hbb]   = gfunc(kx,ky,gbb);[hapap]  = gfunc(kx,ky,gapap);[hbpbp] = gfunc(kx,ky,gbpbp);
[hab]   = ffunc(kx,ky,fab);[hapbp]   = ffunc(kx,ky,fapbp);
        
% ---- Initialize the interlayer H-terms ----

hbbp  = 0 ;  haap  = 0 ;  habp  = 0 ;  hbap  = 0 ;

% ----- Calculate the matrix elements using Mayou -----
nco  =   0  ;
    
for ijk = 1:2:L - 1   %  Loops in different unit cells 
    nco = nco + 1  ;
    m   = number2coordinate(ijk,N) ;
    av  = (a1.*m(1) + a2.*m(2))      ;
    avec(:,ijk)    = av(:)           ;
    avec(:,ijk+1)  = av(:)  + shiftv(:) ;

    norma1 = norm(av)  ;
    norma2 = norm(av+ shiftv)  ;
    normvec1(nco) = norma1  ;
    normvec2(nco) = norma2  ;
    
%    if  (norma1 <= shell) && (norma2 <= shell) %##
    
% --- Starting poitnt di = = (0,0). Endpoint df = d + av. dif = df - di ---    
    dif1a =   av  ;            % 0
    dif1b =   av + shiftv  ;    % +  # suppose top-layer
    dif1c =   av - shiftv  ;    % -
    
    subtA(:,nco) = dif1a;
    subtB(:,nco) = dif1b;
    subtC(:,nco) = dif1c;
    
% --- --- --- --- Interlayer hopping distances ---- ---- ---- ---- ----
    dif2a = dvec +  av ;       
    dif2b = dvec +  av + shiftv  ; % # suppose bot-layer in AB-stack with dvec = (0,aCC)
    dif2c = dvec +  av - shiftv  ;
        
    subbA(:,nco) = dif2a;
    subbB(:,nco) = dif2b;
    subbC(:,nco) = dif2c;

% {
% --- Calculate theta12a, theta21a ---
    u = [dif2a; 0];
    ux = u(1); uy = u(2);
    rxy = sqrt(ux^2+uy^2);
    theta = atan2d(ux,uy);
    
if rxy <= 8
if orient == 0    
t2_jsv = tAbA_td_wtc(dif2a,cz);
elseif orient == 180
t2_jsv = tAbB_td_wtc(dif2a,cz);
end

haap  =  haap  +  (t2_jsv) .* exp( 1j*(dot(kv,[ux,uy])) );

end
% --- Calculate theta12b, theta21b ---
    u = [dif2a; 0]; 
    ux = u(1); uy = u(2);
    rxy = sqrt(ux^2+uy^2); 
    theta = atan2d(ux,uy);
    
if rxy <= 8
if orient == 0
t2_jsv = tBbB_td_wtc(dif2a,cz);
elseif orient == 180
t2_jsv = tAbB_td_wtc(dif2a,cz);
end

hbbp  =  hbbp  +  (t2_jsv) .* exp( 1j*(dot(kv,[ux,uy])) );

end
%--------------------
    u = [dif2b; 0]; 
    ux = u(1); uy = u(2);
    rxy = sqrt(ux^2+uy^2);
    theta = atan2d(ux,uy);
    
if rxy <= 8  
if orient == 0
t2_jsv = tAbB_td_wtc(dif2b,cz);
elseif orient == 180
t2_jsv = tAbA_td_wtc(dif2b,cz);
end

habp  =  habp  +  (t2_jsv) .* exp( 1j*(dot(kv,[ux,uy])) );

end
%--------------------
    u = [dif2c; 0]; 
    ux = u(1); uy = u(2);
    rxy = sqrt(ux^2+uy^2);
    theta = atan2d(ux,uy);
    
if rxy <= 8
if orient == 0
t2_jsv = tAbB_td_wtc(dif2c,cz);
elseif orient == 180
t2_jsv = tBbB_td_wtc(dif2c,cz);
end

hbap  =  hbap  +  (t2_jsv) .* exp( 1j*(dot(kv,[ux,uy])) );

end
%}
%   end %##

end

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 H1 = [    real(haa)          hab         haap      habp         ;   ...
           conj(hab)    real(hbb)         hbap      hbbp         ; ...
           conj(haap)  conj(hbap)   real(hapap)     hapbp        ; ...
           conj(habp)  conj(hbbp)  conj(hapbp)  real(hbpbp)    ] ;
  
                         
        H = zeros(ndim) ;
        H(1:ndim,1:ndim) = H1 ; 
        
        emod(ii,:)        =   eig(H1) ;    

    end
%%
% {
figure1 = figure;
line([0,0],[-6,6])
hold on
line([-8,8],[0,0])
Oshift= 0;%0.5774;
h1 = scatter(subtA(1,:),subtA(2,:)-Oshift,'bo','Displayname','A');
hold on
h2 = scatter(subtB(1,:),subtB(2,:)-Oshift,'ro','Displayname','B');
if orient ==0
h3 = scatter(subbA(1,:),subbA(2,:)-Oshift,'b.','Linewidth',2,'Displayname','A''');
hold on
h4 = scatter(subbB(1,:),subbB(2,:)-Oshift,'r.','Linewidth',2,'Displayname','B''');
else
h3 = scatter(subbA(1,:),subbA(2,:)-Oshift,'r.','Linewidth',2,'Displayname','A''');
hold on
h4 = scatter(subbB(1,:),subbB(2,:)-Oshift,'b.','Linewidth',2,'Displayname','B''');    
end
circle(0,0,shell);
legend([h1,h2,h3,h4])
if stack ==1 && orient == 0
title('AA-stacking')   
elseif stack ==2 && orient == 0
title('AB-stacking')  
elseif stack ==3 && orient == 0
title('BA-stacking')   
elseif stack ==1 && orient == 180
title('AA''-stacking')   
elseif stack ==2 && orient == 180
title('AB''-stacking')  
elseif stack ==3 && orient == 180
title('BA''-stacking')
end
xlim([-3 3])
ylim([-3 3])
box on
annotation(figure1,'textarrow',[0.218213058419244 0.285223367697594],...
    [0.849467289719627 0.852336448598132],'String',{'|r| =3a'});

%}

%%       

const = 1/(2*pi); %*2.5346;

load AA_DFT_bands.dat
load AB_DFT_bands.dat
load BA_DFT_bands.dat
load AAp_DFT_bands.dat
load ABp_DFT_bands.dat
load BAp_DFT_bands.dat

AA = AA_DFT_bands;
AB = AB_DFT_bands;
BA = BA_DFT_bands;
AAp = AAp_DFT_bands;
ABp = ABp_DFT_bands;
BAp = BAp_DFT_bands;

KAA = reshape(AA(:,1),size(AA,1)/36,36);
EAA = reshape(AA(:,2),size(AA,1)/36,36);

KAB = reshape(AB(:,1),size(AB,1)/36,36);
EAB = reshape(AB(:,2),size(AB,1)/36,36);

KBA = reshape(BA(:,1),size(BA,1)/36,36);
EBA = reshape(BA(:,2),size(BA,1)/36,36);

KAAp = reshape(AAp(:,1),size(AAp,1)/36,36);
EAAp = reshape(AAp(:,2),size(AAp,1)/36,36);

KABp = reshape(ABp(:,1),size(ABp,1)/36,36);
EABp = reshape(ABp(:,2),size(ABp,1)/36,36);

KBAp = reshape(BAp(:,1),size(BAp,1)/36,36);
EBAp = reshape(BAp(:,2),size(BAp,1)/36,36);

% EF = 0.4030; %0.004 ;

%%

figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
hold on
if stack ==1 && orient == 0
disp('AA-stacking')
hold on
Efn = 0.030; 
% Efn = 0.3850;% FP
h2 = plot(KAA(:,:)*2.4795,EAA(:,:)-Efn,'-k','Linewidth',3);
title('AA-stacking')   
elseif stack ==2 && orient == 0
disp('AB-stacking')
hold on
Efn = 0.0720;%+0.0960;
% Efn = 0.4400;% FP
h2 = plot(KAB(:,:)*2.4795,EAB(:,:)-Efn,'-k','Linewidth',3);
title('AB-stacking')  
elseif stack ==3 && orient == 0
disp('BA-stacking')
hold on
Efn = 0.0720;%+0.0960;
% Efn = 0.4400;% FP
h2 = plot(KBA(:,:)*2.4795,EBA(:,:)-Efn,'-k','Linewidth',3);
title('BA-stacking')   
elseif stack ==1 && orient == 180
disp('AAp-stacking')
hold on
Efn = 0.0130;%+0.1140;
% Efn = 0.4160;% FP
h2 = plot(KAAp(:,:)*2.4795,EAAp(:,:)-Efn,'-k','Linewidth',3);
title('AA''-stacking')   
elseif stack ==2 && orient == 180
disp('ABp-stacking')
hold on
Efn = 0;
% Efn = 0.3840;% FP
h2 = plot(KABp(:,:)*2.4795,EABp(:,:)-Efn ,'-k','Linewidth',3);
title('AB''-stacking')  
elseif stack ==3 && orient == 180
disp('BAp-stacking')
hold on
Efn = 0.0150;
% Efn = 0.3830;% FP
h2 = plot(KBAp(:,:)*2.4795,EBAp(:,:)-Efn ,'-k','Linewidth',3);
title('BA''-stacking')
end
plot(kvec(3,:)*const,emod(:,ndim/2-1) ,"r-",'Linewidth',2);
   plot(kvec(3,:)*const,emod(:,ndim/2) ,"r-",'Linewidth',2)
   plot(kvec(3,:)*const,emod(:,ndim/2+1) ,"r-",'Linewidth',2) 
   plot(kvec(3,:)*const,emod(:,ndim/2+2) ,"r-",'Linewidth',2)

   h1 = plot(kvec(3,1)*const,emod(1,ndim/2-1) ,"r-",'Linewidth',2);
   h2 = plot(KBAp(1,1)*2.4795,EBAp(1,1) ,'-k','Linewidth',3);

%    legend([h1,h2],'t(r,\theta) model)','DFT')
   legend([h1,h2],'TC : t(d) model)','DFT')

xlabel('\Gamma-M-K-\Gamma');
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',20,'FontWeight','bold','XTick',zeros(1,0));
xlim([0 1.577])
ylim([-9.5 9.5])
F1=550;F2=500;F3=550;F4=500;
figure1.Position=[F1 F2 F3 F4];
hold off
%%
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
hold on
if stack ==1 && orient == 0
disp('AA-stacking')
hold on
Efn = 0.030; 
% Efn = 0.3850;% FP
h2 = plot(KAA(:,:)*2.4795,EAA(:,:)-Efn,'-k','Linewidth',3);
title('AA-stacking')   
elseif stack ==2 && orient == 0
disp('AB-stacking')
hold on
Efn = 0.0720+0.0960;
% Efn = 0.4400;% FP
h2 = plot(KAB(:,:)*2.4795,EAB(:,:)-Efn,'-k','Linewidth',3);
title('AB-stacking')  
elseif stack ==3 && orient == 0
disp('BA-stacking')
hold on
Efn = 0.0720;%+0.0960;
% Efn = 0.4400;% FP
h2 = plot(KBA(:,:)*2.4795,EBA(:,:)-Efn,'-k','Linewidth',3);
title('BA-stacking')   
elseif stack ==1 && orient == 180
disp('AAp-stacking')
hold on
Efn = 0.0130+0.0360+0.1140;
% Efn = 0.4160;% FP
h2 = plot(KAAp(:,:)*2.4795,EAAp(:,:)-Efn,'-k','Linewidth',3);
title('AA''-stacking')   
elseif stack ==2 && orient == 180
disp('ABp-stacking')
hold on
Efn = 0;
% Efn = 0.3840;% FP
h2 = plot(KABp(:,:)*2.4795,EABp(:,:)-Efn ,'-k','Linewidth',3);
title('AB''-stacking')  
elseif stack ==3 && orient == 180
disp('BAp-stacking')
hold on
Efn = 0.0150;
% Efn = 0.3830;% FP
h2 = plot(KBAp(:,:)*2.4795,EBAp(:,:)-Efn ,'-k','Linewidth',3);
title('BA''-stacking')
end
plot(kvec(3,:)*const,emod(:,ndim/2-1) ,"r-",'Linewidth',2);
   plot(kvec(3,:)*const,emod(:,ndim/2) ,"r-",'Linewidth',2)
   plot(kvec(3,:)*const,emod(:,ndim/2+1) ,"r-",'Linewidth',2) 
   plot(kvec(3,:)*const,emod(:,ndim/2+2) ,"r-",'Linewidth',2)

   h1 = plot(kvec(3,1)*const,emod(1,ndim/2-1) ,"r-",'Linewidth',2);
   h2 = plot(KBAp(1,1)*2.4795,EBAp(1,1) ,'-k','Linewidth',3);

%    legend([h1,h2],'t(r,\theta) model)','DFT')
   legend([h1,h2],'TC : t(d) model)','DFT')

xlabel('\Gamma-M-K-\Gamma');
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',20,'FontWeight','bold','XTick',zeros(1,0));
xlim([0.42 1.15])
ylim([-4.3 3.5])
F1=550;F2=500;F3=550;F4=500;
figure1.Position=[F1 F2 F3 F4];
hold off


% end
% end

