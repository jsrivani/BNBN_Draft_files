clear all
clc

%% AB   42x42x1
% {
 a_data = [2.400,2.425,2.450,2.475,2.500,2.525,2.550]./sqrt(3);
 
 % A-B intra-layer process: 
 t1_data = [-2.9097,-2.8607,-2.7558,-2.6991,-2.6375,-2.5774,-2.5165]; % 1 n.n.
 t2_data = [-0.2346,-0.2314,-0.2259,-0.2207,-0.2154,-0.2115,-0.2059 ]; % 2 n.n.
 t3_data = [0.1049,0.0958,0.0888,0.0795,0.0711,0.0652,0.0591 ]; % 3 n.n.
 t4_data = [0.0061,0.0059,0.0046,0.0028,0.0012,0.0000,-0.0011 ]; % 4 n.n.
 t4Eff_data =[-0.0188,-0.0200,-0.0222,-0.0226,-0.0230,-0.0228,-0.0224  ];% 4 n.n.

 
t_1 = -2.6903; % at a = 2.4795 \AA
t_2 = -0.2200; 
t_3 = 0.0779; 
t_4 = 0.0025; 
t_4Eff = -0.0227;

r0ij = 2.4795./sqrt(3);
q_pi =  2.45;%3.4191;%4.1979;%3.13; %3.47; %4.1863; %2.443;%3.47; %3.37; %3.42;%2.443 ;
% t_ij =@(rij) t_0.*exp(-3.37.*(rij-r0ij)./r0ij) ;

rc = 2.5*2.4795;
lc = 0.265; %0.5 a.u;
Fc =@(rij) 1./(1+exp((rij-rc)./lc));

t_ij =@(rij,t_n,qpi) t_n.*exp(-qpi.*(rij-r0ij)./r0ij);%.*Fc(rij) ;


%{
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
plot(a_data,t1_data,'o','LineWidth',3,'DisplayName','Data : t_{AB}(n=1)','MarkerSize',10,...
    'Marker','o',...
    'LineStyle','none');
plot(a_data,t_ij(a_data,t_1,q_pi),'--r','Displayname','Fit(t_{AB}(n=1))', 'Color',[1 0 0],'LineWidth',3);
ylabel('Energy (eV)');
xlabel('r_{ij} (A^o)');
title('BN/BN : AB-stacking');
box(axes1,'on');
set(axes1,'FontSize',20,'FontWeight','bold');
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.253753988065293 0.79319895221058 0.254992319508449 0.117173524150268]);
%%
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
plot(a_data*2,t2_data,'o','LineWidth',3,'DisplayName','Data : t_{AB}(n=2)','MarkerSize',10,...
    'Marker','o',...
    'LineStyle','none');
plot(a_data*2,t_ij(a_data,t_2,q_pi),'--r','Displayname','Fit(t_{AB}(n=2))', 'Color',[1 0 0],'LineWidth',3);
ylabel('Energy (eV)');
xlabel('r_{ij} (A^o)');
title('BN/BN : AB-stacking');
box(axes1,'on');
set(axes1,'FontSize',20,'FontWeight','bold');
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.253753988065293 0.79319895221058 0.254992319508449 0.117173524150268]);
%}

% Create plots.
figure
hold on
subplot(2,1,1)
h1 = plot(a_data,t1_data,'o','LineWidth',3,'DisplayName','Data : t_{AB}(n=1)','MarkerSize',10,...
    'Marker','o',...
    'LineStyle','none');
hold on
h11 = plot(a_data,t_ij(a_data,t_1,q_pi),'--r','Displayname','Fit(t_{AB}(n=1))', 'Color',[1 0 0],'LineWidth',3);

ylabel('Energy (eV)');
title('BN/BN : AB-stacking','FontSize',20,'FontWeight','bold');
legend([h1,h11]);

subplot(2,1,2)
h2 = plot(a_data*2,t2_data,'o','LineWidth',3,'DisplayName','Data : t_{AB}(n=2)','MarkerSize',10,...
    'Marker','o',...
    'LineStyle','none');
hold on
h22 = plot(a_data*2,t_ij(a_data,t_2,q_pi),'--r','Displayname','Fit(t_{AB}(n=2))', 'Color',[1 0 0],'LineWidth',3);

ylabel('Energy (eV)');
xlabel('r_{ij} (A^o)');
box('on')
legend([h2,h22]);

%{
%%
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
plot(a_data*sqrt(7),t3_data,'o','LineWidth',3,'DisplayName','Data : t_{AB}(n=3)','MarkerSize',10,...
    'Marker','o',...
    'LineStyle','none');
plot(a_data*sqrt(7),t_ij(a_data,t_3,q_pi),'--r','Displayname','Fit(t_{AB}(n=3))', 'Color',[1 0 0],'LineWidth',3);
ylabel('Energy (eV)');
xlabel('r_{ij} (A^o)');
title('BN/BN : AB-stacking');
box(axes1,'on');
set(axes1,'FontSize',20,'FontWeight','bold');
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.253753988065293 0.79319895221058 0.254992319508449 0.117173524150268]);

%%
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
plot(a_data*sqrt(13),t4_data,'o','LineWidth',3,'DisplayName','Data : t_{AB}(n=4)','MarkerSize',10,...
    'Marker','o',...
    'LineStyle','none');
plot(a_data*sqrt(13),t_ij(a_data,t_4,q_pi),'--r','Displayname','Fit(t_{AB}(n=4))', 'Color',[1 0 0],'LineWidth',3);
ylabel('Energy (eV)');
xlabel('r_{ij} (A^o)');
title('BN/BN : AB-stacking');
box(axes1,'on');
set(axes1,'FontSize',20,'FontWeight','bold');
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.253753988065293 0.79319895221058 0.254992319508449 0.117173524150268]);

%}

 
 
 
 
 
 
