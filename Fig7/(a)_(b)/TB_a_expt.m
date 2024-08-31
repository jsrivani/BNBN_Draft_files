clear all
clc

a = 2.4795; %A^o
a1 = 2.5114;
a2 = 2.4248; 


load gdata_a.dat
load fdata_a.dat

load gdata_a1.dat
load fdata_a1.dat

load gdata_a2.dat
load fdata_a2.dat

a_g = gdata_a;
a_f = fdata_a;

a1_g = gdata_a1;
a1_f = fdata_a1;

a2_g = gdata_a2;
a2_f = fdata_a2;

% Gl = ['G_0','G_1','G_2','G_2*','G_3','G_4']; 

tgaa = [a1_g(1:6,1),a_g(1:6,1),a2_g(1:6,1)];
tgbb = [a1_g(1:6,2),a_g(1:6,2),a2_g(1:6,2)];
tgapap = [a1_g(1:6,3),a_g(1:6,3),a2_g(1:6,3)];
tgbpbp = [a1_g(1:6,4),a_g(1:6,4),a2_g(1:6,4)];
tgapb = [a1_g(1:6,5),a_g(1:6,5),a2_g(1:6,5)];

tfab = [a1_f(1:4,1),a_f(1:4,1),a2_f(1:4,1)];
tfapbp = [a1_f(1:4,2),a_f(1:4,2),a2_f(1:4,2)];
tfapa = [a1_f(1:4,3),a_f(1:4,3),a2_f(1:4,3)];
tfbpb = [a1_f(1:4,4),a_f(1:4,4),a2_f(1:4,4)];
tfbpa = [a1_f(1:4,5),a_f(1:4,5),a2_f(1:4,5)];


%{
%% AA
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
plot(tgaa(:,1),'DisplayName','a = 2.51 A^o','Marker','*','MarkerSize',8,'LineWidth',3);
plot(tgaa(:,2),'DisplayName','a = 2.48 A^o','Marker','o','LineStyle','--','MarkerSize',8,'LineWidth',3);
plot(tgaa(:,3),'DisplayName','a = 2.42 A^o','Marker','+','LineStyle',':','MarkerSize',8,'LineWidth',3);

ylabel('Energy (eV)');
xlabel('n');
title('t_{AA}');
box(axes1,'on');
set(axes1,'FontSize',24,'FontWeight','bold','XTickLabel',...
    {'0','1','2','2*','3','4'});
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.60635155096012 0.629885993485342 0.225258493353028 0.166938110749186]);


%% BB
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
plot(tgbb(:,1),'DisplayName','a = 2.51 A^o','Marker','*','MarkerSize',8,'LineWidth',3);
plot(tgbb(:,2),'DisplayName','a = 2.48 A^o','Marker','o','LineStyle','--','MarkerSize',8,'LineWidth',3);
plot(tgbb(:,3),'DisplayName','a = 2.42 A^o','Marker','+','LineStyle',':','MarkerSize',8,'LineWidth',3);
ylabel('Energy (eV)');
xlabel('n');
title('t_{BB}');
box(axes1,'on');
set(axes1,'FontSize',24,'FontWeight','bold','XTickLabel',...
    {'0','1','2','2*','3','4'});
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.60635155096012 0.629885993485342 0.225258493353028 0.166938110749186]);

%% A'A'
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
plot(tgapap(:,1),'DisplayName','a = 2.51 A^o','Marker','*','MarkerSize',8,'LineWidth',3);
plot(tgapap(:,2),'DisplayName','a = 2.48 A^o','Marker','o','LineStyle','--','MarkerSize',8,'LineWidth',3);
plot(tgapap(:,3),'DisplayName','a = 2.42 A^o','Marker','+','LineStyle',':','MarkerSize',8,'LineWidth',3);
ylabel('Energy (eV)');
xlabel('n');
title('t_{A''A''}');
box(axes1,'on');
set(axes1,'FontSize',24,'FontWeight','bold','XTickLabel',...
    {'0','1','2','2*','3','4'});
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.60635155096012 0.629885993485342 0.225258493353028 0.166938110749186]);


%% B'B'
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
plot(tgbpbp(:,1),'DisplayName','a = 2.51 A^o','Marker','*','MarkerSize',8,'LineWidth',3);
plot(tgbpbp(:,2),'DisplayName','a = 2.48 A^o','Marker','o','LineStyle','--','MarkerSize',8,'LineWidth',3);
plot(tgbpbp(:,3),'DisplayName','a = 2.42 A^o','Marker','+','LineStyle',':','MarkerSize',8,'LineWidth',3);
ylabel('Energy (eV)');
xlabel('n');
title('t_{B''B''}');
box(axes1,'on');
set(axes1,'FontSize',24,'FontWeight','bold','XTickLabel',...
    {'0','1','2','2*','3','4'});
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.60635155096012 0.629885993485342 0.225258493353028 0.166938110749186]);
%}
%% A'B
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
plot(tgapb(:,1),'DisplayName','a = 2.51 A^o','Marker','*','MarkerSize',8,'LineWidth',3);
plot(tgapb(:,2),'DisplayName','a = 2.48 A^o','Marker','o','LineStyle','--','MarkerSize',8,'LineWidth',3);
plot(tgapb(:,3),'DisplayName','a = 2.42 A^o','Marker','+','LineStyle',':','MarkerSize',8,'LineWidth',3);
ylabel('Energy (eV)');
xlabel('n');
title('t_{A''B}');
box(axes1,'on');
set(axes1,'FontSize',24,'FontWeight','bold','XTickLabel',...
    {'0','1','2','2*','3','4'});
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.60635155096012 0.629885993485342 0.225258493353028 0.166938110749186]);
%}
%% AB
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
plot(tfab(:,1),'DisplayName','a = 2.51 A^o','Marker','*','MarkerSize',8,'LineWidth',3);
plot(tfab(:,2),'DisplayName','a = 2.48 A^o','Marker','o','LineStyle','--','MarkerSize',8,'LineWidth',3);
plot(tfab(:,3),'DisplayName','a = 2.42 A^o','Marker','+','LineStyle',':','MarkerSize',8,'LineWidth',3);

ylabel('Energy (eV)');
xlabel('n');
title('t_{AB}');
box(axes1,'on');
set(axes1,'FontSize',24,'FontWeight','bold','XTickLabel',...
    {'1','2','3','4'});
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.60635155096012 0.629885993485342 0.225258493353028 0.166938110749186]);


figure
plot([a2,a,a1],[tfab(1,3),tfab(1,2),tfab(1,1)],'+-')


%{
%% A'B'
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
plot(tfapbp(:,1),'DisplayName','a = 2.51 A^o','Marker','*','MarkerSize',8,'LineWidth',3);
plot(tfapbp(:,2),'DisplayName','a = 2.48 A^o','Marker','o','LineStyle','--','MarkerSize',8,'LineWidth',3);
plot(tfapbp(:,3),'DisplayName','a = 2.42 A^o','Marker','+','LineStyle',':','MarkerSize',8,'LineWidth',3);

ylabel('Energy (eV)');
xlabel('n');
title('t_{A''B''}');
box(axes1,'on');
set(axes1,'FontSize',24,'FontWeight','bold','XTickLabel',...
    {'1','2','3','4'});
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.60635155096012 0.629885993485342 0.225258493353028 0.166938110749186]);

%% A'A
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
plot(tfapa(:,1),'DisplayName','a = 2.51 A^o','Marker','*','MarkerSize',8,'LineWidth',3);
plot(tfapa(:,2),'DisplayName','a = 2.48 A^o','Marker','o','LineStyle','--','MarkerSize',8,'LineWidth',3);
plot(tfapa(:,3),'DisplayName','a = 2.42 A^o','Marker','+','LineStyle',':','MarkerSize',8,'LineWidth',3);

ylabel('Energy (eV)');
xlabel('n');
title('t_{AA''}');
box(axes1,'on');
set(axes1,'FontSize',24,'FontWeight','bold','XTickLabel',...
    {'1','2','3','4'});
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.60635155096012 0.629885993485342 0.225258493353028 0.166938110749186]);

%% B'B
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
plot(tfbpb(:,1),'DisplayName','a = 2.51 A^o','Marker','*','MarkerSize',8,'LineWidth',3);
plot(tfbpb(:,2),'DisplayName','a = 2.48 A^o','Marker','o','LineStyle','--','MarkerSize',8,'LineWidth',3);
plot(tfbpb(:,3),'DisplayName','a = 2.42 A^o','Marker','+','LineStyle',':','MarkerSize',8,'LineWidth',3);

ylabel('Energy (eV)');
xlabel('n');
title('t_{BB''}');
box(axes1,'on');
set(axes1,'FontSize',24,'FontWeight','bold','XTickLabel',...
    {'1','2','3','4'});
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.60635155096012 0.629885993485342 0.225258493353028 0.166938110749186]);
%% B'A
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
plot(tfbpa(:,1),'DisplayName','a = 2.51 A^o','Marker','*','MarkerSize',8,'LineWidth',3);
plot(tfbpa(:,2),'DisplayName','a = 2.48 A^o','Marker','o','LineStyle','--','MarkerSize',8,'LineWidth',3);
plot(tfbpa(:,3),'DisplayName','a = 2.42 A^o','Marker','+','LineStyle',':','MarkerSize',8,'LineWidth',3);

ylabel('Energy (eV)');
xlabel('n');
title('t_{AB''}');
box(axes1,'on');
set(axes1,'FontSize',24,'FontWeight','bold','XTickLabel',...
    {'1','2','3','4'});
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.60635155096012 0.629885993485342 0.225258493353028 0.166938110749186]);

%}



