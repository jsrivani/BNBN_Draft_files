clear all
clc



load AA_0_d.dat
load AB_0_d.dat
load BA_0_d.dat
load AA_180_d.dat
load AB_180_d.dat
load BA_180_d.dat


figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
h1 = plot(AA_0_d(:,1),(AA_0_d(:,2)-min(AB_0_d(:,2))).*13.605662285137.*1000,'o-','linewidth',3)
hold on
h2 = plot(AB_0_d(:,1),(AB_0_d(:,2)-min(AB_0_d(:,2))).*13.605662285137.*1000,'o-','linewidth',3)
%h3 = plot(BA_0_d(:,1),(BA_0_d(:,2)-min(AB_0_d(:,2))).*13.605662285137.*1000,'-','linewidth',3)
h4 = plot(AA_180_d(:,2),(AA_180_d(:,3)-min(AB_0_d(:,2))).*13.605662285137.*1000,'o-','linewidth',3)
h5 = plot(AB_180_d(:,2),(AB_180_d(:,3)-min(AB_0_d(:,2))).*13.605662285137.*1000,'o-','linewidth',3)
h6 = plot(BA_180_d(:,2),(BA_180_d(:,3)-min(AB_0_d(:,2))).*13.605662285137.*1000,'o-','linewidth',3)

axes1 = gca;
axes1.YAxis.TickLabelInterpreter = 'latex';
axes1.YAxis.TickLabelFormat      = '\\textbf{%g}';
axes1.XAxis.TickLabelInterpreter = 'latex';
axes1.XAxis.TickLabelFormat      = '\\textbf{%g}';

%leg1 = legend([h1,h2,h3,h4,h5,h6],'\bf AA','\bf AB','\bf BA','\bf AA''','\bf AB''','\bf BA''');
leg1 = legend([h1,h2,h4,h5,h6],'\bf AA','\bf AB','\bf AA''','\bf AB''','\bf BA''');

set(leg1,'Interpreter','latex');
set(leg1,'FontSize',20);

ylabel('\bf $\Delta$ T.E (meV)','Interpreter','latex');
xlabel('\bf d (\AA)','Interpreter','latex');

title('Bilayer h-BN','Interpreter','latex');
box(axes1,'on');
set(axes1,'FontSize',24,'FontWeight','bold');
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.455137444489359 0.760942760942761 0.166425470332851 0.142255892255892]);
xlim([3 3.95])
ylim([0 106])
F1=550;F2=600;F3=550;F4=600;
figure1.Position=[F1 F2 F3 F4];
