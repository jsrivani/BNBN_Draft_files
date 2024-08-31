clear all
clc

%   AA  AB  BA  AA'  AB'  BA'

C_EK = [1.449,1.988 ,1.988 ,2.133,2.023,1.467 ];
C_EM = [1.462,1.888 ,1.888 ,1.8  ,1.531,1.86 ];

delta_E = [3.6889,4.4059,4.3940,4.3226,3.7620,3.9885];

diff_EM_EK = C_EM-C_EK;
stack = {'AA','AB','BA','AA''','AB''','BA'''};
% Only AA & BA' are having Direct bandgap 


%%
figure
plot(delta_E,'o-')
set(gca,'xtick',[1:6],'xticklabel',stack)
ylabel('E_g eV')


x = linspace(0,25);
y = delta_E;
yyaxis right
plot(diff_EM_EK,'o-');




