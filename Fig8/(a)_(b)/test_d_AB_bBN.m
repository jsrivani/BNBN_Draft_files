% ~~~~~~ F4G4-parameters as a function of 'd' ~~~~~~

clear all
clc

% optimized 'd' = 3.261 A^o

AA_opt_d = [1.6729,0.0234,0.0273,0.0273,-0.0419,-0.0014,-0.0014 ];
BB_opt_d = [-2.3030,0.1892,0.0210,0.0210,-0.0367,0.0008,0.0008 ];
ApAp_opt_d = [1.7257,0.0106,0.0257,0.0257,-0.0440,-0.0013,-0.0013 ];
BpBp_opt_d = [-2.2097,0.1891,0.0146,0.0146,-0.0373,0.0007,0.0007 ];
ApB_opt_d = [0.3775,-0.0619,-0.0245,-0.0160,-0.0116,-0.0030,-0.0030 ];

AB_opt_d = [-2.6903,-0.2200,0.0779,-0.0227 ];
ApBp_opt_d = [-2.7087,-0.2034,0.0794,-0.0206 ];

BBp_opt_d = [-0.0177,-0.0467,0.0206,0.0072 ];
AAp_opt_d = [0.4835,0.0744,0.0520, 0.0218  ];
ABp_opt_d = [0.1190,-0.1383,-0.0215,-0.0023  ];

% d = 3.1 A^o

AA_d31 = [1.6528,0.0261,0.0240,0.0240,-0.0413,-0.0016,-0.0016   ];
BB_d31 = [-2.3318,0.1925,0.0232, 0.0232,-0.0360,0.0006,0.0006   ];
ApAp_d31 = [1.8491,0.0143,0.0248,0.0248,-0.0438, -0.0016,-0.0016  ];
BpBp_d31 = [-2.1440,0.1967,0.0172,0.0172,-0.0360,0.0010,0.0010  ];
ApB_d31 = [0.5006,-0.0527,-0.0174,-0.0067, -0.0057,-0.0037,-0.0037   ];

AB_d31 = [-2.6864,-0.2259,0.0765,-0.0231  ];
ApBp_d31 = [-2.6968,-0.2092,0.0768,-0.0227   ];

BBp_d31 = [-0.0100,-0.0658,0.0159,0.0061   ];
AAp_d31 = [0.4709,0.0362,0.0388, 0.0199    ];
ABp_d31 = [0.2029,-0.1360,-0.0063,0.0033   ];

% d = 3.2 A^o 
AA_d32 = [1.6567,0.0228,0.0258,0.0258,-0.0412,-0.0016,-0.0016   ];
BB_d32 = [-2.3130,0.1901,0.0222,0.0222,-0.0369,0.0006,0.0006  ];
ApAp_d32 = [1.7925,0.0129, 0.0260,0.0260,-0.0434,-0.0015,-0.0015   ];
BpBp_d32 = [-2.1733,0.1938,0.0168,0.0168,-0.0369,0.0009,0.0009  ];
ApB_d32 = [0.4148,-0.0594,-0.0207,-0.013,-0.0089,-0.0039,-0.003   ];

AB_d32 = [-2.6868,-0.2241,0.0783,-0.0227  ];
ApBp_d32 = [-2.7001,-0.2092,0.0784,-0.0220   ];

BBp_d32 = [-0.0157,-0.0546, 0.0186,0.0068   ];
AAp_d32 = [0.4768,0.0590,0.0463,0.0211  ];
ABp_d32 = [ 0.1558,-0.1363,-0.0146,0.0003   ];

% d = 3.3 A^o

AA_d33 = [1.6659,0.0175,0.0259,0.0259,-0.0422,-0.0018,-0.0018  ];
BB_d33 = [-2.3184,0.1883,0.0208,0.0208,-0.0379,0.0005,0.0005  ];
ApAp_d33 = [1.7539,0.0090,0.0256,0.0256,-0.0441,-0.0017,-0.0017 ];
BpBp_d33 = [-2.2196, 0.1906,0.0160,0.0160,-0.0380,0.0007,0.0007  ];
ApB_d33 = [0.3412,-0.0627,-0.0226,-0.0176,-0.0110,-0.0039,-0.0039  ];

AB_d33 = [-2.6894,-0.2215,0.0805,-0.0223  ];
ApBp_d33 = [-2.7035,-0.2083,0.0807, -0.0215  ];

BBp_d33 = [-0.0165,-0.0443,0.0205,0.0075  ];
AAp_d33 = [0.4628,0.0720,0.0498,0.0218   ];
ABp_d33 = [0.1174,-0.1323,-0.0196,-0.0022  ];

% d = 3.4 A^o

AA_d34 = [1.6789,0.0119,0.0255,0.0255,-0.0434, -0.0020,-0.0020   ];
BB_d34 = [-2.2902,0.1847,0.0191,0.0191, -0.0390,0.0004,0.0004   ];
ApAp_d34 = [1.7303,0.0044,0.0248,0.0248,-0.0451,-0.0017,-0.0017   ];
BpBp_d34 = [-2.2165,0.1857, 0.0147,0.0147,-0.0393,0.0005,0.0005   ];
ApB_d34 = [0.2729,-0.0661,-0.0248,-0.0210,-0.0129,-0.0031,-0.0031  ];

AB_d34 = [-2.6840,-0.2174,0.0827,-0.0219  ];
ApBp_d34 = [-2.6968,-0.2057,0.0831,-0.0208  ];

BBp_d34 = [ -0.0137,-0.0341, 0.0223,0.0084 ];
AAp_d34 = [ 0.4466,0.0829,0.0530,0.0222   ];
ABp_d34 = [0.0815,-0.1277, -0.0242,-0.0046  ];

% d = 3.5 A^o

AA_d35 = [1.6864,0.0066,0.0245,0.0245,-0.0448,-0.0021,-0.0021   ];
BB_d35 = [-2.2867,0.1822,0.0175,0.0175,-0.0400,0.0003,0.0003   ];
ApAp_d35 = [1.7018,-0.0006,0.0234,0.0234,-0.0462,-0.0019,-0.0019   ];
BpBp_d35 = [-2.2415,0.1818,0.0134,0.0134,-0.0405,0.0002,0.0002   ];
ApB_d35 = [0.215,-0.0688, -0.0265,-0.024,-0.0145, -0.0032,-0.0032   ];

AB_d35 = [-2.6833,-0.2138,0.0847,-0.0210  ];
ApBp_d35 = [-2.6953,-0.2033,0.0856,-0.0195   ];

BBp_d35 = [-0.0067,-0.0242,0.0241,0.0091   ];
AAp_d35 = [0.4276,0.0906,0.0551, 0.0223  ];
ABp_d35 = [ 0.0509,-0.1233,-0.0281,-0.0069   ];


%%
% figure
% plot(AA_opt_d,'Displayname','opt-d','linewidth',1)
% hold on
% plot(AA_d31,'Displayname','d=3.1','linewidth',1)
% plot(AA_d32,'Displayname','d=3.2','linewidth',1)
% plot(AA_d33,'Displayname','d=3.3','linewidth',1)
% plot(AA_d34,'Displayname','d=3.4','linewidth',1)
% plot(AA_d35,'Displayname','d=3.5','linewidth',1)
% legend('show')

%%
d = [3.1,3.2,3.3,3.4,3.5]';
AA_d = [AA_d31',AA_d32',AA_d33',AA_d34',AA_d35' ]; 
BB_d = [BB_d31',BB_d32',BB_d33',BB_d34',BB_d35' ]; 
ApAp_d = [ApAp_d31',ApAp_d32',ApAp_d33',ApAp_d34',ApAp_d35' ]; 
BpBp_d = [BpBp_d31',BpBp_d32',BpBp_d33',BpBp_d34',BpBp_d35' ]; 

AB_d = [AB_d31',AB_d32',AB_d33',AB_d34',AB_d35' ]; 
ApBp_d = [ApBp_d31',ApBp_d32',ApBp_d33',ApBp_d34',ApBp_d35' ];

AAp_d = [AAp_d31',AAp_d32',AAp_d33',AAp_d34',AAp_d35' ]; 
BBp_d = [BBp_d31',BBp_d32',BBp_d33',BBp_d34',BBp_d35' ]; 

ApB_d = [ApB_d31',ApB_d32',ApB_d33',ApB_d34',ApB_d35' ]; 
ABp_d = [ABp_d31',ABp_d32',ABp_d33',ABp_d34',ABp_d35' ]; 

%% t_AA
%{
G0 = AA_d(1,:)';
G1 = AA_d(2,:)';
G2 = AA_d(3,:)';G2(1) = 0.025;
G2p = AA_d(4,:)';G2p(1) = 0.025;
G3 = AA_d(5,:)';
G4 = AA_d(6,:)';

     func=@(x,a0,b0,c0,d0) a0.*exp(b0.*x) + c0.*exp(d0*x);

%        a0 =        3395;
%        b0 =      -3.551;
%        c0 =        1.07;
%        d0 =       0.129;
       a0 =       1.402;
       b0 =     0.05266;
       c0 =    -0.09881;
       d0 =      -9.767;

       a1 =       1.418;
       b1 =    -0.03196;
       c1 =      -1.231;
       d1 =    0.006911;
%        a1 =     -0.5444;
%        b1 =      0.2157;
%        c1 =      0.6494;
%        d1 =       0.167;
       
       a2 =      0.5976;
       b2 =     -0.8391;
       c2 =      -54.38;
       d2 =      -2.562;
       
       a3 =           0;
       b3 =      -18.65;
       c3 =    -0.02072 ;
       d3 =      0.2181 ;

       a4 =      -4.332 ;
       b4 =      -7.429 ;
       c4 =  -0.0005356 ;
       d4 =      0.3698 ;


figure
hold on

plot(d,G0,'Displayname','G0','Linewidth',3)
plot(d,func(d,a0,b0,c0,d0),'--ro','Displayname','G0-fit','Markersize',5,'Linewidth',1.5)

% plot(d,G1,'Displayname','G1','Linewidth',3)
% plot(d,func(d,a1,b1,c1,d1),'--ro','Displayname','G1-fit','Markersize',5,'Linewidth',1.5)

% plot(d,G2,'Displayname','G2','Linewidth',3)
% plot(d,func(d,a2,b2,c2,d2),'--ro','Displayname','G2-fit','Markersize',5,'Linewidth',1.5)

% plot(d,G2p,'Displayname','G2p','Linewidth',3)
% plot(d,func(d,a2p,b2p,c2p,d2p),'--ro','Displayname','G2p-fit','Markersize',5,'Linewidth',1.5)

% plot(d,G3,'Displayname','G3','Linewidth',3)
% plot(d,func(d,a3,b3,c3,d3),'--ro','Displayname','G3-fit','Markersize',5,'Linewidth',1.5)
 
% plot(d,G4,'Displayname','G4','Linewidth',3)
% plot(d,func(d,a4,b4,c4,d4),'--ro','Displayname','G4-fit','Markersize',5,'Linewidth',1.5)

legend('show')
       
x = d;
x0 = [a0,b0,c0,d0]; 
x1 = [a1,b1,c1,d1]; 
x2 = [a2,b2,c2,d2]; 
% x2p = [a2p,b2p,c2p,d2p];
x3 = [a3,b3,c3,d3]; 
x4 = [a4,b4,c4,d4]; 

myfit1_BN(x0,G0)
%  myfit1_BN(x1,G1)
% myfit1_BN(x2,G2)
% myfit1_BN(x2p,G2p)
% myfit1_BN(x3,G3)
% myfit1_BN(x4,G4)

%%
figure1 = figure;
axes1 = axes('Parent',figure1,...
    'Position',[0.13 0.680241327300151 0.77 0.244758672699849]);
hold(axes1,'on');
plot(d,AA_d(1,:),'Parent',axes1,'DisplayName','G0','Marker','o','LineWidth',2);
% ylabel('E (eV)');

title('t_{AA} (eV)');

box(axes1,'on');
set(axes1,'FontSize',20,'FontWeight','bold','XTick',zeros(1,0));
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.771028961154493 0.750303835144163 0.0970149253731343 0.0392156862745098]);
axes2 = axes('Parent',figure1,...
    'Position',[0.13 0.11 0.775 0.515367647058824]);
hold(axes2,'on');
plot(d,AA_d(2,:),'o-','Displayname','G1','linewidth',2)
plot(d,AA_d(3,:),'o-','Displayname','G2','linewidth',2)
plot(d,AA_d(4,:),'o-','Displayname','G2*','linewidth',2)
plot(d,AA_d(5,:),'o-','Displayname','G3','linewidth',2)
plot(d,AA_d(6,:),'o-','Displayname','G4','linewidth',2)
% ylabel('E (eV)','FontWeight','bold');
xlabel('d A^o','FontWeight','bold');
box(axes2,'on');
set(axes2,'FontSize',20,'FontWeight','bold');
legend2 = legend(axes2,'show');
set(legend2,...
    'Position',[0.772164426213304 0.20594807893233 0.107462686567164 0.174962292609351]);
F1=550;F2=600;F3=550;F4=600;
figure1.Position=[F1 F2 F3 F4];
%}

%% t_BB
%{
G0 = BB_d(1,:)';
G1 = BB_d(2,:)';
G2 = BB_d(3,:)';
G2p = BB_d(4,:)';
G3 = BB_d(5,:)';
G4 = BB_d(6,:)';%G4(5) = 1.0e-03 *0.7;

     func=@(x,a0,b0,c0,d0) a0.*exp(b0.*x) + c0.*exp(d0*x);

       a0 =     0.01081 ;
       b0 =      0.6491 ;
       c0 =      -2.582 ;
       d0 =    -0.02213 ;

       a1 =    -0.05535 ;
       b1 =      0.1875 ;
       c1 =      0.3135 ;
       d1 =    -0.02322 ;

       a2 =      0.1805 ;
       b2 =     -0.6548 ;
       c2 =    -0.03781 ;
       d2 =      -1.454 ;

       a2p =      0.1805 ;
       b2p =     -0.6548 ;
       c2p =    -0.03781 ;
       d2p =      -1.454 ;
       
       a3 =     -0.0798 ;
       b3 =     -0.8135 ;
       c3 =   -0.007377 ;
       d3 =      0.4485 ;

       a4 =      0.2453;
       b4 =      -1.451;
       c4 =    -0.09688;
       d4 =      -1.252;

       
figure
hold on

% plot(d,G0,'Displayname','G0','Linewidth',3)
% plot(d,func(d,a0,b0,c0,d0),'--ro','Displayname','G0-fit','Markersize',5,'Linewidth',1.5)

% plot(d,G1,'Displayname','G1','Linewidth',3)
% plot(d,func(d,a1,b1,c1,d1),'--ro','Displayname','G1-fit','Markersize',5,'Linewidth',1.5)

% plot(d,G2,'Displayname','G2','Linewidth',3)
% plot(d,func(d,a2,b2,c2,d2),'--ro','Displayname','G2-fit','Markersize',5,'Linewidth',1.5)

% plot(d,G2p,'Displayname','G2p','Linewidth',3)
% plot(d,func(d,a2p,b2p,c2p,d2p),'--ro','Displayname','G2p-fit','Markersize',5,'Linewidth',1.5)

% plot(d,G3,'Displayname','G3','Linewidth',3)
% plot(d,func(d,a3,b3,c3,d3),'--ro','Displayname','G3-fit','Markersize',5,'Linewidth',1.5)
 
plot(d,G4,'Displayname','G4','Linewidth',3)
plot(d,func(d,a4,b4,c4,d4),'--ro','Displayname','G4-fit','Markersize',5,'Linewidth',1.5)

legend('show')
       
x = d;
x0 = [a0,b0,c0,d0]; 
x1 = [a1,b1,c1,d1]; 
x2 = [a2,b2,c2,d2]; 
x2p = [a2p,b2p,c2p,d2p];
x3 = [a3,b3,c3,d3]; 
x4 = [a4,b4,c4,d4]; 

% myfit1_BN(x0,G0)
% myfit1_BN(x1,G1)
% myfit1_BN(x2,G2)
% myfit1_BN(x2p,G2p)
% myfit1_BN(x3,G3)
myfit1_BN(x4,G4)


% figure1 = figure;
% axes1 = axes('Parent',figure1,...
%     'Position',[0.13 0.680241327300151 0.77 0.244758672699849]);
% hold(axes1,'on');
% plot(d',G0,'Parent',axes1,'DisplayName','G0','Marker','o','LineWidth',2);
% % ylabel('E (eV)');
% 
% title('t_{BB} (eV)');
% 
% box(axes1,'on');
% set(axes1,'FontSize',20,'FontWeight','bold','XTick',zeros(1,0));
% legend1 = legend(axes1,'show');
% set(legend1,...
%     'Position',[0.771028961154493 0.750303835144163 0.0970149253731343 0.0392156862745098]);
% axes2 = axes('Parent',figure1,...
%     'Position',[0.13 0.11 0.775 0.515367647058824]);
% hold(axes2,'on');
% plot(d',G1,'o-','Displayname','G1','linewidth',2)
% plot(d',G2,'o-','Displayname','G2','linewidth',2)
% plot(d',G2p,'o-','Displayname','G2*','linewidth',2)
% plot(d',G3,'o-','Displayname','G3','linewidth',2)
% plot(d',G4,'o-','Displayname','G4','linewidth',2)
% % ylabel('E (eV)','FontWeight','bold');
% xlabel('d A^o','FontWeight','bold');
% box(axes2,'on');
% set(axes2,'FontSize',20,'FontWeight','bold');
% legend2 = legend(axes2,'show');
% set(legend2,...
%     'Position',[0.772164426213304 0.20594807893233 0.107462686567164 0.174962292609351]);
% F1=550;F2=600;F3=550;F4=600;
% figure1.Position=[F1 F2 F3 F4];
%}

%% t_ApAp
%{
G0 = ApAp_d(1,:)';
G1 = ApAp_d(2,:)';
G2 = ApAp_d(3,:)';
G2p = ApAp_d(4,:)';
G3 = ApAp_d(5,:)';
G4 = ApAp_d(6,:)';

     func=@(x,a0,b0,c0,d0) a0.*exp(b0.*x) + c0.*exp(d0*x);

%        a0 =       239.8;
%        b0 =      -1.974;
%        c0 =      0.5514;
%        d0 =      0.2821;       
       a0 =       2.514;
       b0 =      -6.322;
       c0 =       3.671;
       d0 =     -0.2225;

       a1 =     -0.1604;
       b1 =      0.0707;
       c1 =      0.2551;
       d1 =    -0.05683;

       a2 =     -0.1176  ;
       b2 =      0.3172  ;
       c2 =      0.1424  ;
       d2 =       0.281  ;

       a2p =     -0.1176  ;
       b2p =      0.3172  ;
       c2p =      0.1424  ;
       d2p =       0.281  ; 
       
       a3 =           0;
       b3 =      -13.68;
       c3 =    -0.02733;
       d3 =      0.1478;
      
       a4 =      -4.028 ;
       b4 =       -3.73 ;
       c4 =  -0.0003268 ;
       d4 =      0.4917 ;

       
figure
hold on

plot(d,G0,'Displayname','G0','Linewidth',3)
plot(d,func(d,a0,b0,c0,d0),'--ro','Displayname','G0-fit','Markersize',5,'Linewidth',1.5)

% plot(d,G1,'Displayname','G1','Linewidth',3)
% plot(d,func(d,a1,b1,c1,d1),'--ro','Displayname','G1-fit','Markersize',5,'Linewidth',1.5)

% plot(d,G2,'Displayname','G2','Linewidth',3)
% plot(d,func(d,a2,b2,c2,d2),'--ro','Displayname','G2-fit','Markersize',5,'Linewidth',1.5)

% plot(d,G2p,'Displayname','G2p','Linewidth',3)
% plot(d,func(d,a2p,b2p,c2p,d2p),'--ro','Displayname','G2p-fit','Markersize',5,'Linewidth',1.5)
% 
% plot(d,G3,'Displayname','G3','Linewidth',3)
% plot(d,func(d,a3,b3,c3,d3),'--ro','Displayname','G3-fit','Markersize',5,'Linewidth',1.5)
%  
% plot(d,G4,'Displayname','G4','Linewidth',3)
% plot(d,func(d,a4,b4,c4,d4),'--ro','Displayname','G4-fit','Markersize',5,'Linewidth',1.5)

legend('show')
       
x = d;
x0 = [a0,b0,c0,d0]; 
x1 = [a1,b1,c1,d1]; 
x2 = [a2,b2,c2,d2]; 
x2p = [a2p,b2p,c2p,d2p];
x3 = [a3,b3,c3,d3]; 
x4 = [a4,b4,c4,d4]; 

myfit1_BN(x0,G0)
% myfit1_BN(x1,G1)
% myfit1_BN(x2,G2)
% myfit1_BN(x2p,G2p)
% myfit1_BN(x3,G3)
% myfit1_BN(x4,G4)      
      
% figure1 = figure;
% axes1 = axes('Parent',figure1,...
%     'Position',[0.13 0.680241327300151 0.77 0.244758672699849]);
% hold(axes1,'on');
% plot(d,ApAp_d(1,:),'Parent',axes1,'DisplayName','G0','Marker','o','LineWidth',2);
% % ylabel('E (eV)');
% 
% title('t_{ApAp} (eV)');
% 
% box(axes1,'on');
% set(axes1,'FontSize',20,'FontWeight','bold','XTick',zeros(1,0));
% legend1 = legend(axes1,'show');
% set(legend1,...
%     'Position',[0.771028961154493 0.750303835144163 0.0970149253731343 0.0392156862745098]);
% axes2 = axes('Parent',figure1,...
%     'Position',[0.13 0.11 0.775 0.515367647058824]);
% hold(axes2,'on');
% plot(d,ApAp_d(2,:),'o-','Displayname','G1','linewidth',2)
% plot(d,ApAp_d(3,:),'o-','Displayname','G2','linewidth',2)
% plot(d,ApAp_d(4,:),'o-','Displayname','G2*','linewidth',2)
% plot(d,ApAp_d(5,:),'o-','Displayname','G3','linewidth',2)
% plot(d,ApAp_d(6,:),'o-','Displayname','G4','linewidth',2)
% % ylabel('E (eV)','FontWeight','bold');
% xlabel('d A^o','FontWeight','bold');
% box(axes2,'on');
% set(axes2,'FontSize',20,'FontWeight','bold');
% legend2 = legend(axes2,'show');
% set(legend2,...
%     'Position',[0.772164426213304 0.20594807893233 0.107462686567164 0.174962292609351]);
% F1=550;F2=600;F3=550;F4=600;
% figure1.Position=[F1 F2 F3 F4];
%}

%% t_BpBp
%{
G0 = BpBp_d(1,:)';
G1 = BpBp_d(2,:)';
G2 = BpBp_d(3,:)';
G2p = BpBp_d(4,:)';
G3 = BpBp_d(5,:)';
G4 = BpBp_d(6,:)';

func=@(x,a0,b0,c0,d0) a0.*exp(b0.*x) + c0.*exp(d0*x);

       a0 =     0.07838;
       b0 =      0.9919;
       c0 =      -0.726;
       d0 =      0.5373;
       
       a1 =      0.3536;
       b1 =     -0.1884;
       c1 =      0.8459;
       d1 =      -4.853;

       a2 =  -0.0003819;
       b2 =       1.497;
       c2 =    0.002532;
       d2 =       1.003;
       
       a3 =     -0.2003 ;
       b3 =     -0.9921 ;
       c3 =   -0.0036591 ;
       d3 =      0.6417 ;

       a4 =     0.04211;
       b4 =      0.1121;
       c4 =    -0.03772;
       d4 =       0.142;

       
figure
hold on

% plot(d,G0,'Displayname','G0','Linewidth',3)
% plot(d,func(d,a0,b0,c0,d0),'--ro','Displayname','G0-fit','Markersize',5,'Linewidth',1.5)

% plot(d,G1,'Displayname','G1','Linewidth',3)
% plot(d,func(d,a1,b1,c1,d1),'--ro','Displayname','G1-fit','Markersize',5,'Linewidth',1.5)

% plot(d,G2,'Displayname','G2','Linewidth',3)
% plot(d,func(d,a2,b2,c2,d2),'--ro','Displayname','G2-fit','Markersize',5,'Linewidth',1.5)

% plot(d,G2p,'Displayname','G2p','Linewidth',3)
% plot(d,func(d,a2p,b2p,c2p,d2p),'--ro','Displayname','G2p-fit','Markersize',5,'Linewidth',1.5)

% plot(d,G3,'Displayname','G3','Linewidth',3)
% plot(d,func(d,a3,b3,c3,d3),'--ro','Displayname','G3-fit','Markersize',5,'Linewidth',1.5)
 
plot(d,G4,'Displayname','G4','Linewidth',3)
plot(d,func(d,a4,b4,c4,d4),'--ro','Displayname','G4-fit','Markersize',5,'Linewidth',1.5)

legend('show')
       
x = d;
x0 = [a0,b0,c0,d0]; 
x1 = [a1,b1,c1,d1]; 
x2 = [a2,b2,c2,d2]; 
% x2p = [a2p,b2p,c2p,d2p];
x3 = [a3,b3,c3,d3]; 
x4 = [a4,b4,c4,d4]; 

% myfit1_BN(x0,G0)
% myfit1_BN(x1,G1)
% myfit1_BN(x2,G2)
% myfit1_BN(x2p,G2p)
% myfit1_BN(x3,G3)
myfit1_BN(x4,G4)      
      
% figure1 = figure;
% axes1 = axes('Parent',figure1,...
%     'Position',[0.13 0.680241327300151 0.77 0.244758672699849]);
% hold(axes1,'on');
% plot(d',G0,'Parent',axes1,'DisplayName','G0','Marker','o','LineWidth',2);
% % ylabel('E (eV)');
% 
% title('t_{BpBp} (eV)');
% 
% box(axes1,'on');
% set(axes1,'FontSize',20,'FontWeight','bold','XTick',zeros(1,0));
% legend1 = legend(axes1,'show');
% set(legend1,...
%     'Position',[0.771028961154493 0.750303835144163 0.0970149253731343 0.0392156862745098]);
% axes2 = axes('Parent',figure1,...
%     'Position',[0.13 0.11 0.775 0.515367647058824]);
% hold(axes2,'on');
% plot(d',G1,'o-','Displayname','G1','linewidth',2)
% plot(d',G2,'o-','Displayname','G2','linewidth',2)
% plot(d',G2p,'o-','Displayname','G2*','linewidth',2)
% plot(d',G3,'o-','Displayname','G3','linewidth',2)
% plot(d',G4,'o-','Displayname','G4','linewidth',2)
% % ylabel('E (eV)','FontWeight','bold');
% xlabel('d A^o','FontWeight','bold');
% box(axes2,'on');
% set(axes2,'FontSize',20,'FontWeight','bold');
% legend2 = legend(axes2,'show');
% set(legend2,...
%     'Position',[0.772164426213304 0.20594807893233 0.107462686567164 0.174962292609351]);
% F1=550;F2=600;F3=550;F4=600;
% figure1.Position=[F1 F2 F3 F4];
%}

%% t_ApB
% {
G0 = ApB_d(1,:)';
G1 = ApB_d(2,:)';
G2 = ApB_d(3,:)';
G2p = ApB_d(4,:)';
G3 = ApB_d(5,:)';
G4 = ApB_d(6,:)';

func=@(x,a0,b0,c0,d0) a0.*exp(b0.*x) + c0.*exp(d0*x);

%        a0 =    -0.04005 ;
%        b0 =      0.3796 ;
%        c0 =       42.64 ;
%        d0 =      -1.359 ;
       
       a0 =     -0.0418 ;%
       b0 =      0.3708 ;%
       c0 =       42.41 ;%
       d0 =      -1.357 ;%

%        a1 =   -0.007954;
%        b1 =      0.6208;
%        c1 =           0;
%        d1 =      -38.65;
%        a1 =   -0.005896;
%        b1 =      0.7141;
%        c1 =       4.225;
%        d1 =      -8.981;
       
       a1 =   -0.006561;
       b1 =      0.6784;
       c1 =       4.225;
       d1 =      -8.981;

       a2 =     -0.5243;
       b2 =     -0.5917;
       c2 =       5.096;
       d2 =      -1.395;
       
%        a2 =   -0.006561;
%        b2 =      0.6784;
%        c2 =       4.225;
%        d2 =      -8.981;
       
       a2p =      0.4423;
       b2p =     -0.2127;
       c2p =     -0.2044;
       d2p =     0.04448;
       
%        a2p =      0.4334;
%        b2p =     -0.2172;
%        c2p =     -0.1945;
%        d2p =     0.04846;
       a2p =      0.4334 ;% (-1.224e+04, 1.224e+04)
       b2p =     -0.2173 ;% (-6679, 6678)
       c2p =     -0.1944 ;% (-1.304e+04, 1.304e+04)
       d2p =     0.04843 ;% (-6196, 6196)

       a3 =      -0.778;
       b3 =     -0.6396;
       c3 =        1.66;
       d3 =     -0.9052;
       
       a3 =     -0.7499 ;
       b3 =     -0.6357 ;
       c3 =       1.633 ;
       d3 =     -0.9094 ;
       
       a4 =      0.0155;
       b4 =      -10.72;
       c4 =    -0.01329;
       d4 =     -0.4014;
       a4 =      0.0155;
       b4 =      -10.72;
       c4 =    -0.01439 ;
       d4 =     -0.4246 ;


%figure
%hold on
% plot(d,G0,'Displayname','G0','Linewidth',3)
% plot(d,func(d,a0,b0,c0,d0),'--ro','Displayname','G0-fit','Markersize',5,'Linewidth',1.5)
% plot(d,G1,'Displayname','G1','Linewidth',3)
% plot(d,func(d,a1,b1,c1,d1),'--ro','Displayname','G1-fit','Markersize',5,'Linewidth',1.5)
% plot(d,G2,'Displayname','G2','Linewidth',3)
% plot(d,func(d,a2,b2,c2,d2),'--ro','Displayname','G2-fit','Markersize',5,'Linewidth',1.5)
% plot(d,G2p,'Displayname','G2p','Linewidth',3)
% plot(d,func(d,a2p,b2p,c2p,d2p),'--ro','Displayname','G2p-fit','Markersize',5,'Linewidth',1.5)
% plot(d,G3,'Displayname','G3','Linewidth',3)
% plot(d,func(d,a3,b3,c3,d3),'--ro','Displayname','G3-fit','Markersize',5,'Linewidth',1.5)
%plot(d,G4,'Displayname','G4','Linewidth',3)
%plot(d,func(d,a4,b4,c4,d4),'--ro','Displayname','G4-fit','Markersize',5,'Linewidth',1.5)
%legend('show')
       
x = d;
x0 = [a0,b0,c0,d0]; 
x1 = [a1,b1,c1,d1]; 
x2 = [a2,b2,c2,d2]; 
x2p = [a2p,b2p,c2p,d2p];
x3 = [a3,b3,c3,d3]; 
x4 = [a4,b4,c4,d4]; 

% myfit1_BN(x0,G0)
%myfit1_BN(x1,G1)
%myfit1_BN(x2,G2)
%  myfit1_BN(x2p,G2p)
% myfit1_BN(x3,G3)
%myfit1_BN(x4,G4)   
       
% d = [3.1,3.2,3.3,3.4]';
dd = linspace(3.1,3.5,100);

figure1 = figure;
axes1 = axes('Parent',figure1,...
    'Position',[0.13 0.680241327300151 0.77 0.244758672699849]);
hold(axes1,'on');
% plot(d',G0,'Parent',axes1,'DisplayName','G0','Marker','o','LineWidth',2);
h1 = plot(d',G0,'o','Displayname','G0','linewidth',2,'Markersize',10);
plot(dd,func(dd,a0,b0,c0,d0),'-','Markersize',5,'Linewidth',3)
% ylabel('E (eV)');

title('t_{ApB} (eV)');

box(axes1,'on');
set(axes1,'FontSize',20,'FontWeight','bold','XTick',zeros(1,0));
legend1 = legend(axes1,'show',h1);
set(legend1,...
    'Position',[0.771028961154493 0.750303835144163 0.0970149253731343 0.0392156862745098]);
axes2 = axes('Parent',figure1,...
    'Position',[0.13 0.11 0.775 0.515367647058824]);
hold(axes2,'on');
h2 = plot(d',G1,'o','Displayname','G1','linewidth',2,'Markersize',10);
plot(dd,func(dd,a1,b1,c1,d1),'-','Markersize',5,'Linewidth',3)
h3 = plot(d',G2,'o','Displayname','G2','linewidth',2,'Markersize',10);
plot(dd,func(dd,a2,b2,c2,d2),'-','Markersize',5,'Linewidth',3)
h4 = plot(d',G2p,'o','Displayname','G2*','linewidth',2,'Markersize',10);
plot(dd,func(dd,a2p,b2p,c2p,d2p),'-','Markersize',5,'Linewidth',3)
h5 = plot(d',G3,'o','Displayname','G3','linewidth',2,'Markersize',10);
plot(dd,func(dd,a3,b3,c3,d3),'-','Markersize',5,'Linewidth',3)
h6 = plot(d',G4,'o','Displayname','G4','linewidth',2,'Markersize',10);
plot(dd,func(dd,a4,b4,c4,d4),'-','Markersize',5,'Linewidth',3)
% ylabel('E (eV)','FontWeight','bold');
xlabel('d A^o','FontWeight','bold');
box(axes2,'on');
set(axes2,'FontSize',20,'FontWeight','bold');
legend2 = legend(axes2,'show',[h2,h3,h4,h5,h6]);


set(legend2,...
    'Position',[0.772164426213304 0.20594807893233 0.107462686567164 0.174962292609351]);
F1=550;F2=600;F3=550;F4=600;
figure1.Position=[F1 F2 F3 F4];
%}

%% t_AB
% {
F1 = AB_d(1,:)';F1(3) = -2.6846;F1(2) = -2.6858; 
F2 = AB_d(2,:)';
F3 = AB_d(3,:)';
F4 = AB_d(4,:)';

func=@(x,a0,b0,c0,d0) a0.*exp(b0.*x) + c0.*exp(d0*x);
       a1 =      -2.805;
       b1 =   -0.003352;
       c1 =     0.09248;
       d1 =    -0.01006;

       a2 =       -2.89;
       b2 =      -7.455;
       c2 =     -0.3349;
       d2 =     -0.1262;

       a3 =     0.03819;
       b3 =      0.2253;
       c3 =      -1.991;
       d3 =       -7.39;

       a4 =      0.2825;
       b4 =    0.006668;
       c4 =      -0.318;
       d4 =   -0.006668;


       
% figure
% hold on

% plot(d,F1,'Displayname','F1','Linewidth',3)
% plot(d,func(d,a1,b1,c1,d1),'--r','Displayname','F1-fit','Markersize',5,'Linewidth',1.5)
% 
% plot(d,F2,'Displayname','F2','Linewidth',3)
% plot(d,func(d,a2,b2,c2,d2),'--r','Displayname','F2-fit','Markersize',5,'Linewidth',1.5)

% plot(d,F3,'Displayname','F3','Linewidth',3)
% plot(d,func(d,a3,b3,c3,d3),'--r','Displayname','F3-fit','Markersize',5,'Linewidth',1.5)
%  
% plot(d,F4,'Displayname','F4','Linewidth',3)
% plot(d,func(d,a4,b4,c4,d4),'--r','Displayname','F4-fit','Markersize',5,'Linewidth',1.5)
% 
% legend('show')

x = d;
x1 = [a1,b1,c1,d1]; 
x2 = [a2,b2,c2,d2]; 
x3 = [a3,b3,c3,d3]; 
x4 = [a4,b4,c4,d4]; 

% myfit1_BN(x1,F1)
% myfit1_BN(x2,F2)
% myfit1_BN(x3,F3)
% myfit1_BN(x4,F4)
       
figure1 = figure;
axes1 = axes('Parent',figure1,...
    'Position',[0.13 0.680241327300151 0.77 0.244758672699849]);
hold(axes1,'on');
% plot(d',F1,'Parent',axes1,'DisplayName','F1','Marker','o','LineWidth',2);
h1 = plot(d',F1,'o','Displayname','F1','linewidth',2);
plot(d,func(d,a1,b1,c1,d1),'-','Markersize',5,'Linewidth',1.5)
legend1 = legend(axes1,'show',h1);
% ylabel('E (eV)');

title('t_{AB} (eV)');

box(axes1,'on');
set(axes1,'FontSize',20,'FontWeight','bold','XTick',zeros(1,0));
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.771028961154493 0.750303835144163 0.0970149253731343 0.0392156862745098]);
axes2 = axes('Parent',figure1,...
    'Position',[0.13 0.11 0.775 0.515367647058824]);
hold(axes2,'on');
h2 = plot(d',F2,'o','Displayname','F2','linewidth',2)
plot(d,func(d,a2,b2,c2,d2),'-','Markersize',5,'Linewidth',1.5)
h3 = plot(d',F3,'o','Displayname','F3','linewidth',2)
plot(d,func(d,a3,b3,c3,d3),'-','Markersize',5,'Linewidth',1.5)
h4 = plot(d',F4,'o','Displayname','F4','linewidth',2)
plot(d,func(d,a4,b4,c4,d4),'-','Markersize',5,'Linewidth',1.5)
% ylabel('E (eV)','FontWeight','bold');
xlabel('d A^o','FontWeight','bold');
box(axes2,'on');
set(axes2,'FontSize',20,'FontWeight','bold');
legend2 = legend(axes2,'show',[h2,h3,h4]);
set(legend2,...
    'Position',[0.772164426213304 0.20594807893233 0.107462686567164 0.174962292609351]);
Ff1=550;Ff2=600;Ff3=550;Ff4=600;
figure1.Position=[Ff1 Ff2 Ff3 Ff4];
%}

%% t_ApBp
%{
F1 = ApBp_d(1,:)';
F2 = ApBp_d(2,:)';
F3 = ApBp_d(3,:)';
F4 = ApBp_d(4,:)';

func=@(x,a0,b0,c0,d0) a0.*exp(b0.*x) + c0.*exp(d0*x);

       a1 =    -0.02221;
       b1 =       -4.32;
       c1 =      -2.688;
       d1 =    0.001259;

       a2 =       0.787;
       b2 =      -0.273;
       c2 =     -0.9677;
       d2 =      -0.184;

       a3 =     0.03648;
       b3 =      0.2405;
       c3 =      -4.233;
       d3 =      -5.304;
       
       a4 =     -0.5826;
       b4 =      0.4629;
       c4 =      0.5664;
       d4 =       0.469;

   
figure
hold on

% plot(d,F1,'Displayname','F1','Linewidth',3)
% plot(d,func(d,a1,b1,c1,d1),'--ro','Displayname','F1-fit','Markersize',5,'Linewidth',1.5)

% plot(d,F2,'Displayname','F2','Linewidth',3)
% plot(d,func(d,a2,b2,c2,d2),'--ro','Displayname','F2-fit','Markersize',5,'Linewidth',1.5)

plot(d,F3,'Displayname','F3','Linewidth',3)
plot(d,func(d,a3,b3,c3,d3),'--ro','Displayname','F3-fit','Markersize',5,'Linewidth',1.5)
 
% plot(d,F4,'Displayname','F4','Linewidth',3)
% plot(d,func(d,a4,b4,c4,d4),'--ro','Displayname','F4-fit','Markersize',5,'Linewidth',1.5)

legend('show')

x = d;
x1 = [a1,b1,c1,d1]; 
x2 = [a2,b2,c2,d2]; 
x3 = [a3,b3,c3,d3]; 
x4 = [a4,b4,c4,d4]; 

% myfit1_BN(x1,F1)
% myfit1_BN(x2,F2)
myfit1_BN(x3,F3)
% myfit1_BN(x4,F4)

       
% figure1 = figure;
% axes1 = axes('Parent',figure1,...
%     'Position',[0.13 0.680241327300151 0.77 0.244758672699849]);
% hold(axes1,'on');
% plot(d',F1,'Parent',axes1,'DisplayName','F1','Marker','o','LineWidth',2);
% % ylabel('E (eV)');
% 
% title('t_{ApBp} (eV)');
% 
% box(axes1,'on');
% set(axes1,'FontSize',20,'FontWeight','bold','XTick',zeros(1,0));
% legend1 = legend(axes1,'show');
% set(legend1,...
%     'Position',[0.771028961154493 0.750303835144163 0.0970149253731343 0.0392156862745098]);
% axes2 = axes('Parent',figure1,...
%     'Position',[0.13 0.11 0.775 0.515367647058824]);
% hold(axes2,'on');
% plot(d',F2,'o-','Displayname','F2','linewidth',2)
% plot(d',F3,'o-','Displayname','F3','linewidth',2)
% plot(d',F4,'o-','Displayname','F4','linewidth',2)
% % ylabel('E (eV)','FontWeight','bold');
% xlabel('d A^o','FontWeight','bold');
% box(axes2,'on');
% set(axes2,'FontSize',20,'FontWeight','bold');
% legend2 = legend(axes2,'show');
% set(legend2,...
%     'Position',[0.772164426213304 0.20594807893233 0.107462686567164 0.174962292609351]);
% Ff1=550;Ff2=600;Ff3=550;Ff4=600;
% figure1.Position=[Ff1 Ff2 Ff3 Ff4];
%}

%% t_AAp
%{
F1 = AAp_d(1,:)';
F2 = AAp_d(2,:)';
F3 = AAp_d(3,:)';
F4 = AAp_d(4,:)';F4(4) = 0.0216;


func=@(x,a0,b0,c0,d0) a0.*exp(b0.*x) + c0.*exp(d0*x);

       a1 =     -0.5973;
       b1 =      0.9163;
       c1 =      0.7005;
       d1 =      0.8792;

       a2 =     -0.6613;
       b2 =      0.8112;
       c2 =      0.6552;
       d2 =      0.8161;

       a3 =   -0.008976;
       b3 =      0.6663;
       c3 =      0.0128;
       d3 =      0.7016;
       
       a4 =    0.009705;
       b4 =      0.2387;
       c4 =      0.7201;
       d4 =      -3.841;

figure
hold on

% plot(d,F1,'Displayname','F1','Linewidth',3)
% plot(d,func(d,a1,b1,c1,d1),'--ro','Displayname','F1-fit','Markersize',5,'Linewidth',1.5)

% plot(d,F2,'Displayname','F2','Linewidth',3)
% plot(d,func(d,a2,b2,c2,d2),'--ro','Displayname','F2-fit','Markersize',5,'Linewidth',1.5)

% plot(d,F3,'Displayname','F3','Linewidth',3)
% plot(d,func(d,a3,b3,c3,d3),'--ro','Displayname','F3-fit','Markersize',5,'Linewidth',1.5)

plot(d,F4,'Displayname','F4','Linewidth',3)
plot(d,func(d,a4,b4,c4,d4),'--ro','Displayname','F4-fit','Markersize',5,'Linewidth',1.5)

legend('show')

x = d;
x1 = [a1,b1,c1,d1]; 
x2 = [a2,b2,c2,d2]; 
x3 = [a3,b3,c3,d3]; 
x4 = [a4,b4,c4,d4]; 

% myfit1_BN(x1,F1)
% myfit1_BN(x2,F2)
% myfit1_BN(x3,F3)
myfit1_BN(x4,F4)
      
% figure1 = figure;
% axes1 = axes('Parent',figure1,...
%     'Position',[0.13 0.680241327300151 0.77 0.244758672699849]);
% hold(axes1,'on');
% plot(d',F1,'Parent',axes1,'DisplayName','F1','Marker','o','LineWidth',2);
% % ylabel('E (eV)');
% 
% title('t_{AAp} (eV)');
% 
% box(axes1,'on');
% set(axes1,'FontSize',20,'FontWeight','bold','XTick',zeros(1,0));
% legend1 = legend(axes1,'show');
% set(legend1,...
%     'Position',[0.771028961154493 0.750303835144163 0.0970149253731343 0.0392156862745098]);
% axes2 = axes('Parent',figure1,...
%     'Position',[0.13 0.11 0.775 0.515367647058824]);
% hold(axes2,'on');
% plot(d',F2,'o-','Displayname','F2','linewidth',2)
% plot(d',F3,'o-','Displayname','F3','linewidth',2)
% plot(d',F4,'o-','Displayname','F4','linewidth',2)
% % ylabel('E (eV)','FontWeight','bold');
% xlabel('d A^o','FontWeight','bold');
% box(axes2,'on');
% set(axes2,'FontSize',20,'FontWeight','bold');
% legend2 = legend(axes2,'show');
% set(legend2,...
%     'Position',[0.772164426213304 0.20594807893233 0.107462686567164 0.174962292609351]);
% Ff1=550;Ff2=600;Ff3=550;Ff4=600;
% figure1.Position=[Ff1 Ff2 Ff3 Ff4];
%}

%% t_BBp
%{
F1 = BBp_d(1,:)';
F2 = BBp_d(2,:)';
F3 = BBp_d(3,:)';
F4 = BBp_d(4,:)';


func=@(x,a0,b0,c0,d0) a0.*exp(b0.*x) + c0.*exp(d0*x);

       a1 =     -0.1576;
       b1 =     -0.3161;
       c1 =       1.357;
       d1 =      -1.068;
       
       a2 =   0.0008165;
       b2 =      0.8209;
       c2 =     -0.1054;
       d2 =     -0.1697;

       a3 =      0.1372;
       b3 =      0.2613;
       c3 =     -0.1557;
       d3 =      0.2034;

       a4 =    0.002648;
       b4 =      0.3089;
       c4 =      -1.814;
       d4 =      -5.524;


figure
hold on

% plot(d,F1,'Displayname','F1','Linewidth',3)
% plot(d,func(d,a1,b1,c1,d1),'--ro','Displayname','F1-fit','Markersize',5,'Linewidth',1.5)

% plot(d,F2,'Displayname','F2','Linewidth',3)
% plot(d,func(d,a2,b2,c2,d2),'--ro','Displayname','F2-fit','Markersize',5,'Linewidth',1.5)

% plot(d,F3,'Displayname','F3','Linewidth',3)
% plot(d,func(d,a3,b3,c3,d3),'--ro','Displayname','F3-fit','Markersize',5,'Linewidth',1.5)

plot(d,F4,'Displayname','F4','Linewidth',3)
plot(d,func(d,a4,b4,c4,d4),'--ro','Displayname','F4-fit','Markersize',5,'Linewidth',1.5)

legend('show')

x = d;
x1 = [a1,b1,c1,d1]; 
x2 = [a2,b2,c2,d2]; 
x3 = [a3,b3,c3,d3]; 
x4 = [a4,b4,c4,d4]; 

% myfit1_BN(x1,F1)
% myfit1_BN(x2,F2)
% myfit1_BN(x3,F3)
myfit1_BN(x4,F4)
       
% figure1 = figure;
% axes1 = axes('Parent',figure1,...
%     'Position',[0.13 0.680241327300151 0.77 0.244758672699849]);
% hold(axes1,'on');
% plot(d',F1,'Parent',axes1,'DisplayName','F1','Marker','o','LineWidth',2);
% % ylabel('E (eV)');
% 
% title('t_{BBp} (eV)');
% 
% box(axes1,'on');
% set(axes1,'FontSize',20,'FontWeight','bold','XTick',zeros(1,0));
% legend1 = legend(axes1,'show');
% set(legend1,...
%     'Position',[0.771028961154493 0.750303835144163 0.0970149253731343 0.0392156862745098]);
% axes2 = axes('Parent',figure1,...
%     'Position',[0.13 0.11 0.775 0.515367647058824]);
% hold(axes2,'on');
% plot(d',F2,'o-','Displayname','F2','linewidth',2)
% plot(d',F3,'o-','Displayname','F3','linewidth',2)
% plot(d',F4,'o-','Displayname','F4','linewidth',2)
% % ylabel('E (eV)','FontWeight','bold');
% xlabel('d A^o','FontWeight','bold');
% box(axes2,'on');
% set(axes2,'FontSize',20,'FontWeight','bold');
% legend2 = legend(axes2,'show');
% set(legend2,...
%     'Position',[0.772164426213304 0.20594807893233 0.107462686567164 0.174962292609351]);
% Ff1=550;Ff2=600;Ff3=550;Ff4=600;
% figure1.Position=[Ff1 Ff2 Ff3 Ff4];
%}

%% t_ABp
%{
F1 = ABp_d(1,:)';
F2 = ABp_d(2,:)';
F3 = ABp_d(3,:)';
F4 = ABp_d(4,:)';

func=@(x,a0,b0,c0,d0) a0.*exp(b0.*x) + c0.*exp(d0*x);
      
       a1 =    -0.03866;
       b1 =      0.3746;
       c1 =       9.974;
       d1 =      -1.113;

       a2 =      -0.268;
       b2 =     -0.2156;
       c2 =     -0.9209;
       d2 =      -6.499;

       a3 =     -0.5453;
       b3 =     -0.2354;
       c3 =       1.415;
       d3 =     -0.5496;
       
       a4 =      -0.639;
       b4 =      -2.757;
       c4 =       4.099;
       d4 =      -3.934;


figure
hold on

plot(d,F1,'Displayname','F1','Linewidth',3)
plot(d,func(d,a1,b1,c1,d1),'--ro','Displayname','F1-fit','Markersize',5,'Linewidth',1.5)

plot(d,F2,'Displayname','F2','Linewidth',3)
plot(d,func(d,a2,b2,c2,d2),'--ro','Displayname','F2-fit','Markersize',5,'Linewidth',1.5)

plot(d,F3,'Displayname','F3','Linewidth',3)
plot(d,func(d,a3,b3,c3,d3),'--ro','Displayname','F3-fit','Markersize',5,'Linewidth',1.5)

plot(d,F4,'Displayname','F4','Linewidth',3)
plot(d,func(d,a4,b4,c4,d4),'--ro','Displayname','F4-fit','Markersize',5,'Linewidth',1.5)

legend('show')

x = d;
x1 = [a1,b1,c1,d1]; 
x2 = [a2,b2,c2,d2]; 
x3 = [a3,b3,c3,d3]; 
x4 = [a4,b4,c4,d4]; 

myfit1_BN(x1,F1)
% myfit1_BN(x2,F2)
% myfit1_BN(x3,F3)
%myfit1_BN(x4,F4)

figure1 = figure;
axes1 = axes('Parent',figure1,...
    'Position',[0.13 0.680241327300151 0.77 0.244758672699849]);
hold(axes1,'on');
% plot(d',F1,'Parent',axes1,'DisplayName','F1','Marker','o','LineWidth',2);
h1 = plot(d',F1,'o','Displayname','F1','linewidth',2);
plot(d,func(d,a1,b1,c1,d1),'-','Linewidth',1.5)
% ylabel('E (eV)');

title('t_{ABp} (eV)');

box(axes1,'on');
set(axes1,'FontSize',20,'FontWeight','bold','XTick',zeros(1,0));
legend1 = legend(axes1,'show',h1);
set(legend1,...
    'Position',[0.771028961154493 0.750303835144163 0.0970149253731343 0.0392156862745098]);
axes2 = axes('Parent',figure1,...
    'Position',[0.13 0.11 0.775 0.515367647058824]);
hold(axes2,'on');
h2 = plot(d',F2,'o','Displayname','F2','linewidth',2);
plot(d,func(d,a2,b2,c2,d2),'-','Linewidth',1.5)
h3 = plot(d',F3,'o','Displayname','F3','linewidth',2);
plot(d,func(d,a3,b3,c3,d3),'-','Linewidth',1.5)
h4 = plot(d',F4,'o','Displayname','F4','linewidth',2);
plot(d,func(d,a4,b4,c4,d4),'-','Linewidth',1.5)
% ylabel('E (eV)','FontWeight','bold');
xlabel('d A^o','FontWeight','bold');
box(axes2,'on');
set(axes2,'FontSize',20,'FontWeight','bold');
legend2 = legend(axes2,'show',[h2,h3,h4]);
set(legend2,...
    'Position',[0.772164426213304 0.20594807893233 0.107462686567164 0.174962292609351]);
Ff1=550;Ff2=600;Ff3=550;Ff4=600;
figure1.Position=[Ff1 Ff2 Ff3 Ff4];
%}







