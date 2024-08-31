function [fitted_curve,gof] = myfit_BN(x0,y)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
dis = [3.1,3.2,3.3,3.4,3.5]';
x = dis;

fitfun0 = fittype( @(a,b,c,d,x) a.*exp(b.*x) + c.*exp(d*x) );   % V0
[fitted_curve,gof] = fit(x,y,fitfun0,'StartPoint',x0);
coeffvals = coeffvalues(fitted_curve);
figure
plot(x, y, 'r-+','Displayname','data')
hold on
plot(x,fitted_curve(x),'Displayname','fit')
legend('show')
hold off

% outputArg1 = inputArg1;
% outputArg2 = inputArg2;
end

