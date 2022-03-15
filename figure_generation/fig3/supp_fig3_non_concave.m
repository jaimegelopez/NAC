clear;clc

%This script plots examples of functions that are not everywhere concave
%but still have negative Jensen gaps

%Define functions
f1 = @(x) 1./(1+exp(-5*x + 2.1)) - 1/(1 + exp(2.1));
f2 = @(x) 0.02*sin(100*x) + x./(0.3 + x);

%Simulate to estimate Jensen gaps
a =  0;
b = 1; 
x = linspace(0,1,1000);
mu = (a+b)/2;
X = (b - a)*rand(1e5,1) + a;
true_gap1 = mean(f1(X)) - f1(mean(X));
true_gap2 = mean(f2(X)) - f2(mean(X));

%Plot
fontsize = 4;
linewidth = 0.4;
newfigure(0.6*3.42/2, (1/3)*3.42/3*(3.5/2));
hold on
plot(x,f1(x),'r-','LineWidth',linewidth);
plot(x,f2(x),'k-','LineWidth',linewidth);
xlabel('Metabolite level','Interpreter','latex')
ylabel('Growth rate','Interpreter','latex')
set(gca,'FontSize',fontsize);
set(gca,'TickLabelInterpreter','latex')

print(gcf, '-dpng','supp_fig3_non_concave.png','-r1500');
