%This script generates figure 2

clear;clc

newfigure(3.42/2, (2.3/3)*3.42/3*(3.5/2));

gap = [0.23,0.18];
marg_h = [0.15,0.05];

Fig2ax = tight_subplot(2,2,gap,marg_h,[0.18,0.05]);

for i = 1:2
    axes(Fig2ax(i));
    axis off
end

fontsize = 4;
linewidth = 0.4;

text(0.5,0.5,'\textbf{A}','Interpreter','latex','Units','normalized','FontSize',4)

%% Parameter sweep
fig2_par_sweep(Fig2ax,linewidth,fontsize)


%% Save and close

print(gcf, '-dpng','fig2.png','-r1200');
close all