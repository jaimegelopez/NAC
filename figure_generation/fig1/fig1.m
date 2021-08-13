clear;clc

newfigure(3.42/2, (4/3)*3.42/3*(3.5/2));

gap = [0.12,0.11];
marg_h = [0.1,0.07];

Fig1ax = tight_subplot(4,2,gap,marg_h,[0.18,0.05]);

for i = 1:4
    axes(Fig1ax(i));
    axis off
end

axes(Fig1ax(3));
text(-0.25,0,'\textbf{A}','Interpreter','latex','Units','normalized','FontSize',4)
text(-0.25,-0.5,'\textbf{B}','Interpreter','latex','Units','normalized','FontSize',4)


fontsize = 4;
linewidth = 0.4;

%% Look at typical behavior for good and bad regulation
seed=666;
cutoff = 1;
fig1_regulation(Fig1ax,linewidth,fontsize,seed,cutoff);

%% Look at growth as a function of parameters
fig1_par_sweep(Fig1ax,linewidth,fontsize)

%% Save

print(gcf, '-dpng','fig1.png','-r1500');
close all
