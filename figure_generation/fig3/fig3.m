%This script plots figure 3

clear;clc

newfigure(3.42/2, (2/3)*3.42/3*(3.5/2));

gap = [0.17*1.35,0.11];
marg_h = [0.15,0.07];

Fig3ax = tight_subplot(2,2,gap,marg_h,[0.175,0.05]);

axes(Fig3ax(1));
axis off
text(-0.45,1.1,'\textbf{A}','Interpreter','latex','Units','normalized','FontSize',4)
axes(Fig3ax(1))
pos = get(gca,'Position');
pos(2) = 0.97*pos(2);
set(gca,'Position',pos);

fontsize = 4;
linewidth = 0.4;

load('fig3_data.mat')


%% Plot pdfs of varying number of co-limiting metabolites
axes(Fig3ax(2))
pos = get(gca,'Position');
pos(2) = 0.97*pos(2);
set(gca,'Position',pos);


hold on
colors= parula(n1);
ind = 5;
for i = 1:n1
    plot(x,pdf{i,4},'-','LineWidth',linewidth, 'Color',colors(i,:),'DisplayName',['n=',num2str(i)]);
end
xlim([0,2*prot_mean]);
ylim([0,0.12]);
yticks([])
xticks([])
xlabel('Growth rate','FontSize',fontsize,'Interpreter','latex')
ylabel('Probability','FontSize',fontsize,'Interpreter','latex')
text(-0.15,1.1,'\textbf{B}','Interpreter','latex','Units','normalized','FontSize',4)

colormap(Fig3ax(4),parula);
c = colorbar('TickLabels',{'$1$' '$20$'},'Ticks',[0 1],'TickLabelInterpreter','latex');
c.Location = 'North';
c.Position(1) = 0.67;
c.Position(2) = 0.945;
c.Position(3) = 0.25;
c.Position(4) = 0.05;
c.Ticks=[];
c.AxisLocationMode='manual';
c.AxisLocation = 'in';
c.Label.String = '$\#$ of metabolites';
c.Label.Interpreter = 'latex';
c.Label.Position = [0.5,0.94,0];
c.Label.Rotation = 0;
c.Label.Color='w';
c.Label.FontSize = 4;
text(0.1,1.065+0.15,'1','Interpreter','latex','Units','normalized','FontSize',4)
text(0.94,1.065+0.15,'20','Interpreter','latex','Units','normalized','FontSize',4)


%% Plot growth rates for different CVs

axes(Fig3ax(3))
hold on
colors= cool(length(cvs));
for i = 1:length(cvs)
    if i == 4
        linetype = ':';
        mult = 1.7;
    else
        linetype = '-';
        mult = 1;
    end
    plot(1:n1,dim_means(:,i)/prot_mean,linetype,'LineWidth',mult*linewidth, 'Color',colors(i,:),'DisplayName',['CV=',num2str(cvs(i))]);
end

ylim([0,1.1]);
yticks([0,0.5,1]);
xlim([1,20])
xticks([1,10,20])
xticklabels({'1','10','20'})
yticklabels({'0','0.5','1'})
ylabel({'Mean', 'growth rate'},'FontSize',fontsize,'Interpreter','latex')
xlabel('Number of metabolites','FontSize',fontsize,'Interpreter','latex')
text(-0.45,1.1,'\textbf{C}','Interpreter','latex','Units','normalized','FontSize',4)
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',4)

colormap(Fig3ax(3),cool);
c = colorbar('TickLabels',{'$1$' '$20$'},'Ticks',[0 1],'TickLabelInterpreter','latex');
c.Location = 'North';
c.Position(1) = 0.2455;
c.Position(2) = 0.459;
c.Position(3) = 0.25;
c.Position(4) = 0.05;
c.Ticks=[];
c.AxisLocationMode='manual';
c.AxisLocation = 'in';
c.Label.String = 'Metab. CV';
c.Label.Interpreter = 'latex';
c.Label.Position = [0.5,1,0];
c.Label.Rotation = 0;
c.Label.Color='w';
c.Label.FontSize = 4;
text(0.02,1.065+0.15,'0.05','Interpreter','latex','Units','normalized','FontSize',4)
text(0.99,1.065+0.15,'1','Interpreter','latex','Units','normalized','FontSize',4)
ax = gca;
ax.Clipping = 'off';

%% Plot growth rates for different metabolite correlations

axes(Fig3ax(4))
hold on
colors= summer(length(rho_vec));
for i = 1:length(rho_vec)
    if i == 8
        linetype = ':';
        mult = 1.7;
    else
        linetype = '-';
        mult = 1;
    end
    plot(n_vec,results(:,i)/prot_mean,linetype,'LineWidth',mult*linewidth, 'Color',colors(i,:),'DisplayName',['cov=',num2str(rho_vec(i))]);
end

ylim([0,1.1]);
xlim([1,20])
xticks([1,10,20])
xticklabels({'1','10','20'})
xlabel('Number of metabolites','FontSize',fontsize,'Interpreter','latex')
text(-0.15,1.1,'\textbf{D}','Interpreter','latex','Units','normalized','FontSize',4)
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',4)

colormap(Fig3ax(4),summer);
c = colorbar('TickLabels',{'$1$' '$20$'},'Ticks',[0 1],'TickLabelInterpreter','latex');
c.Location = 'North';
c.Position(1) = 0.67;
c.Position(2) = 0.459;
c.Position(3) = 0.25;
c.Position(4) = 0.05;
c.Ticks=[];
c.AxisLocationMode='manual';
c.AxisLocation = 'in';
c.Label.String = 'Metab. corr.';
c.Label.Interpreter = 'latex';
c.Label.Position = [0.5,1,0];
c.Label.Rotation = 0;
c.Label.Color='w';
c.Label.FontSize = 4;
text(0.1,1.065+0.15,'0','Interpreter','latex','Units','normalized','FontSize',4)
text(0.94,1.065+0.15,'0.9','Interpreter','latex','Units','normalized','FontSize',4)
ax = gca;
ax.Clipping = 'off';

print(gcf,'-dpng','fig3.png','-r1200');
close all
