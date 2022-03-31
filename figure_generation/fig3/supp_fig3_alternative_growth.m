%This script looks at mean growth behavior for an alternative growth
%function. 

clear;clc

%Define function and generate growth data

altfun = @(x,pstar) 1./((1/size(x,1)).*sum((x+pstar)./x,1));

prot_mean = (5.2*6.4);
pstar = prot_mean;
n1 = 20;
cvs = linspace(0.05,1,9);
cvs(4) = 0.4;

n2 = 50000;

for j = 1:length(cvs)
    k = (1/cvs(j))^2;
    theta = prot_mean/k;
    
    for i = 1:n1
        data = gamrnd(k,theta,i,n2);
        g(i,j) = mean(altfun(data,pstar));
        
    end
end

max_growth = altfun(prot_mean,pstar);

%% Plot growth rates for different CVs

newfigure(1/2*3.42/2, 1/2*(2/3)*3.42/3*(3.5/2));
fontsize = 4;
linewidth = 0.4;

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
    plot(1:n1,g(:,i)/max_growth,linetype,'LineWidth',mult*linewidth, 'Color',colors(i,:),'DisplayName',['CV=',num2str(cvs(i))]);
end

ylim([0,1.1]);
yticks([0,0.5,1]);
xlim([1,20])
xticks([1,10,20])
xticklabels({'1','10','20'})
yticklabels({'0','0.5','1'})
ylabel({'Mean', 'growth rate'},'FontSize',fontsize,'Interpreter','latex')
xlabel('Number of metabolites','FontSize',fontsize,'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',4)

colormap(cool);
c = colorbar('TickLabels',{'$1$' '$20$'},'Ticks',[0 1],'TickLabelInterpreter','latex');
c.Location = 'North';
c.Position(1) = 0.1755*2;
c.Position(2) = 0.445*2;
c.Position(3) = 0.25*2;
c.Position(4) = 0.05*2;
c.Ticks=[];
c.AxisLocationMode='manual';
c.AxisLocation = 'in';
c.Label.String = 'Metab. CV';
c.Label.Interpreter = 'latex';
c.Label.Position = [0.5,1,0];
c.Label.Rotation = 0;
c.Label.Color='w';
c.Label.FontSize = 4;
text_y = 1.125;
text_x_off = 0.08;
text(0.01-text_x_off,text_y,'0.05','Interpreter','latex','Units','normalized','FontSize',4)
text(1.02-text_x_off,text_y,'1','Interpreter','latex','Units','normalized','FontSize',4)
ax = gca;
ax.Clipping = 'off';

pos = get(gca,'Position');
pos(4) = 0.9*pos(4);
set(gca,'Position',pos)

print(gcf,'-dpng','supp_fig3_alternative_growth.png','-r1200');
