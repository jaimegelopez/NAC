function fig1_par_sweep(Fig1ax,linewidth,fontsize)

load('large_single_cell_sweep.mat')

max_growth = 25*100/(1+1);

colors = [0 0 0;...
    230 159 0;...
    86 180 233;...
    0 158 115; ...
    240 228 66;...
    213 94 0;...
    204 121 167]/255;

%% Plot

par_cell = {'P','beta'};
label_cell = {{'Cell permeability, $P$'},'Burst size, $\beta$'};
inds = [1,5];
n1 = 40;
n2 = 2;
n3 = [4,3];

for i = 1:n2
    axes(Fig1ax(i+6))
    hold on
    for j = 1:n3(i)
        truei = inds(i)+j-1;
        curr_index = (truei-1)*n1 + 1:truei*n1;
        curr_mean = results_table.mean_growth(curr_index);
        curr_par = results_table{curr_index,par_cell{i}};
        
        plot(curr_par,curr_mean./max_growth,'k-','LineWidth',linewidth,'Color',colors(j+4*(i-1),:));
    end
    
    xlabel(label_cell{i},'Interpreter','latex')
    set(gca,'XScale','log')
    set(gca,'FontSize',fontsize);
    xlim([1e-1,max(curr_par)])
    ylim([-0.03,1])
    yticks([0,0.5,1])
    
    if i == 1
        ylabel('Growth rate','Interpreter','latex')
        yticklabels({'0','0.5','1'})
    end
    
    if i == 1
        text(-0.25,1.25,'\textbf{E}','Interpreter','latex','Units','normalized','FontSize',4)
    else
        text(-0.1,1.25,'\textbf{F}','Interpreter','latex','Units','normalized','FontSize',4)
    end
    
    xticks([1e-1,1e1,1e3])
    xticklabels({'$10^{-1}$','$10^1$','$10^3$'})
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'XMinorTick','Off')
    
    leg_top = 1;
    leg_left = 30;
    mylines = {'-','-','-','-'};
    if i == 1
        labels = {'$\beta = 2$','$\beta = 10$','$\beta = 20$','$\beta = 100$'};
    else
        labels = {'$P = 0$','$P = 1$','$P = 10$'};
    end
    line_length = 30;
    spacing_x = 0.3;
    curr_y = leg_top;
    for k = 1:length(labels)
        plot([leg_left,leg_left + line_length],[curr_y,curr_y],mylines{k},'Color',colors(k+4*(i-1),:))
        text(leg_left + line_length + spacing_x,curr_y,labels{k},'Interpreter','latex','FontSize',4)
        curr_y = curr_y-0.18;
    end
    
    
end

end

