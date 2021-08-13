clear;clc

bdata = 'aa_analysis_data.mat';

newfigure(3.42/2, (1/3)*3.42/3*(3.5/2));
hold on

%Load data
load(bdata);

%Select best loss estimation, direct > logP regression
combo_vec = data.aa_loss;
combo_vec(isnan(combo_vec)) = data.aa_est_loss(isnan(combo_vec));
data.best_loss = combo_vec;
bvec = data{:,'best_loss'};
ind = ~isnan(bvec);
bvec = bvec(ind);
gene_names = data.aa_triple(ind);

%Plot
bar(bvec,'FaceColor',[0.5,0.5,0.5])
set(gca,'YScale','log')
xticks(1:length(bvec))
xticklabels(gene_names)
ylim([1e-6,1e-2])
yticks([1e-6,1e-4,1e-2])
ylabel('Frac. prod. lost','Interpreter','latex')
pos = get(gcf,'Position');
pos(3) = (0.2+1.1*(length(bvec)/20))*pos(3);
set(gcf,'Position',pos);
xtickangle(60)
set(gca,'YMinorTick','Off')
set(gca,'FontSize',4)
set(gca,'TickLabelInterpreter','latex')
yticklabels({'$10^{-6}$','$10^{-4}$','$10^{-2}$'})
print(gcf,'-dpng','supp_fig_budget_loss.png','-r1200');

