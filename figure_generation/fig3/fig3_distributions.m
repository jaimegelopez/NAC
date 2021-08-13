
%% Generate gaussians
clear;clc

n1 = 2;

mean = 100;
stddev = 100;
init_x = linspace(mean - 8*stddev,mean + 8*stddev,2000);
cdf = normcdf(init_x,mean,stddev);
for i = 1:n1
    [pdf{i},x,dim_means(i),dim_stds(i),areas(i)] = compute_multipdf(init_x,cdf,i);
end

growth_cvs = dim_stds./dim_means;



%% Plot

dim = 1.9;
FontSize = 8;
fig_alpha = 0.2;
ref = zeros(size(x));

gap = [0.1,0.35];
marg_h = [0.1,0.1];
newfigure(dim,dim);
Fig3ax = tight_subplot(2,2,gap,marg_h,[0.02,0.02]);

axes(Fig3ax(2))
axis off

color1 = [206 37 123]/255;
color2 = [15 104 194]/255;

gam_mean = trapz(x,(x).*pdf{1});
[~,minpoint] = min(abs(x - gam_mean));

ymax = 1.1*max(pdf{n1});
xmax = mean + 4*stddev;
xmin = mean - 4*stddev;

axes(Fig3ax(1))
hold on
line_width = 1.5;
plot(x,pdf{1},'k-','LineWidth',line_width);
fill([x,fliplr(x)],[pdf{1},fliplr(ref)],color1,'FaceAlpha',fig_alpha)
plot([gam_mean,gam_mean],[0,pdf{1}(minpoint)],'k-','LineWidth',line_width/2)
set(gca,'LineWidth',line_width)
xticks([])
yticks([])
xlim([xmin,xmax])
ylim([0,ymax])

axes(Fig3ax(3))
hold on
plot(x,pdf{1},'k-','LineWidth',line_width);
fill([x,fliplr(x)],[pdf{1},fliplr(ref)],color2,'FaceAlpha',fig_alpha)
plot([gam_mean,gam_mean],[0,pdf{1}(minpoint)],'k-','LineWidth',line_width/2)
set(gca,'LineWidth',line_width)
xticks([])
yticks([])
xlim([xmin,xmax])
ylim([0,ymax])

min_mean = trapz(x,x.*pdf{n1});
[~,minpoint2] = min(abs(x - min_mean));

axes(Fig3ax(4))
pos = get(gca,'Position');
pos(2) = 0.35;
set(gca,'Position',pos)
hold on
plot(x,pdf{n1},'k-','LineWidth',line_width);
fill([x,fliplr(x)],[pdf{n1},fliplr(ref)],[0.3,0.3,0.3],'FaceAlpha',fig_alpha)
plot([min_mean,min_mean],[0,pdf{n1}(minpoint2)],'k-','LineWidth',line_width/2)
set(gca,'LineWidth',line_width)
%plot(x,pdf{1},'k:','LineWidth',line_width,'Color',[0,0,0]);
plot([gam_mean,gam_mean],[0,pdf{1}(minpoint)],'k-','LineWidth',line_width/2)
xticks([])
yticks([])
xlim([xmin,xmax])
ylim([0,ymax])


%print(gcf,'-dpng','fig3_distribution_template.png','-r1200');


%% 

dim = 1.9;
FontSize = 8;
fig_alpha = 0.3;
ref = zeros(size(x));

gap = [0.1,0.35];
marg_h = [0.1,0.1];
newfigure(dim,0.6*dim);
Fig3ax = tight_subplot(1,2,gap,marg_h,[0.02,0.02]);


color1 = [206 37 123]/255;
color2 = [15 104 194]/255;


axes(Fig3ax(1))
hold on
line_width = 2;
plot(x,pdf{1,3},'k-','LineWidth',line_width);
fill([x,fliplr(x)],[pdf{1,3},fliplr(ref)],color1,'FaceAlpha',fig_alpha)
set(gca,'LineWidth',line_width)
xticks([])
yticks([])
xlim([0,220])
ylim([0,0.022])

axes(Fig3ax(2))
hold on
line_width = 2;
plot(x,pdf{1,3},'k-','LineWidth',line_width);
fill([x,fliplr(x)],[pdf{1,3},fliplr(ref)],color1,'FaceAlpha',fig_alpha)
set(gca,'LineWidth',line_width)
xticks([])
yticks([])
xlim([0,220])
ylim([0,0.022])

print(gcf,'-dpng','fig3a_distribution_template.png','-r1200');

