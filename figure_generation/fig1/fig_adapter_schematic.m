%This script generates a schematic detailing the dynamics of the adapter
close all

fig_x = 3.42;
fig_y = 1.5;
newfigure(fig_x,fig_y);
%Set up master axis
masterax = axes('Position',[0 0 1 1]);
xlim([0 1])
ylim([0 1])
axis off
hold on

LineWidth = 1;

%Make central diff equation
FontSize = 11; 
fig_x = 1;
fig_y = 1;
xlim([0,fig_x])
ylim([0,fig_y])
eq_frac_x = 0.35;
eq_frac_y = 0.5;
text(eq_frac_x*fig_x,0.98*eq_frac_y*fig_y,{'Which metabolite','is low?'},...
    'Interpreter','latex','HorizontalAlignment','center',...
    'VerticalAlignment','middle','FontSize',FontSize)

%Make central rectangle
rec_x = 1.3/3.42;
rec_y = 0.48/2;
rec_position = [eq_frac_x*fig_x - 0.5*rec_x,eq_frac_y*fig_y-0.5*rec_y,rec_x,rec_y];
rectangle('Position',rec_position,'Curvature',0,'LineWidth',LineWidth)

%Make c_1 and c_2 inputs
c_x = 0.1/3.42;
c_offset = 0.5/2;
c_pad = 0.03;
text(c_x,0.5*fig_y + c_offset,'$m_1$','FontSize',FontSize,...
    'Interpreter','latex','HorizontalAlignment','center',...
    'VerticalAlignment','middle')
text(c_x,0.5*fig_y - c_offset,'$m_2$','FontSize',FontSize,...
    'Interpreter','latex','HorizontalAlignment','center',...
    'VerticalAlignment','middle')
rec_left = eq_frac_x*fig_x - eq_frac_y*rec_x;
rec_center = 0.5*fig_y;

annotation('arrow',[c_x+c_pad, rec_left],[0.5*fig_y + c_offset-c_pad,rec_center+c_pad],...
    'Units','normalized','HeadLength',5,'HeadWidth',5,'LineWidth',LineWidth);
annotation('arrow',[c_x+c_pad, rec_left],[0.5*fig_y - c_offset+c_pad,rec_center-c_pad],...
    'Units','normalized','HeadLength',5,'HeadWidth',5,'LineWidth',LineWidth);


%Make three outcomee
ellip_x = 0.71;
box_strings = {'Make enzyme 2','No change','Make enzyme 1'};
rec_right = eq_frac_x*fig_x + eq_frac_y*rec_x;
new_rec_x = 0.75*rec_x;
new_rec_y = 0.8*rec_y;
ellip_y_vec = [0.31,0.5,0.69]-0.5*new_rec_y;
theta = zeros([1,3]);
eqn_mid = mean([rec_right,ellip_x]);

center_offsets = [-0.03,0,0.03];
for i = 1:3
    if sum(i == [1,3]) == 1
    rec_position = [ellip_x,ellip_y_vec(i),new_rec_x,new_rec_y];
    rectangle('Position',rec_position,'Curvature',1,'LineWidth',LineWidth)
    text(ellip_x+0.5*+new_rec_x,ellip_y_vec(i)+0.5*new_rec_y,box_strings(i),...
        'HorizontalAlignment','center',...
        'VerticalAlignment','middle',...
        'FontName','serif','FontSize',FontSize-2,'Interpreter','latex');
    annotation('arrow',[rec_right,ellip_x],...
        [rec_center+center_offsets(i),ellip_y_vec(i)+0.5*new_rec_y],...
    'Units','normalized','HeadLength',5,'HeadWidth',5,'LineWidth',LineWidth);
    else
    end
end


mid_offset = -0.00;
text_offset = 0.0;
t_offset = 0.06;
text(eqn_mid-mid_offset,ellip_y_vec(3)+0.5*new_rec_y-text_offset+t_offset,'$m_1 < m_2$','HorizontalAlignment','center',...
    'VerticalAlignment','middle','FontSize',FontSize,'Interpreter','latex')

b_offset = -0.04;
text(eqn_mid-mid_offset,ellip_y_vec(1)+0.5*new_rec_y-text_offset+b_offset,'$m_2 < m_1$','HorizontalAlignment','center',...
    'VerticalAlignment','middle','FontSize',FontSize,'Interpreter','latex')

set(gca,'FontName','sans-serif')

print(gcf,'-dpng','regulation_schematic.png','-r600')
