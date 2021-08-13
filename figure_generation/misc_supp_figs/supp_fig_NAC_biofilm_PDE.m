clear;clc

D = 500; %microns^2/s
rho = 0.99; %cells/micron^3
R = 0.55; %microns
a = 5; %microns
mu = rho*4*pi*R*D; %1/s
invnu = sqrt(mu/D);
xi = 2.5e6; %molecules/um^3 s

nondim_a = a*invnu;


Jout_denom = mu.*(1./sqrt(mu*D)).*((besselk(0,nondim_a)./besselk(1,nondim_a))...
    + (besseli(0,nondim_a)./besseli(1,nondim_a)));
Jout = xi./Jout_denom;

n1 = 1e4;
r = linspace(0,25*invnu,n1);

m = biofilm_PDE_soln(r,mu,D,xi,Jout,nondim_a);

gap = [0.15,0.15];
marg_h = [0.2,0.1];
newfigure(2*0.8*3.42/2, (1/2)*3.42/3*(3.5/2));
Figax = tight_subplot(1,2,gap,marg_h,[0.1,0.1]);
axes(Figax(1))
plot([a,a],[0,2*max(m)],'k-.')
hold on
plot(r*(1/invnu),m,'r-')
ylim([0,800])
xlim([0,25])
xlabel('$\hat{r}$ ($\mu$m)','Interpreter','latex')
ylabel('$m$ (molecules/$\mu$m$^3$)','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',6)

a = 20; %microns
nondim_a = a*invnu;
Jout_denom = mu.*(1./sqrt(mu*D)).*((besselk(0,nondim_a)./besselk(1,nondim_a))...
    + (besseli(0,nondim_a)./besseli(1,nondim_a)));
Jout = xi./Jout_denom;
r = linspace(0,25*invnu,n1);
m = biofilm_PDE_soln(r,mu,D,xi,Jout,nondim_a);
axes(Figax(2))
plot([a,a],[0,2*max(m)],'k-.')
hold on
plot(r*(1/invnu),m,'r-')
ylim([0,800])
yticks([])
xlim([0,25])
xlabel('$\hat{r}$ ($\mu$m)','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',6)
print(gcf,'-dpng','supp_fig_NAC_biofilm_PDE.png','-r1200');

