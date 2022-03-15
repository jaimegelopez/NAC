%This script generates the data for figure 3

clear;clc

%% Generate distributions for varying CV and n

init_x = linspace(0,750,20000);
prot_mean = (5.2*6.4);
n1 = 20;
cvs = linspace(0.05,1,9);
cvs(4) = 0.4;

for j = 1:length(cvs)
    k = (1/cvs(j))^2;
    theta = prot_mean/k;
    cdf = gamcdf(init_x,k,theta);
    for i = 1:n1
        [pdf{i,j},x,dim_means(i,j),dim_stds(i,j),areas(i,j)] = compute_multipdf(init_x,cdf,i);
    end
end

growth_cvs = dim_stds./dim_means;


%% Generate correlated gamma distributions

n_vec = 1:20;
rho_vec = linspace(0,0.9,10);
m = 50000;
k = (1/0.4)^2;
theta = prot_mean/k;

results = zeros(length(n_vec),length(rho_vec));

for i = 1:length(n_vec)
    for j = 1:length(rho_vec)
        [X,growth] = correlated_gamma(n_vec(i),m,k,theta,rho_vec(j));
        results(i,j) = mean(growth);
    end
end


save('fig3_data.mat')
