function [pdf,x,dim_mean,dim_std,area] = compute_multipdf(init_x,cdf,n)

temp_cdf = 1 - (1 - cdf).^n;
pdf = diff(temp_cdf)./diff(init_x);
x = init_x(1:end-1) ;
dim_mean = trapz(x+x(2)/2,(x+x(2)/2).*pdf);
dim_std = sqrt(trapz((x-dim_mean).^2.*pdf));
area = trapz(x,pdf);
end

