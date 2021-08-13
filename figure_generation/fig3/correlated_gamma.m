function [X,growth] = correlated_gamma(n,m,k,theta,rho)

Rho = ones(n,n)*rho;
Rho(1:n+1:end) = 1;

Z = mvnrnd(zeros(1,n), Rho, m);
U = normcdf(Z,0,1);

X = zeros(m,n);
for i = 1:n
    X(:,i) = gaminv(U(:,i),k,theta);
end

growth = min(X,[],2);

end

