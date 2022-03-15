function m = biofilm_PDE_soln(r,mu,D,xi,Jout,nondim_a)

%The script computes the biofilm PDF soln

m = (r>nondim_a).*(Jout./sqrt(mu*D)).*(besselk(0,r)./besselk(1,nondim_a)) ...
    + (r<=nondim_a).*(xi/mu - (Jout./(sqrt(mu.*D))).*(besseli(0,r)./besseli(1,nondim_a)));

end

