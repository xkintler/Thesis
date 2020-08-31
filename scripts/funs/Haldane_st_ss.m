function ist = Haldane_st_ss(D,param)
    mu_max = param(1); %h-1
    Ks = param(3); %g/l
    Ki = param(4); %g/l
    
    ist = -(Ki*(D - mu_max + ((D^2*Ki - 4*D^2*Ks + Ki*mu_max^2 - 2*D*Ki*mu_max)/Ki)^(1/2)))/(2*D);    
end