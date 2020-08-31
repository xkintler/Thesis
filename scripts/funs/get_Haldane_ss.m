function [xs,ss,ps] = get_Haldane_ss(D,params)
    mu_max = params(1);
    nu = params(2);
    Ks = params(3);
    Ki = params(4);
    Yx = params(5);
    Yp = params(6);
    s_in = params(7);
    
    ss = -(Ki*(D - mu_max + ((D^2*Ki - 4*D^2*Ks + Ki*mu_max^2 - 2*D*Ki*mu_max)/Ki)^(1/2)))/(2*D);
    xs = -(D*(ss - s_in))/(D/Yx + nu/Yp);
    ps = xs*nu/D;
end