function [xt,st,pt] = generate_Monod_data(D,params,init_cond,tf,Ts)
% params = [mu_max nu Ks Ki Yx Yp s_in ];
% init_c = [x0; s0; p0];
% tf -> final time
% Ts -> number of samples (sampling time)

    x0 = init_cond(1);
    s0 = init_cond(2);
    p0 = init_cond(3);
    
    mu_max = params(1);
    nu = params(2);
    Ks = params(3);
    Ki = params(4);
    Yx = params(5);
    Yp = params(6);
    s_in = params(7);
    
    t = linspace(0,tf,Ts);
    
    st = [];
    xt = [];
    pt = [];
    
    for i = 1:1:length(D)
        param = [mu_max nu Ks Ki Yx Yp s_in D(i) tf];
        S = ode45(@(t,y) myOdes(y,param),[0 tf],[s0;x0;p0]);
        st = [st deval(S,t,1)];
        xt = [xt deval(S,t,2)];
        pt = [pt deval(S,t,3)];

        s0 = st(end); x0 = xt(end); p0 = pt(end);
    end
end