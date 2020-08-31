function J = getValueOfCostFun_H(a,m_sub,init,param)
    n = param{10};
    s0 = init(2);
    x0 = init(1);
    p0 = init(3);
    tf = param{9};
    D = param{8}; %h-1
    
    xt = [];
    st = [];
    pt = []; 
    t = linspace(0,tf,n);
    
    for k = 1:1:length(D)
        S = ode45(@(t,y) myInnerOdes(y,a,[1 param{2} 1 1 param{5} param{6} param{7} D(k) tf]),[0 tf],[s0;x0;p0]);
        st = [st deval(S,t,1)];
        xt = [xt deval(S,t,2)];
        pt = [pt deval(S,t,3)];
        x0 = xt(end);
        s0 = st(end);
        p0 = pt(end);
    end
    J = sum((m_sub - st).^2);
    
    function f = myInnerOdes(y,a,param)    
    
    nu = param(2); %h-1 
    Yx = param(5);
    Yp = param(6);
    s_in = param(7); %g/l
    d = param(8); %h-1
    
    f = [- ((a(1)*y(1)/(a(2) + y(1)))*y(2)/Yx) - (nu*y(2)/Yp) + d*(s_in - y(1));...
           (a(1)*y(1)/(a(2) + y(1)))*y(2) - d*y(2);...
            nu*y(2) - d*y(3)];
    end
end