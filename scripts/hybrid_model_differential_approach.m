clear all
close all
clc

addpath('e:\skola\semestralny projekt\1')
addpath('funs')

F = 1; %l/h
V = 3.33; %l
mu_max = 0.53; %h-1
nu = 0.5; %h-1 
Ks = 1.2; %g/l
Ki = 70; %g/l
Yx = .4;
Yp = 1;
s_in = 20; %g/l
tf = 200;
s0 = 2.9296; x0 = 4.4571; p0 = 5.9277;
is0 = 3.3120; ix0 = 4.3574; ip0 = 5.7945;

param = [mu_max nu Ks Ki Yx Yp s_in [] tf];
init_c = [x0;s0;p0];
i_init_c = [ix0;is0;ip0];

ist = []; ixt = []; ipt = [];
st = []; xt = []; pt = [];
Ts = 10; % number of measured samples through the simulation time
t = linspace(0,tf,Ts);

kmax = 3;

%Real (Monod) model optimum
alpha = .5;

a = sdpvar(1,1);
obj = monod_obj(a,alpha,param);
cons = [a >= 0; a <= mu_max];
assign(a,0.1);
ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
optimize(cons,obj,ops);
real_d = value(a);
real_j = obj;


% Nominal (Haldane) model optimum
a = sdpvar(1,1);
obj = haldane_obj(a,alpha,param);
cons = [a >= 0; a <= mu_max];
assign(a,0.1);
ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
optimize(cons,obj,ops);
nom_d = value(a);
nom_j = obj;

D = [nom_d real_d];
dD = D(2) - D(1);
param = [mu_max nu Ks Ki Yx Yp s_in D(1,1) tf];

for j = 2:1:kmax
    % hybrid model based on biomass concentration
    % Real data measurement
    [xt,~,~] = generate_Monod_data(D(1,:),param,init_c,tf,Ts);
    
    % noise addition
    rng(0)
    power = 0.01;
    noise_bio = power*(2*rand(size(xt)) - 1);
    
    m_xt = xt + noise_bio;
    
    % nominal model simulation
    [ixt,~,~] = generate_Haldane_data(D(1,:),param,i_init_c,tf,Ts);
    
    % FIR model design
    dx = m_xt - ixt;
    dx = dx(Ts+1:end);
    
    u = ones([Ts 1])*dD(1,:);
    u = u(:)';
    dx = dx';
    model_error = power;    
    rad_p = 1;
    C = zeros(rad_p);
    v = 1;
%     options = optimoptions('linprog','Algorithm','dual-simplex','Display','off');

    while v==1
       C = get_matica_vstupov(u',rad_p);
       s = size(C);
       n = zeros(1,s(2));
       v = isempty(linprog(n,[-C;C],[model_error-dx; model_error+dx]));
       if v == 1
           rad_p = rad_p + 1;
       else
           rad_p;   
       end
       if rad_p == 1000
           v = 0;
       end
    end

    P_bio = inv(C'*C)*C'*dx;
    
     % hybrid model optimum
    a = sdpvar(1,1);
    obj = df_hybrid_obj_bio(a,alpha,param,P_bio,D(1,j));
    cons = [a >= 0.1;...
            a <= mu_max];
    assign(a,0.1);
    ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
    optimize(cons,obj,ops);
    hyb_d_bio = value(a)
    hyb_bio_j = value(obj);
    
    D(1,j+1) = hyb_d_bio;
    dD(1,j) = D(1,j+1) - D(1,j);
end