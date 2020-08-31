clear all
close all
clc

addpath('funs')

F = 1; %l/h
V = 3.33; %l
D = F/V;
mu_max = 0.53; %h-1
nu = 0.5; %h-1 
Ks = 1.2; %g/l
Ki = 70; %g/l
Yx = .4;
Yp = 1;
s_in = 20; %g/l

alpha = .5;
Ts = 30; % number of measured samples through the simulation time
tf = 50;

s0 = 2.9291; x0 = 4.4573; p0 = 5.9277;
is0 = 3.3115; ix0 = 4.3574; ip0 = 5.7949;

init_c = [x0;s0;p0];
param = [mu_max nu Ks Ki Yx Yp s_in 0 tf];

ist = []; ixt = []; ipt = [];
st = []; xt = []; pt = [];

%Real (Monod) model optimum
a = sdpvar(1,1);
obj = monod_obj(a,alpha,param);
cons = [a >= 0; a <= mu_max];
assign(a,0.1);
ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
optimize(cons,obj,ops);
real_d = value(a);
real_j = value(obj);

% Nominal (Haldane) model optimum
a = sdpvar(1,1);
obj = haldane_obj(a,alpha,param);
cons = [a >= 0; a <= mu_max];
assign(a,0.1);
ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
optimize(cons,obj,ops);
nom_d = value(a);
nom_j = obj;

Ts = [10];
for i = 1:1:length(Ts)
    D = nom_d;
    
    kmax = 10;
    est_mu = .5;
    est_Ks = 1;
    est_Ki = 50;
    
    for j = 1:1:kmax
        % Real data measurement
        [xt,st,~] = generate_Monod_data(D,param,[x0;s0;p0],tf,Ts(i));
        % noise addition
        rng(0)
        power = 0.05;
        noise_bio = power*(2*rand(size(xt)) - 1);
        rng(1)
        power = 0.025;
        noise_sub = power*(2*rand(size(xt)) - 1);
        
        m_xt = xt + noise_bio;
        m_st = st + noise_sub;
        
        % parameter estimation approach
        p1 = {mu_max nu Ks Ki Yx Yp s_in D tf Ts(i)};
        opts = optimset('MaxFunEvals',1e5,'MaxIter',1e3);
        [X,FVAL,EXITFLAG,OUTPUT] = fminsearch(@(a) getValueOfCostFun(a,m_st,init_c,p1),[est_mu;est_Ks;est_Ki],opts);
        est_mu = X(1);
        est_Ks = X(2);
        est_Ki = X(3);
        
        est_param = [est_mu nu est_Ks est_Ki Yx Yp s_in D(j) tf];
        [ex0,es0,ep0] = get_Haldane_ss(nom_d,est_param);
        
        [~,est,~] = generate_Haldane_data(D,est_param,[ex0;es0;ep0],tf,Ts(i));
        
        % Optimal solution based on parameter estimation
        a = sdpvar(1,1);
        obj = haldane_obj(a,alpha,est_param);
        cons = [a >= 0; a <= mu_max; a <= est_mu];
        assign(a,0.3);
        ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
        optimize(cons,obj,ops);
        pe_d = value(a);
        pe_j = value(obj);
        fprintf(['\nParameter estimation approach (Haldane) optimal D : ', num2str(pe_d), '\n'])
        
        param_warehouse(:,j) = [est_mu; est_Ks; est_Ki];
        D = [D pe_d];
        
        figure(5)
        hold on
        plot(m_st)
        plot(est)
    end
    cell_data{i} = [D(1:end-1); param_warehouse];
end
%% Figures
close all
clr = {'b','r','g'};

for i = 1:1:length(Ts)
    D = cell_data{i}(1,:);
    param_warehouse = cell_data{i}(2:end,:);

    figure(1)
    iter = 1:1:j;
    hold on
    plot([1 j],[mu_max mu_max],':k','LineWidth',1.5)
    stairs(iter,param_warehouse(1,:),clr{i},'LineWidth',2)
    xlabel('Iterácia')
    ylabel('\mu_{max} [hod^{-1}]')
    xlim([1 kmax])
    set(gca,'FontSize',18)
    box on
    hold off
    
    figure(2)
    hold on
    plot([1 j],[Ks Ks],':k','LineWidth',1.5)
    stairs(iter,param_warehouse(2,:),clr{i},'LineWidth',2)
    xlabel('Iterácia')
    ylabel('K_M [gL^{-1}]')
    xlim([1 kmax])
    set(gca,'FontSize',18)
    box on
    hold off
    
    figure(3)
    hold on
    plot([1 j],[Ki Ki],':k','LineWidth',1.5)
    stairs(iter,param_warehouse(3,:),clr{i},'LineWidth',2)
    xlabel('Iterácia')
    ylabel('K_I [gL^{-1}]')
    xlim([1 kmax])
    set(gca,'FontSize',18)
    box on
    hold off
    
    figure(4)
    hold on
    fd = monod_obj(D,alpha,param);
    plot([iter(1) iter(end)],[monod_obj(real_d,alpha,param) monod_obj(real_d,alpha,param)],':k','LineWidth',1.5)
    stairs(iter,fd,clr{i},'LineWidth',2)
    xlabel('Iterácia')
    ylabel('Hodnota účelovej funkcie J_{Monod}')
    axis([1 kmax -0.475 -0.445])
    set(gca,'FontSize',12)
    box on
    hold off
    
    if i == 1 || i == 3     
       est_mu = param_warehouse(1,end);
       est_Ks = param_warehouse(2,end);
       est_Ki = param_warehouse(3,end);
       est_param =  [est_mu nu est_Ks est_Ki Yx Yp s_in D(1) tf];
       
       a = sdpvar(1,1);
       obj = haldane_obj(a,alpha,est_param);
       cons = [a >= 0; a <= mu_max; a <= est_mu];
       assign(a,0.3);
       ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
       opt=optimize(cons,obj,ops);
       pe_d = value(a);
       pe_j = value(obj);
       
       dd = 0:0.00001:0.53;
       figure(7)
       hold on
       plot(dd,monod_obj(dd,alpha,param),'-b')
       plot(real_d,real_j,'ok','MarkerFaceColor','k')
       if i == 1
           plot(dd,haldane_obj(dd,alpha,est_param),'-r')
           plot(pe_d,pe_j,'ok','MarkerFaceColor','k')
       else
           plot(dd,haldane_obj(dd,alpha,est_param),'-g')
           plot(pe_d,pe_j,'ok','MarkerFaceColor','k')
       end
       
    end
end

figure(6)
hold on
for i = 1:1:length(Ts)
    D = cell_data{i}(1,:);
    fd = monod_obj(D,alpha,param);
    stairs(iter,fd,clr{i},'LineWidth',2)
end
plot([iter(1) iter(end)],[monod_obj(real_d,alpha,param) monod_obj(real_d,alpha,param)],':k','LineWidth',1.5)
xlabel('Iterácia')
ylabel('Hodnota účelovej funkcie J_{Monod}')
legend('Ts = 5 hod', 'Ts = 2.5 hod', 'Ts = 1.67 hod','Location','Best')
% set(gca,'YScale','log')
box on
hold off

%% derivation approximation
clear all

addpath('funs')

F = 1; %l/h
V = 3.33; %l
D = F/V;
mu_max = 0.53; %h-1
nu = 0.5; %h-1 
Ks = 1.2; %g/l
Ki = 70; %g/l
Yx = .4;
Yp = 1;
s_in = 20; %g/l

alpha = .5;
Ts = 10; % number of measured samples through the simulation time
tf = 50;

s0 = 2.9291; x0 = 4.4573; p0 = 5.9277;
is0 = 3.3115; ix0 = 4.3574; ip0 = 5.7949;
es0 = 2.9178; ex0 = 4.4602; ep0 = 5.9316;

init_c = [ix0;is0;ip0];
param = [mu_max nu Ks Ki Yx Yp s_in 0 tf];

ist = []; ixt = []; ipt = [];
st = []; xt = []; pt = [];

% Nominal (Haldane) model optimum
a = sdpvar(1,1);
obj = haldane_obj(a,alpha,param);
cons = [a >= 0; a <= mu_max];
assign(a,0.1);
ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
optimize(cons,obj,ops);
nom_d = value(a);
nom_j = obj;

D = [nom_d 0.38 0.385 0.39 0.395];

% Real data measurement
[xt,st,~] = generate_Monod_data(D,param,[x0;s0;p0],tf,Ts);
% noise addition
rng(0)
power = 0.05;
noise_bio = power*(2*rand(size(xt)) - 1);
rng(1)
power = 0.025;
noise_sub = power*(2*rand(size(xt)) - 1);

m_xt = xt + noise_bio;
m_st = st + noise_sub;

tt = linspace(0,length(D)*tf,length(m_xt));

% derivation approx
dx = diff(m_xt);
dt = diff(tt);

der = dx./dt;
der = [der der(end)];
Ds = D.*ones([Ts length(D)]);
Ds = Ds(:)';

obj = @(a) sum((der - ((a(1)*m_st./(a(2) + m_st + (m_st.^2)/a(3))) - Ds).*m_xt).^2);
opts = optimset('MaxFunEvals',1e5,'MaxIter',1e5);
sol = fminsearch(obj,[0.5;1;50],opts);
est_mu = sol(1);
est_Ks = sol(2);
est_Ki = sol(3);

est_param = [est_mu nu est_Ks est_Ki Yx Yp s_in 0 tf];

[ext,est,~] = generate_Haldane_data(D,est_param,[ex0;es0;ep0],tf,Ts);

figure
hold on
stairs(tt,m_st,':k','MarkerFaceColor','r','LineWidth',1.5)
stairs(tt,st,'-r','MarkerFaceColor','r','LineWidth',2)
stairs(tt,est,'-b','LineWidth',2)
xlabel('Čas [hod]')
ylabel('Koncentrácia substrátu [gL^{-1}]')
legend('Zariadenie (Monod)','Bez šumu (Monod)','Odhad (Haldane)','Location','Best')
set(gca,'FontSize',15)
box on
hold off

figure
hold on
stairs(tt,m_xt,':k','MarkerFaceColor','b','LineWidth',1.5)
stairs(tt,xt,'-r','MarkerFaceColor','r','LineWidth',2)
stairs(tt,ext,'-b','LineWidth',2)
xlabel('Čas [hod]')
ylabel('Koncentrácia biomasy [gL^{-1}]')
legend('Zariadenie (Monod)','Bez šumu (Monod)','Odhad (Haldane)','Location','Best')
box on
hold off

%% Estimation of diff eq parameters
clear all

addpath('funs')

F = 1; %l/h
V = 3.33; %l
D = F/V;
mu_max = 0.53; %h-1
nu = 0.5; %h-1 
Ks = 1.2; %g/l
Ki = 70; %g/l
Yx = .4;
Yp = 1;
s_in = 20; %g/l

alpha = .5;
Ts = 10; % number of measured samples through the simulation time
tf = 50;

s0 = 2.9291; x0 = 4.4573; p0 = 5.9277;
is0 = 3.3115; ix0 = 4.3574; ip0 = 5.7949;
es0 = 2.9263; ex0 = 4.4580; ep0 = 5.9286;
est_mu = 0.5; est_Ks = 1; est_Ki = 50;

init_c = [ix0;is0;ip0];
param = [mu_max nu Ks Ki Yx Yp s_in 0 tf];

ist = []; ixt = []; ipt = [];
st = []; xt = []; pt = [];

% Nominal (Haldane) model optimum
a = sdpvar(1,1);
obj = haldane_obj(a,alpha,param);
cons = [a >= 0; a <= mu_max];
assign(a,0.1);
ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
optimize(cons,obj,ops);
nom_d = value(a);
nom_j = obj;

D = [nom_d 0.38 0.385 0.39 0.395];

% Real data measurement
[~,st,~] = generate_Monod_data(D,param,[x0;s0;p0],tf,Ts);
% noise addition

rng(1)
power = 0.025;
noise_sub = power*(2*rand(size(st)) - 1);

m_st = st + noise_sub;

tt = linspace(0,length(D)*tf,length(m_st));

% parameter estimation approach
p1 = {mu_max nu Ks Ki Yx Yp s_in D tf Ts};
opts = optimset('MaxFunEvals',1e5);
[X,FVAL,EXITFLAG,OUTPUT] = fminsearch(@(a) getValueOfCostFun(a,m_st,init_c,p1),[est_mu;est_Ks;est_Ki],opts);
est_mu = X(1);
est_Ks = X(2);
est_Ki = X(3);

est_param = [est_mu nu est_Ks est_Ki Yx Yp s_in D(1) tf];

[ext,est,ept] = generate_Haldane_data(D,est_param,[ex0;es0;ep0],tf,Ts);

figure
hold on
stairs(tt,m_st,':k','MarkerFaceColor','r','LineWidth',1.5)
stairs(tt,st,'-r','MarkerFaceColor','r','LineWidth',2)
stairs(tt,est,'-b','LineWidth',2)
xlabel('Čas [hod]')
ylabel('Koncentrácia substrátu [gL^{-1}]')
legend('Zariadenie (Monod)','Bez šumu (Monod)','Odhad (Haldane)','Location','Best')
set(gca,'FontSize',15)
box on
hold off

%% situacia ked Haldane predstavuje zariadenie
addpath('funs')

F = 1; %l/h
V = 3.33; %l
D = F/V;
mu_max = 0.53; %h-1
nu = 0.5; %h-1 
Ks = 1.2; %g/l
Ki = 70; %g/l
Yx = .4;
Yp = 1;
s_in = 20; %g/l

alpha = .5;
Ts = 10; % number of measured samples through the simulation time
tf = 50;

s0 = 2.9291; x0 = 4.4573; p0 = 5.9277;
is0 = 3.3115; ix0 = 4.3574; ip0 = 5.7949;

init_c = [x0;s0;p0];
param = [mu_max nu Ks Ki Yx Yp s_in 0 tf];

ist = []; ixt = []; ipt = [];
st = []; xt = []; pt = [];

%Real (Monod) model optimum
a = sdpvar(1,1);
obj = monod_obj(a,alpha,param);
cons = [a >= 0; a <= mu_max];
assign(a,0.1);
ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
optimize(cons,obj,ops);
real_d = value(a);
real_j = value(obj);

% Nominal (Haldane) model optimum
a = sdpvar(1,1);
obj = haldane_obj(a,alpha,param);
cons = [a >= 0; a <= mu_max];
assign(a,0.1);
ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
optimize(cons,obj,ops);
nom_d = value(a);
nom_j = obj;

D = real_d;

[ix0,is0,ip0] = get_Haldane_ss(D,param);
    
kmax = 10;
est_mu = .5;
est_Ks = 1;
est_Ki = 50;

for j = 1:1:kmax
    
    % Real data measurement
    [xt,st,~] = generate_Haldane_data(D,param,[ix0;is0;ip0],tf,Ts);
    % noise addition
    rng(0)
    power = 0.05;
    noise_bio = power*(2*rand(size(xt)) - 1);
    rng(1)
    power = 0.025;
    noise_sub = power*(2*rand(size(xt)) - 1);
    
    m_xt = xt + noise_bio;
    m_st = st + noise_sub;
    
    % parameter estimation approach
    p1 = {mu_max nu Ks Ki Yx Yp s_in D tf Ts};
    opts = optimset('MaxFunEvals',1e5,'MaxIter',1e3);
    [X,FVAL,EXITFLAG,OUTPUT] = fminsearch(@(a) getValueOfCostFun_H(a,m_st,[ix0;is0;ip0],p1),[est_mu;est_Ks;est_Ki],opts);
    est_mu = X(1);
    est_Ks = X(2);
    est_Ki = X(3);
    
    est_param = [est_mu nu est_Ks est_Ki Yx Yp s_in D(j) tf];
    [ex0,es0,ep0] = get_Monod_ss(real_d,est_param);
    
    [~,est,~] = generate_Monod_data(D,est_param,[ex0;es0;ep0],tf,Ts);
    
    % Optimal solution based on parameter estimation
    a = sdpvar(1,1);
    obj = monod_obj(a,alpha,est_param);
    cons = [a >= 0; a <= mu_max; a <= est_mu];
    assign(a,0.3);
    ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
    optimize(cons,obj,ops);
    pe_d = value(a);
    pe_j = value(obj);
    fprintf(['\nParameter estimation approach (Haldane) optimal D : ', num2str(pe_d), '\n'])
    
    param_warehouse(:,j) = [est_mu; est_Ks; est_Ki];
    D = [D pe_d];
    
    figure(5)
    hold on
    plot(m_st)
    plot(est)
end