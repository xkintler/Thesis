clc 
clear all
close all

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
tf = 50;
s0 = 1.5688; x0 = 4.4252; p0 = 7.3680; 
is0 = 1.6177; ix0 = 4.4135; ip0 = 7.3485; 

param = [mu_max nu Ks Ki Yx Yp s_in [] tf];

Ts = 20; % number of measured samples through the simulation time
t = linspace(0,tf,Ts);

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

D = nom_d;

[x0,s0,p0] = get_Monod_ss(D,param);
[ix0,is0,ip0] = get_Haldane_ss(D,param);

init_c = [x0;s0;p0];
i_init_c = [ix0;is0;ip0];

kmax = 8; % number of step changes
c = .8; % weigh constant

for j = 1:1:kmax
    % hybrid model based on biomass concentration
    % Real data measurement
    [xt,st,~] = generate_Monod_data(D,param,init_c,tf,Ts);
    
    % noise addition
    rng(0)
    power = 0.01;
    noise_bio = power*(2*rand(size(xt)) - 1);
    
    rng(0)
    power = 0.025;
    noise_sub = power*(2*rand(size(xt)) - 1);
    
    m_xt = xt + noise_bio;
    m_st = st + noise_sub;
    
    % modifier adaptation scheme approach
    if j == 1
        lam = -0.07;
    else
        J0 = D(j).*(1 - alpha*m_xt(end));
        J1 = D(j-1).*(1 - alpha*m_xt(end-Ts));
        real_grad_app = (J0 - J1)/(D(j) - D(j-1));            
        temp_lam = real_grad_app - grad_haldane_obj(D(j),alpha,param);
        lam = c*lam + (1 - c)*temp_lam;
    end
    
    a = sdpvar(1,1);
    obj = modifier_adapt(a,alpha,param,lam);
    cons = [a >= 0; a <= mu_max];
    assign(a,0.1);
    ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
    optimize(cons,obj,ops);
    mas_d = value(a)
    
    D = [D mas_d];
end

D = D(1:kmax);

% nominal model simulation
[~,ist,~] = generate_Haldane_data(D,param,i_init_c,tf,Ts);

% FIR model design
ds = m_st - ist;
u = ones([Ts 1])*D;
u = u(:);
ds = ds';

u = u - nom_d;
dsss = s0 - is0;
ds = ds - dsss;

model_error = 2*power;    

[r_sub, P_sub] = gpe_fir_min_order(u,ds,model_error,1000,0,0);
[P_sub_min, P_sub_max] = gpe_fir_min_max_bound(r_sub,u,ds,model_error,0);

C = get_matica_vstupov(u,r_sub);
tt = linspace(0,kmax*tf,length(ds));

figure
hold on
stairs(tt,ds+dsss,'-r','MarkerFaceColor','r','LineWidth',1.5)
stairs(tt,C*P_sub+dsss,'-b','LineWidth',2)
xlabel('Čas [hod]')
ylabel('Rozdiel koncentrií substrátu \Delta_{s} [gL^{-1}]')
legend('Namerané údaje','Výstup FIR modelu','Location','Best')
set(gca,'FontSize',15)
box on
hold off

% Comparison figures
D = [nom_d 0.385 real_d]; % step change
% Plant data
[~,st,~] = generate_Monod_data(D,param,init_c,tf,Ts);
    % noise addition
    rng(1)
    power = 0.025;
    noise_sub = power*(2*rand(size(st)) - 1);
    m_st = st + noise_sub;
% Nominal model data
[~,ist,~] = generate_Haldane_data(D,param,i_init_c,tf,Ts);
% Hybrid model data
u = ones([Ts 1])*D;
u = u(:)-nom_d;
% Ak sa ma zacat z roznych ustalenych stavov, je nutne zmenit maticu
% vstupov tak, aby uz v case u(t-0) bol zaznamenany nejaky vstupny signal.
C = get_matica_vstupov(u,r_sub);

hst = ist + (C*P_sub)' + dsss;
hst_min = ist + (C*P_sub_min)' + dsss;
hst_max = ist + (C*P_sub_max)' + dsss;

figure
nd = length(D);
time = linspace(0,nd*tf,nd*Ts);

hold on
stairs(time,m_st,'-c','LineWidth',2)
stairs(time,ist,'-m','LineWidth',2)
stairs(time,hst,'-g','LineWidth',2)
stairs(time,hst_min,'-r','LineWidth',2)
stairs(time,hst_max,'-b','LineWidth',2)
xlabel('Čas [hod]')
ylabel('Koncentrácia substrátu [gL^{-1}]')
legend('Zariadenie','Nominálny model','Hybridný model','Location','Best')
set(gca,'FontSize',15)
box on


