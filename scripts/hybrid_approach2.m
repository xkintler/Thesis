%% iteracny pristup
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

tf = 50;
Ts = 10; % number of measured samples through the simulation time

param = [mu_max nu Ks Ki Yx Yp s_in [] tf];

kmax = 20;

%Real (Monod) model optimum
alpha = .5;

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
nom_j = value(obj);

D = [nom_d; nom_d];

[x0,s0,p0] = get_Monod_ss(D(1),param);
[ix0,is0,ip0] = get_Haldane_ss(D(1),param);
init_c = [x0;s0;p0];
i_init_c = [ix0;is0;ip0];

param = [mu_max nu Ks Ki Yx Yp s_in D(1,1) tf];

for j = 1:1:kmax
    % hybrid model based on substrate 
    % Real data measurement
    [~,st,~] = generate_Monod_data(D(1,:),param,init_c,tf,Ts);
        
    % noise addition
    rng(1)
    power = 0.025;
    noise_sub = power*(2*rand(size(st)) - 1);

    m_st = st + noise_sub;
    
    % nominal model simulation
    [~,ist,~] = generate_Haldane_data(D(1,:),param,i_init_c,tf,Ts);
        
    %FIR model design
    ds = m_st - ist;
    ds = ds';
    u_sub = ones([Ts 1])*D(1,:);
    u_sub = u_sub(:)';
    
    u_sub = u_sub - nom_d;
    ds0 = s0 - is0;
    ds = ds - ds0;
    
    model_error = 2*power;
    
    [r_sub(j),P_sub] = gpe_fir_min_order(u_sub,ds,model_error,1000,1,0);
        
    % hybrid model optimum 
    a = sdpvar(1,1);
    obj = hybrid_obj_sub_corr(a,nom_d,ds0,alpha,param,P_sub);
    cons = [a >= 0.1;... 
            a <= mu_max];
    assign(a,0.1);
    ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
    optimize(cons,obj,ops);
    hyb_d_sub = value(a);
    hyb_sub_j = value(obj);
    
    fprintf(['\nHybrid model (Substrate) optimal D : ', num2str(hyb_d_sub), '\n'])
    
    % hybrid model based on biomass concentration
    % Real data measurement
    [xt,~,~] = generate_Monod_data(D(2,:),param,init_c,tf,Ts);
    
    % noise addition
    rng(0)
    power = 0.05;
    noise_bio = power*(2*rand(size(xt)) - 1);
    
    m_xt = xt + noise_bio;
    
    % nominal model simulation
    [ixt,~,~] = generate_Haldane_data(D(2,:),param,i_init_c,tf,Ts);
    
    % FIR model design
    dx = m_xt - ixt;
    u_bio = ones([Ts 1])*D(2,:);
    u_bio = u_bio(:)';
    dx = dx';
    
    u_bio = u_bio - nom_d;
    dx0 = x0 - ix0;
    dx = dx - dx0;
    
    model_error = 2*power;    
    
    [r_bio(j),P_bio] = gpe_fir_min_order(u_bio,dx,model_error,1000,0,0);
    
    % hybrid model optimum
    a = sdpvar(1,1);
    obj = hybrid_obj_bio_corr(a,nom_d,dx0,alpha,param,P_bio);
    cons = [a >= 0.1;...
            a <= mu_max];
    assign(a,0.1);
    ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
    optimize(cons,obj,ops);
    hyb_d_bio = value(a);
    hyb_bio_j = value(obj);
    
    fprintf(['\nHybrid model (Substrate) optimal D : ', num2str(hyb_d_bio), '\n'])
    
    D(:,j+1) = [hyb_d_sub; hyb_d_bio];
    P_us(:,j) = [sum(P_sub)*(D(1,j)-nom_d) + ds0; sum(P_bio)*(D(2,j)-nom_d) + dx0];
    optimal_d(:,j) = [hyb_d_sub; hyb_d_bio];
end

% steady state of P_sub for optimal D*
P_bio_limit = Monod_xt_ss(real_d,param) - Haldane_xt_ss(real_d,param);
P_sub_limit = Monod_st_ss(real_d,param) - Haldane_st_ss(real_d,param);

%% Figures 
close all

figure(1)
C = get_matica_vstupov(u_sub,r_sub(j));
fir_out = C*P_sub + ds0;
tt = linspace(0,kmax*tf,length(st));

hold on
stairs(tt,ds + ds0,'-r','LineWidth',1.5)
stairs(tt,st-ist,'-g','LineWidth',1.5)
stairs(tt,fir_out,'-b','LineWidth',2)
xlabel('Cas [hod]')
ylabel('Rozdiel koncentracie substratu \Delta_{s} [gL^{-1}]')
box on
hold off
close(1)

figure(2)
C = get_matica_vstupov(u_bio,r_bio(j));
fir_out = C*P_bio + dx0;

hold on 
stairs(tt,dx + dx0,'-b','LineWidth',1.5)
stairs(tt,xt-ixt,'-g','LineWidth',1.5)
stairs(tt,fir_out,'-r','LineWidth',2)
xlabel('Cas [hod]')
ylabel('Rozdiel koncentracie biomasy \Delta_{x} [gL^{-1}]')
box on
hold off
close(2)

figure(3)
iter = 1:1:j;
hold on
plot(iter,P_sub_limit*ones(size(iter)),':k','LineWidth',1.5)
stairs(iter,P_us(1,:),'-r','LineWidth',2)
xlabel('Iterácia')
ylabel('Ustálený stav FIR modelu')
set(gca,'FontSize',12)
xlim([iter(1) iter(end)])
box on
hold off

figure(4)
fs = monod_obj(optimal_d(1,:),alpha,param);

hold on
stairs(iter,fs,'-r','LineWidth',2)
plot(iter,monod_obj(real_d,alpha,param)*ones(size(iter)),':k','LineWidth',1.5)
xlabel('Iterácia')
ylabel('Hodnota účelovej funkcie J_{Monod}')
set(gca,'FontSize',12)
xlim([iter(1) iter(end)])
ylim([-0.471 -0.462])
box on
hold off

figure(5)
stairs(iter,r_sub,'-r','LineWidth',2)
xlabel('Iterácia')
ylabel('Rád FIR modelu')
set(gca,'FontSize',12)
xlim([iter(1) iter(end)])
ylim([0 2.5])
box on

figure(6)
fx = monod_obj(optimal_d(2,:),alpha,param);

hold on
stairs(iter,fx,'-b','LineWidth',2)
plot(iter,monod_obj(real_d,alpha,param)*ones(size(iter)),':k','LineWidth',1.5)
xlabel('Iterácia')
ylabel('Hodnota účelovej funkcie J_{Monod}')
set(gca,'FontSize',12)
xlim([iter(1) iter(end)])
box on
hold off

figure(7)
hold on
plot(iter,P_bio_limit*ones(size(iter)),':k','LineWidth',1.5)
stairs(iter,P_us(2,:),'-b','LineWidth',2)
xlabel('Iterácia')
ylabel('Ustálený stav FIR modelu')
set(gca,'FontSize',12)
xlim([iter(1) iter(end)])
box on
hold off

figure(8)
stairs(iter,r_bio,'-b','LineWidth',2)
xlabel('Iterácia')
ylabel('Rád FIR modelu')
set(gca,'FontSize',12)
xlim([iter(1) iter(end)])
ylim([0 3.5])
box on

%% mensi rozptyl sumu pre koncentraciu biomasy
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

tf = 50;
Ts = 10; % number of measured samples through the simulation time

param = [mu_max nu Ks Ki Yx Yp s_in [] tf];

kmax = 20;

%Real (Monod) model optimum
alpha = .5;

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
nom_j = value(obj);

D = [nom_d; nom_d];

[x0,s0,p0] = get_Monod_ss(D(1),param);
[ix0,is0,ip0] = get_Haldane_ss(D(1),param);
init_c = [x0;s0;p0];
i_init_c = [ix0;is0;ip0];

param = [mu_max nu Ks Ki Yx Yp s_in D(1,1) tf];

for j = 1:1:kmax
    % hybrid model based on substrate 
    % Real data measurement
    [~,st,~] = generate_Monod_data(D(1,:),param,init_c,tf,Ts);
        
    % noise addition
    rng(1)
    power = 0.025;
    noise_sub = power*(2*rand(size(st)) - 1);

    m_st = st + noise_sub;
    
    % nominal model simulation
    [~,ist,~] = generate_Haldane_data(D(1,:),param,i_init_c,tf,Ts);
        
    %FIR model design
    ds = m_st - ist;
    ds = ds';
    u_sub = ones([Ts 1])*D(1,:);
    u_sub = u_sub(:)';
    
    u_sub = u_sub - nom_d;
    ds0 = s0 - is0;
    ds = ds - ds0;
    
    model_error = 2*power;
    
    [r_sub(j),P_sub] = gpe_fir_min_order(u_sub,ds,model_error,1000,1,0);
        
    % hybrid model optimum 
    a = sdpvar(1,1);
    obj = hybrid_obj_sub_corr(a,nom_d,ds0,alpha,param,P_sub);
    cons = [a >= 0.1;... 
            a <= mu_max];
    assign(a,0.1);
    ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
    optimize(cons,obj,ops);
    hyb_d_sub = value(a);
    hyb_sub_j = value(obj);
    
    fprintf(['\nHybrid model (Substrate) optimal D : ', num2str(hyb_d_sub), '\n'])
    
    % hybrid model based on biomass concentration
    % Real data measurement
    [xt,~,~] = generate_Monod_data(D(2,:),param,init_c,tf,Ts);
    
    % noise addition
    rng(0)
    power = 0.01;
    noise_bio = power*(2*rand(size(xt)) - 1);
    
    m_xt = xt + noise_bio;
    
    % nominal model simulation
    [ixt,~,~] = generate_Haldane_data(D(2,:),param,i_init_c,tf,Ts);
    
    % FIR model design
    dx = m_xt - ixt;
    u_bio = ones([Ts 1])*D(2,:);
    u_bio = u_bio(:)';
    dx = dx';
    
    u_bio = u_bio - nom_d;
    dx0 = x0 - ix0;
    dx = dx - dx0;
    
    model_error = 2*power;    
    
    [r_bio(j),P_bio] = gpe_fir_min_order(u_bio,dx,model_error,1000,0,0);
    
    % hybrid model optimum
    a = sdpvar(1,1);
    obj = hybrid_obj_bio_corr(a,nom_d,dx0,alpha,param,P_bio);
    cons = [a >= 0.1;...
            a <= mu_max];
    assign(a,0.1);
    ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
    optimize(cons,obj,ops);
    hyb_d_bio = value(a);
    hyb_bio_j = value(obj);
    
    fprintf(['\nHybrid model (Substrate) optimal D : ', num2str(hyb_d_bio), '\n'])
    
    D(:,j+1) = [hyb_d_sub; hyb_d_bio];
    P_us(:,j) = [sum(P_sub)*(D(1,j)-nom_d) + ds0; sum(P_bio)*(D(2,j)-nom_d) + dx0];
    optimal_d(:,j) = [hyb_d_sub; hyb_d_bio];
end

% steady state of P_sub for optimal D*
P_bio_limit = Monod_xt_ss(real_d,param) - Haldane_xt_ss(real_d,param);
P_sub_limit = Monod_st_ss(real_d,param) - Haldane_st_ss(real_d,param);

%% Figures
close all

figure(4)
iter = 1:1:kmax;
fs = monod_obj(optimal_d(1,:),alpha,param);
fx = monod_obj(optimal_d(2,:),alpha,param);

hold on
stairs(iter,fs,'-r','LineWidth',2)
stairs(iter,fx,'-b','LineWidth',2)
plot(iter,monod_obj(real_d,alpha,param)*ones(size(iter)),':k','LineWidth',1.5)
xlabel('Iterácia')
ylabel('Hodnota účelovej funkcie J_{Monod}')
xlim([iter(1) iter(end)])
box on
hold off


%% model natrenovany na viacerych skokovych zmenach
close all
clear all
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

tf = 100;
Ts = 10; % number of measured samples through the simulation time
t = linspace(0,tf,Ts);

param = [mu_max nu Ks Ki Yx Yp s_in [] tf];

%Real (Monod) model optimum
alpha = .5;

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
nom_j = value(obj);

D = linspace(nom_d,0.41,5);

[x0,s0,p0] = get_Monod_ss(D(1),param);
[ix0,is0,ip0] = get_Haldane_ss(D(1),param);
init_c = [x0;s0;p0];
i_init_c = [ix0;is0;ip0];

% Real data measurement
[~,st,~] = generate_Monod_data(D(1,:),param,init_c,tf,Ts);

% noise addition
rng(1)
power = 0.025;
noise_sub = power*(2*rand(size(st)) - 1);

m_st = st + noise_sub;

% nominal model simulation
[~,ist,~] = generate_Haldane_data(D(1,:),param,i_init_c,tf,Ts);

%FIR model design
ds = m_st - ist;
ds = ds';
u = ones([Ts 1])*D(1,:);
u = u(:)';

u = u - nom_d;
ds0 = s0 - is0;
ds = ds - ds0;

model_error = 7.97*power;

[r_sub,P_sub] = gpe_fir_min_order(u,ds,model_error,1000,1,0);

% hybrid model optimum
a = sdpvar(1,1);
obj = hybrid_obj_sub_corr(a,nom_d,ds0,alpha,param,P_sub);
cons = [a >= 0.1;...
    a <= mu_max];
assign(a,0.1);
ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
optimize(cons,obj,ops);
hyb_d_sub = value(a);
hyb_sub_j = value(obj);

fprintf(['\nHybrid model (Substrate) optimal D : ', num2str(hyb_d_sub), '\n'])

% steady state of P_sub for optimal D*
P_sub_limit = Monod_st_ss(real_d,param) - Haldane_st_ss(real_d,param);
P_us = sum(P_sub)*(hyb_d_sub - nom_d) + ds0;
%% Figures 
close all

figure(1)
C = get_matica_vstupov(u,r_sub);
fir_out = C*P_sub + ds0;
tt = linspace(0,length(D)*tf,length(st));

hold on
stairs(tt,ds + ds0,'or','MarkerFaceColor','r','LineWidth',1.5)
stairs(tt,fir_out,'-b','LineWidth',2)
xlabel('Cas [hod]')
ylabel('Rozdiel koncentracie substratu [gL^{-1}]')
legend('Namerané údaje','Výstup FIR','Location','Best')
set(gca,'FontSize',12)
box on
hold off

figure(2)
iter = 1:1:length(D);
hold on
plot(iter,monod_obj(nom_d,alpha,param)*ones(size(iter)),':b','LineWidth',1.5)
plot(iter,monod_obj(real_d,alpha,param)*ones(size(iter)),':k','LineWidth',1.5)
plot(iter,monod_obj(hyb_d_sub,alpha,param)*ones(size(iter)),'-r','LineWidth',2)
ylabel('Hodnota účelovej funkcie J_{Monod}')
legend('Nominálny model','Zariadenie','Hybridný model','Location','Best')
ylim([-0.472 -0.46])
set(gca,'XTick',[])
set(gca,'FontSize',12)
box on
hold off

%% Skokova zmena priamo na optimalne D zariadenia
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

tf = 50;
Ts = 10; % number of measured samples through the simulation time
t = linspace(0,tf,Ts);

param = [mu_max nu Ks Ki Yx Yp s_in [] tf];

%Real (Monod) model optimum
alpha = .5;

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
nom_j = value(obj);

D = [nom_d real_d];

[x0,s0,p0] = get_Monod_ss(D(1),param);
[ix0,is0,ip0] = get_Haldane_ss(D(1),param);
init_c = [x0;s0;p0];
i_init_c = [ix0;is0;ip0];

% Real data measurement
[~,st,~] = generate_Monod_data(D(1,:),param,init_c,tf,Ts);

% noise addition
rng(1)
power = 0.025;
noise_sub = power*(2*rand(size(st)) - 1);

m_st = st + noise_sub;

% nominal model simulation
[~,ist,~] = generate_Haldane_data(D(1,:),param,i_init_c,tf,Ts);

%FIR model design
ds = m_st - ist;
ds = ds';
u = ones([Ts 1])*D(1,:);
u = u(:)';

u = u - nom_d;
ds0 = s0 - is0;
ds = ds - ds0;

model_error = 2*power;

[r_sub,P_sub] = gpe_fir_min_order(u,ds,model_error,1000,1,0);

% hybrid model optimum
a = sdpvar(1,1);
obj = hybrid_obj_sub_corr(a,nom_d,ds0,alpha,param,P_sub);
cons = [a >= 0.1;...
    a <= mu_max];
assign(a,0.1);
ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
optimize(cons,obj,ops);
hyb_d_sub = value(a);
hyb_sub_j = value(obj);

fprintf(['\nHybrid model (Substrate) optimal D : ', num2str(hyb_d_sub), '\n'])

% steady state of P_sub for optimal D*
P_sub_limit = Monod_st_ss(real_d,param) - Haldane_st_ss(real_d,param);
P_us = sum(P_sub)*(real_d - nom_d) + ds0;
%% Figures 
close all

figure(1)
C = get_matica_vstupov(u,r_sub);
fir_out = C*P_sub + ds0;
tt = linspace(0,length(D)*tf,length(st));

hold on
plot(tt,ds + ds0,'or','MarkerFaceColor','r','LineWidth',1.5)
stairs(tt,fir_out,'-b','LineWidth',2)
xlabel('Čas [hod]')
ylabel('Rozdiel koncentrácie substrátu [gL^{-1}]')
legend('Namerané údaje','Výstup FIR','Location','Best')
set(gca,'FontSize',12)
box on
hold off

figure(2)
iter = 1:1:length(D);
hold on
plot(iter,monod_obj(nom_d,alpha,param)*ones(size(iter)),':b','LineWidth',1.5)
plot(iter,monod_obj(real_d,alpha,param)*ones(size(iter)),':k','LineWidth',1.5)
plot(iter,monod_obj(hyb_d_sub,alpha,param)*ones(size(iter)),'-r','LineWidth',2)
ylabel('Hodnota účelovej funkcie J_{Monod}')
legend('Nominálny model','Zariadenie','Hybridný model','Location','Best')
ylim([-0.472 -0.46])
set(gca,'XTick',[])
set(gca,'FontSize',12)
box on
hold off

dd = linspace(0,0.42,300);
fm = [];
fh = [];
for i = 1:1:length(dd)
   fm(i) = monod_obj(dd(i),alpha,param);
   fh(i) = hybrid_obj_sub_corr(dd(i),nom_d,ds0,alpha,param,P_sub);
   fn(i) = haldane_obj(dd(i),alpha,param);
end

figure
hold on
plot(dd,fm,'-b','LineWidth',2)
plot(dd,fn,'-r','LineWidth',2)
plot(dd,fh,'-g','LineWidth',2)
plot(real_d,real_j,'ko','LineWidth',2,'MarkerFaceColor','k')
plot(nom_d,nom_j,'ko','LineWidth',2,'MarkerFaceColor','k')
plot(hyb_d_sub,hyb_sub_j,'ko','LineWidth',2,'MarkerFaceColor','k')
xlabel('Rýchlosť riedenia D [hod^{-1}]')
ylabel('Hodnota účelovej funkcie J')
legend('Zariadenie','Nominálny model','Hybridný model')
xlim([dd(1) dd(end)])
box on
