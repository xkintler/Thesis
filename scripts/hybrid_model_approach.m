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

kmax = 30;

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
nom_j = obj;

D = [nom_d; nom_d];

[x0,s0,p0] = get_Monod_ss(D(1),param);
[ix0,is0,ip0] = get_Haldane_ss(D(1),param);
init_c = [x0;s0;p0];
i_init_c = [ix0;is0;ip0];

param = [mu_max nu Ks Ki Yx Yp s_in D(1,1) tf];

for j = 1:1:kmax
    % hybrid model based on biomass concentration
    % Real data measurement
    [xt,~,~] = generate_Monod_data(D(1,:),param,init_c,tf,Ts);
    
    % noise addition
    rng(0)
    power = 0.05;
    noise_bio = power*(2*rand(size(xt)) - 1);
    
    m_xt = xt + noise_bio;
    
    % nominal model simulation
    [ixt,~,~] = generate_Haldane_data(D(1,:),param,i_init_c,tf,Ts);
    
    % FIR model design
    dx = m_xt - ixt;
    u = ones([Ts 1])*D(1,:);
    u = u(:)';
    dx = dx';
    model_error = power;    
    
    [r_bio(j),P_bio] = gpe_fir_min_order(u,dx,model_error,1000,0,0);
    
    % hybrid model optimum
    a = sdpvar(1,1);
    obj = hybrid_obj_bio(a,alpha,param,P_bio);
    cons = [a >= 0.1;...
            a <= mu_max];
    assign(a,0.1);
    ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
    optimize(cons,obj,ops);
    hyb_d_bio = value(a);
    hyb_bio_j = value(obj);
    
    % figures
    if j == kmax
        tt = linspace(0,kmax*tf,kmax*Ts);
        C = get_matica_vstupov(u,r_bio(j));
        
        figure
        e = model_error.*ones(size(dx));
        hold on
        errorbar(tt,dx,-e,e,'LineStyle','none','Marker','o','Color','b','MarkerFaceColor','b','LineWidth',2)
        stairs(tt,C*0.2743,'-r','LineWidth',2)
        xlabel('Čas [hod]')
        ylabel('Rozdiel koncentrácií biomasy [gL^{-1}]')
        legend('Namerané údaje','Výstupy FIR modelu','Location','Best')
        set(gca,'FontSize',15)
        box on
        hold off

    end
    
    d = 0:0.001:0.4;
    for k = 1:1:length(d)
        fm(k) = monod_obj(d(k),alpha,param);
        fh(k) = hybrid_obj_bio(d(k),alpha,param,P_bio);
    end
    figure(2)
    hold on
    plot(d,fm,'b')
    plot(d,fh)
    plot(hyb_d_bio,hyb_bio_j,'o')
    box on
    hold off 
    
    % hybrid model based on substrate 
    % Real data measurement
    [~,st,~] = generate_Monod_data(D(2,:),param,init_c,tf,Ts);
        
    % noise addition
    rng(1)
    power = 0.025;
    noise_sub = power*(2*rand(size(st)) - 1);

    m_st = st + noise_sub;
    
    % nominal model simulation
    [~,ist,~] = generate_Haldane_data(D(2,:),param,i_init_c,tf,Ts);
        
    %FIR model design
    ds = m_st - ist;
    ds = ds';
    u = ones([Ts 1])*D(2,:);
    u = u(:)';
    model_error = power;
    
    [r_sub(j),P_sub] = gpe_fir_min_order(u,ds,model_error,1000,0,0);

    % figure
    if j == kmax
        C = get_matica_vstupov(u,r_sub(j));
        
        figure
        hold on
        plot(tt,ds)
        plot(tt,st-ist)
        stairs(tt,C*P_sub)
        box on
        hold off
        
        figure
        e = model_error.*ones(size(ds));
        hold on
        errorbar(tt,ds,-e,e,'LineStyle','none','Marker','o','Color','r','MarkerFaceColor','r','LineWidth',2)
        stairs(tt,C*P_sub,'-b','LineWidth',2)
        xlabel('Čas [hod]')
        ylabel('Koncentrácia substrátu [gL^{-1}]')
        legend('Namerané údaje','Výstupy FIR modelu','Location','Best')
        set(gca,'FontSize',15)
        box on
        hold off
    end
    
    % hybrid model optimum 
    a = sdpvar(1,1);
    obj = hybrid_obj_sub_corr(a,alpha,param,P_sub);
    cons = [a >= 0.1;... 
            a <= mu_max];
    assign(a,0.1);
    ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
    optimize(cons,obj,ops);
    hyb_d_sub = value(a);
    hyb_sub_j = obj;

%     % ARX model design
%     % arx identification toolbox
%     io_data = iddata(ds,u');
%     na = 3;
%     nb = 2;
%     nk = 1;
%     orders = [na nb nk];
%     sys = arx(io_data,orders);
%     opt_A = sys.A;
%     opt_B = sys.B;
    
    D(:,j+1) = [hyb_d_bio; hyb_d_sub];
%     if j >= 2
%         if sum(diff(D(1,:))) < 0.005
%             D(1,j+1) = 0.005 + D(1,j);
%         else
%             D(1,j+1) = hyb_d_bio;
%         end
%         
%         if sum(diff(D(2,:))) < 0.005
%             D(2,j+1) = 0.005 + D(2,j);
%         else
%             D(2,j+1) = hyb_d_bio;
%         end
%     end
   
    optimal_d(:,j) = [real_d; nom_d; hyb_d_bio; hyb_d_sub];
    opt_j(:,j) = [real_j; nom_j; hyb_bio_j; hyb_sub_j];
    P_us(:,j) = [sum(P_bio)*hyb_d_bio sum(P_sub)*(hyb_d_sub-nom_d)];

end

fprintf(['\nReal (Monod) model optimal D : ', num2str(real_d), '\n'])
fprintf(['\nNominal (Haldane) model optimal D : ', num2str(nom_d), '\n'])
fprintf(['\nHybrid model (Biomass) optimal D : ', num2str(hyb_d_bio), '\n'])
fprintf(['\nHybrid model (Substrate) optimal D : ', num2str(hyb_d_sub), '\n'])

% steady state of P_sub for optimal D*
D = real_d;
sin = s_in;
sm = -(D.*Ks)./(D - mu_max);
sn = -(Ki.*(D - mu_max + (((D.^2).*Ki - 4.*(D.^2).*Ks + Ki*mu_max^2 - 2.*D.*Ki.*mu_max)./Ki).^(1/2)))./(2.*D);
P_sub_limit = sm - sn;

% steady state of P_bio for optimal D*
s = sm;
xm = -(D.*(s - sin))./(D./Yx + nu/Yp);
s = sn;
xn = -(D.*(s - sin))./(D./Yx + nu/Yp);
P_bio_limit = xm - xn;

figure
iter = 1:1:j;
hold on
plot(iter,P_sub_limit*ones(size(iter)),':k','LineWidth',1.5)
stairs(iter,P_us(2,:),'-r','LineWidth',2)
xlabel('Iterácia')
ylabel('Ustálený stav FIR_{sub}')
set(gca,'FontSize',15)
xlim([iter(1) iter(end)])
box on
hold off

figure
stairs(iter,r_sub,'-r','LineWidth',2)
xlabel('Iterácia')
ylabel('Rád FIR_{sub} modelu')
set(gca,'FontSize',15)
xlim([iter(1) iter(end)])
box on

figure
hold on
stairs(tt,m_xt,'-b','LineWidth',2)
stairs(tt,m_st,'-r','LineWidth',2)
xlabel('Čas [hod]')
ylabel('Koncentrácia x(t), s(t) [hod]')
legend('Biomasa','Substrát','Location','Best')
box on
hold off
close

figure
hold on
stairs(tt,ixt,'-b','LineWidth',2)
stairs(tt,ist,'-r','LineWidth',2)
xlabel('Čas [hod]')
ylabel('Koncentrácia x(t), s(t) [hod]')
legend('Biomasa','Substrát','Location','Best')
box on
hold off
close 

figure
hold on
plot(iter,P_bio_limit*ones(size(iter)),':k','LineWidth',1.5)
stairs(iter,P_us(1,:),'b','LineWidth',2)
xlabel('Iterácia')
ylabel('Ustálený stav FIR_{bio}')
box on
hold off

figure
fb = monod_obj(optimal_d(3,:),alpha,param);
fs = monod_obj(optimal_d(4,:),alpha,param);
hold on
stairs(iter,fb,'-b','LineWidth',2)
stairs(iter,fs,'-r','LineWidth',2)
plot(iter,monod_obj(real_d,alpha,param)*ones(size(iter)),':k','LineWidth',1.5)
xlabel('Iterácia')
ylabel('Hodnota účelovej funkcie J_{Monod}')
legend('Biomasa','Substrát','Location','Best')
box on
hold off

figure
hxt = [];
hold on
plot(xt)
plot(ixt)
for i = 1:1:kmax
    hxt = [hxt ixt((1+(i-1)*Ts):(i*Ts))+P_us(1,i)];
end
plot(hxt)
close

figure
hold on
stairs(tt,m_xt,'-b','LineWidth',2)
stairs(tt,ixt,'-c','LineWidth',2)
xlabel('Čas [hod]')
ylabel('Koncentrácia biomasy [gL^{-1}]')
legend('Zariadenie (Monod)','Nominálny model (Haldane)','Location','Best')
axis([tt(1) tt(end) 4.30 4.51])
set(gca,'FontSize',15)
box on
close

figure
hold on
stairs(tt,m_st,'-r','LineWidth',2)
stairs(tt,ist,'-m','LineWidth',2)
xlabel('Čas [hod]')
ylabel('Koncentrácia substrátu [gL^{-1}]')
legend('Zariadenie (Monod)','Nominálny model (Haldane)','Location','Best')
axis([tt(1) tt(end) 2.85 3.55])
set(gca,'FontSize',15)
box on
close



%% Trenovanie na niekolkych skokovych zmenach
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
nom_j = obj;

D = linspace(nom_d,0.41,5);
tf = 100;
optimal_d = [];

[x0,s0,p0] = get_Monod_ss(D(1),param);
[ix0,is0,ip0] = get_Haldane_ss(D(1),param);
init_c = [x0;s0;p0];
i_init_c = [ix0;is0;ip0];

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
u = ones([Ts 1])*D(1,:);
u = u(:)';
model_error = power;

[r_sub,P_sub] = gpe_fir_min_order(u,ds,model_error,1000,0,0);

% hybrid model optimum
a = sdpvar(1,1);
obj = hybrid_obj_sub(a,alpha,param,P_sub);
cons = [a >= 0.1;...
    a <= mu_max];
assign(a,0.1);
ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
optimize(cons,obj,ops);
hyb_d_sub = value(a);
hyb_sub_j = obj;

tt = linspace(0,length(D)*tf,length(st));

figure
e = model_error.*ones(size(ds));
C = get_matica_vstupov(u,r_sub);
hold on
plot(tt,ds,'Color','r','Marker','o','MarkerFaceColor','r','LineWidth',1.5,'LineStyle','none')
stairs(tt,C*P_sub,'-b','LineWidth',2)
xlabel('Čas [hod]')
ylabel('Rozdiel koncentrácií substrátu \Delta_{s} [gL^{-1}]')
legend('Trénovacie údaje','FIR_{sub}','Location','Best')
set(gca,'FontSize',15)
box on
hold off

D = [D hyb_d_sub];
n = length(D);
kmax = 10;

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
    u = ones([Ts 1])*D(1,:);
    u = u(:)';
    model_error = power;
    
    [r_sub(j),P_sub] = gpe_fir_min_order(u,ds,model_error,1000,0,0);

   % hybrid model optimum 
    a = sdpvar(1,1);
    obj = hybrid_obj_sub(a,alpha,param,P_sub);
    cons = [a >= 0.1;... 
            a <= mu_max];
    assign(a,0.1);
    ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
    optimize(cons,obj,ops);
    hyb_d_sub = value(a);
    hyb_sub_j = obj;

    D(:,j+n) = hyb_d_sub;
  
    optimal_d(:,j) = [real_d; nom_d; hyb_d_sub];
    opt_j(:,j) = [real_j; nom_j; hyb_sub_j];
    P_us(:,j) = sum(P_sub)*hyb_d_sub;

end

figure
iter = 0:1:kmax;
fs = monod_obj([D(6) optimal_d(3,:)],alpha,param);
hold on
stairs(iter,fs,'-r','LineWidth',2)
plot(iter,monod_obj(real_d,alpha,param)*ones(size(iter)),':k','LineWidth',1.5)
xlabel('Iterácia')
ylabel('Hodnota účelovej funkcie J_{Monod}')
set(gca,'FontSize',15)
box on
hold off

%% Skokova zmena priamo na optimalnu hodnotu D zariadenia
% porovnanie ucelovych funkcii
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

kmax = 10;

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
nom_j = obj;

D = [nom_d real_d];

[x0,s0,p0] = get_Monod_ss(D(1),param);
[ix0,is0,ip0] = get_Haldane_ss(D(1),param);
init_c = [x0;s0;p0];
i_init_c = [ix0;is0;ip0];

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
u = ones([Ts 1])*D(1,:);
u = u(:)';
model_error = power;

[r_sub,P_sub] = gpe_fir_min_order(u,ds,model_error,1000,0,0);

% hybrid model optimum
a = sdpvar(1,1);
obj = hybrid_obj_sub(a,alpha,param,P_sub);
cons = [a >= 0.1;...
    a <= mu_max];
assign(a,0.1);
ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
optimize(cons,obj,ops);
hyb_d_sub = value(a);
hyb_sub_j = value(obj);

P_us = [0 sum(P_sub)*hyb_d_sub];

dd = linspace(0,0.42,300);
fm = [];
fh = [];
for i = 1:1:length(dd)
   fm(i) = monod_obj(dd(i),alpha,param);
   fh(i) = hybrid_obj_sub(dd(i),alpha,param,P_sub);
end

% steady state of P_sub for optimal D*
D = real_d;
sin = s_in;
sm = -(D.*Ks)./(D - mu_max);
sn = -(Ki.*(D - mu_max + (((D.^2).*Ki - 4.*(D.^2).*Ks + Ki*mu_max^2 - 2.*D.*Ki.*mu_max)./Ki).^(1/2)))./(2.*D);
P_sub_limit = sm - sn;

figure
hold on
plot(dd,fm,'-b','LineWidth',2)
plot(dd,fh,'-r','LineWidth',2)
plot(real_d,real_j,'ko','LineWidth',2,'MarkerFaceColor','k')
plot(hyb_d_sub,hyb_sub_j,'ko','LineWidth',2,'MarkerFaceColor','k')
xlabel('Rýchlosť riedenia D [hod^{-1}]')
ylabel('Hodnota účelovej funkcie J')
legend('Zariadenie','Hybridný model','Optimum')
xlim([dd(1) dd(end)])
set(gca,'FontSize',15)
box on

iter = [1 2];

figure
hold on
plot(iter,monod_obj(nom_d,alpha,param)*ones(size(iter)),':b','LineWidth',1.5)
plot(iter,monod_obj(real_d,alpha,param)*ones(size(iter)),':k','LineWidth',1.5)
plot(iter,monod_obj(hyb_d_sub,alpha,param)*ones(size(iter)),'-r','LineWidth',2)
ylabel('Hodnota účelovej funkcie J_{Monod}')
legend('Nominálny model','Zariadenie','Hybridný model','Location','Best')
set(gca,'XTick',[])
set(gca,'FontSize',15)
box on
hold off

figure
hold on
plot(iter,P_sub_limit*ones(size(iter)),':k','LineWidth',1.5)
stairs(iter,P_us(1,:),'b','LineWidth',2)
xlabel('Iterácia')
ylabel('Ustálený stav FIR_{sub}')
box on
hold off
