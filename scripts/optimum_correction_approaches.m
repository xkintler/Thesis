clear all
close all
clc


% statistika 


addpath('e:\skola\semestralny projekt\1')
addpath('funs')
addpath('data')

F = 1; %l/h
V = 3.33; %l
D = [F/V; F/V; F/V; F/V];
mu_max = 0.53; %h-1
nu = 0.5; %h-1 
Ks = 1.2; %g/l
Ki = 70; %g/l
Yx = .4;
Yp = 1;
s_in = 20; %g/l

Ts = 10; % number of measured samples through the simulation time
tf = 50;
alpha = .5;

params = [mu_max nu Ks Ki Yx Yp s_in [] tf];

%Real (Monod) model optimum
a = sdpvar(1,1);
obj = monod_obj(a,alpha,params);
cons = [a >= 0; a <= mu_max];
assign(a,0.1);
ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
optimize(cons,obj,ops);
real_d = value(a);
real_j = value(obj);


% Nominal (Haldane) model optimum
a = sdpvar(1,1);
obj = haldane_obj(a,alpha,params);
cons = [a >= 0; a <= mu_max];
assign(a,0.1);
ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
optimize(cons,obj,ops);
nom_d = value(a);
nom_j = value(obj);

kmax = 20;
imax = 20;
for i = 1:1:imax
    D = [nom_d; nom_d; nom_d];
    
    est_mu = 0.5;
    est_Ks = 1;
    est_Ki = 50;
    opts = optimset('MaxFunEvals',1e5,'MaxIter',1e3);
    
    grad_store = grad_monod_obj(nom_d,alpha,params);
    
    [x0,s0,p0] = get_Monod_ss(D(1),params);
    [ix0,is0,ip0] = get_Haldane_ss(D(1),params);
    tic
    for j = 1:1:kmax
        % Parameter estimation approach (2-step approach)
        % Real data measurement
        [~,st,~] = generate_Monod_data(D(1,:),params,[x0;s0;p0],tf,Ts);

        % Noise generation
        rng(i-1)
        power = 0.025;
        noise_sub = power*(2*rand(size(st)) - 1);

        m_st = st + noise_sub;
        p1 = {mu_max nu Ks Ki Yx Yp s_in D(1,:) tf Ts};

        [X,FVAL,EXITFLAG,OUTPUT] = fminsearch(@(a) getValueOfCostFun(a,m_st,[x0;s0;p0],p1),[est_mu;est_Ks;est_Ki],opts);
        est_mu = X(1);
        est_Ks = X(2);
        est_Ki = X(3);
        est_params = [est_mu nu est_Ks est_Ki Yx Yp s_in D(j) tf];

        % Optimal solution based on parameter estimation
        a = sdpvar(1,1);
        obj = haldane_obj(a,alpha,est_params);
        cons = [a >= 0; a <= mu_max; a <= est_mu];
        assign(a,0.1);
        ops = sdpsettings('usex0',1,'solver','baron','verbose',1);
        optimize(cons,obj,ops);
        pe_d = value(a);
        pe_j = value(obj);

        % Modifier adaptation scheme approach
        [xt,~,~] = generate_Monod_data(D(2,:),params,[x0;s0;p0],tf,Ts);

        rng(i)
        power = 0.03;
        noise_bio = power*(2*rand(size(xt)) - 1);

        m_xt = xt + noise_bio;

        c = 0.8;
       
        if j == 1
            lam = -0.07;
        else
            a_xt0 = mean(m_xt(1,end-3:end));
            a_xt1 = mean(m_xt(1,end-3-Ts:end-Ts));
            
            J0 = D(2,j).*(1 - alpha*a_xt0);
            J1 = D(2,j-1).*(1 - alpha*a_xt1);
            real_grad_app = (J0 - J1)/(D(2,j) - D(2,j-1));

            % ohranicenie odhadu gradientu
            if j > 1 && abs((real_grad_app - grad_store(j-1))/grad_store(j-1)) > 2
                if real_grad_app > 0
                    real_grad_app = 1;
                elseif real_grad_app < 0
                    real_grad_app = -1;
                end
            end
            grad_store(j) = real_grad_app;

            temp_lam = real_grad_app - grad_haldane_obj(D(2,j),alpha,params);
            lam = c*lam + (1 - c)*temp_lam;

            figure(4)
            hold on
            plot(j,grad_monod_obj(D(2,j),alpha,params),'ok')
            plot(j,grad_store(j),'ro')
            box on
            hold off
        end

        a = sdpvar(1,1);
        obj = modifier_adapt(a,alpha,params,lam);
        cons = [a >= 0.1; a <= mu_max];
        assign(a,0.1);
        ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
        optimize(cons,obj,ops);
        mas_d = value(a);
        mas_j = value(obj);

        % Hybrid model based on substrate concetration
        % plant data generation
        [~,st,~] = generate_Monod_data(D(3,:),params,[x0;s0;p0],tf,Ts);
        m_st = st + noise_sub;

        [~,ist,~] = generate_Haldane_data(D(3,:),params,[ix0;is0;ip0],tf,Ts);

        ds = m_st - ist;
        ds = ds';
        u = ones([Ts 1])*D(3,:);
        u = u(:)';
        
        u = u - nom_d;
        ds0 = s0 - is0;
        ds = ds - ds0;
        
        model_error = 2*0.025;

        [r_sub(j),P_sub] = gpe_fir_min_order(u,ds,model_error,1000,0,0);

        a = sdpvar(1,1);
        obj = hybrid_obj_sub_corr(a,nom_d,ds0,alpha,params,P_sub);
        cons = [a >= 0; a <= mu_max];
        assign(a,0.1);
        ops = sdpsettings('usex0',1,'solver','baron','verbose',0);
        optimize(cons,obj,ops);
        hyb_d_sub = value(a);
        hyb_sub_j = value(obj);

        D(:,j+1) = [pe_d; mas_d; hyb_d_sub];
        optimal_d(:,j) = [real_d; nom_d; pe_d; mas_d; hyb_d_sub];
        opt_j(:,j) = [real_j; nom_j; pe_j; mas_j; hyb_sub_j];

    end
    cell_data{i} = [optimal_d; opt_j];
    toc
end
%% data conditioning
d_opt_pe = [];
d_opt_mas = [];
d_opt_hyb = [];
j_opt_pe = [];
j_opt_mas = [];
j_opt_hyb = [];
for i = 1:1:imax
   d_opt_pe = [d_opt_pe; cell_data{i}(3,:)];
   d_opt_mas = [d_opt_mas; cell_data{i}(4,:)];
   d_opt_hyb = [d_opt_hyb; cell_data{i}(5,:)];
   
   j_opt_pe = [j_opt_pe; cell_data{i}(8,:)];
   j_opt_mas = [j_opt_mas; cell_data{i}(9,:)];
   j_opt_hyb = [j_opt_hyb; cell_data{i}(10,:)];
end

optimal_d = [real_d.*ones([1 kmax]); nom_d.*ones([1 kmax]); mean(d_opt_pe); mean(d_opt_mas); mean(d_opt_hyb)];
optimal_j = [real_j.*ones([1 kmax]); nom_j.*ones([1 kmax]); mean(j_opt_pe); mean(j_opt_mas); mean(j_opt_hyb)];

%% 
clear
load('opt_approaches_data.mat')
%% figure
close all

figure(1)
iter = 1:1:kmax;
hold on
plot(iter,monod_obj(optimal_d(1,:),alpha,params),':k','LineWidth',1.5)
plot(iter,monod_obj(optimal_d(2,:),alpha,params),':r','LineWidth',1.5)
stairs(iter,monod_obj(optimal_d(3,:),alpha,params),'b','LineWidth',2)
stairs(iter,monod_obj(optimal_d(4,:),alpha,params),'g','LineWidth',2)
stairs(iter,monod_obj(optimal_d(5,:),alpha,params),'m','LineWidth',2)
ylabel('Hodnota účelovej funkcie J_{Monod}')
xlabel('Iterácia')
%legend('Zariadenie','Nominálny model','PE','MAS','HYB','Location','Best')
xlim([iter(1) iter(end)])
ylim([-0.471 -0.458])
set(gca,'FontSize',12)
box on

figure(2)
hold on
plot(iter,optimal_d(1,:),':k','LineWidth',1.5)
plot(iter,optimal_d(2,:),':r','LineWidth',1.5)
stairs(iter,optimal_d(3,:),'b','LineWidth',2)
stairs(iter,optimal_d(4,:),'g','LineWidth',2)
stairs(iter,optimal_d(5,:),'c','LineWidth',2)
ylabel('Rýchlosť riedenia D [hod^{-1}]')
xlabel('Iterácia')
% legend('Zariadenie','Nominálny model','PE','MAS','HYB','Location','Best')
xlim([iter(1) iter(end)])
set(gca,'FontSize',15)
box on

figure(3)
hold on
plot(iter,opt_j(1,:),':k')
plot(iter,opt_j(2,:),':r')
stairs(iter,opt_j(3,:),'b')
stairs(iter,opt_j(4,:),'g')
stairs(iter,opt_j(5,:),'c')
xlabel('Iteration')
ylabel('Objective function value J*')
close(3)

figure(5)
hold on
yyaxis left
plot(iter,optimal_d(1,:),':','LineWidth',1.5)
plot(iter,optimal_d(2,:),'.-.','LineWidth',1.5)
stairs(iter,optimal_d(4,:),'-','LineWidth',2)
stairs(iter,optimal_d(5,:),'--','LineWidth',2)
ylabel('Rýchlosť riedenia D [hod^{-1}]')
xlabel('Iterácia')

yyaxis right
stairs(iter,optimal_d(3,:),'-','LineWidth',2)
% legend('Zariadenie','Nominálny model','PE','MAS','HYB','Location','Best')
xlim([iter(1) iter(end)])
set(gca,'FontSize',15)
box on