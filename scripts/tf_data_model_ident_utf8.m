clear all
close all
clc

addpath('e:\skola\semestralny projekt\1')
addpath('funs')

%% measured data
K = 8;
T = 2;
Gs = tf(K,[T 1]);

tf = 20;
Ts = 0.5;
t = linspace(0,tf,tf/Ts);

input = ones(size(t));
output = step(Gs,t)';

% noise
rng(0)
power = 0.8;
noise = (2*rand(size(t)) - 1);

y_measured = power.*noise + output;figure
hold on
stairs(t,input,'-k','LineWidth',2,'MarkerFaceColor','k')
xlabel('Čas [s]')
ylabel('Vstupný signál')
set(gca,'FontSize',15)
box on 

figure
hold on
plot(t,output,':','LineWidth',1.5,'Color',[105/255 105/255 105/255])
plot(t,y_measured,'ok','LineWidth',1.3,'MarkerFaceColor','k')
xlabel('Čas [s]')
ylabel('Výstupný signál')
legend('Bez šumu','So šumom','Location','Best')
set(gca,'FontSize',15)
box on 

%% FIR model identification
model_error = 0.8;
max_iter = 100;
[r_min,theta] = gpe_fir_min_order(input,y_measured,model_error,max_iter,1,0);
[theta_min, theta_max] = gpe_fir_min_max_bound(r_min, input, y_measured, model_error,1);

fir_min = generate_FIR_output(theta_min,input);
fir_max = generate_FIR_output(theta_max,input);

%% ARX model identification
model_error = 0.8;
y0 = 0;
max_iter = 10;
param_bound = 20;

gpe_out = gpe_arx_min_order(input,y_measured,model_error,y0,max_iter,param_bound,0);
na_min = gpe_out.na;
nb_min = gpe_out.nb;

[A_min,B_min,A_max,B_max] = gpe_arx_min_max_bound(na_min,nb_min,input,y_measured,y0,model_error,param_bound);

arx_min = generate_ARX_output(A_min,B_min,input,y0);
arx_max = generate_ARX_output(A_max,B_max,input,y0);


%% Results
e = model_error*ones(size(t));

figure
hold on
errorbar(t,y_measured,-e,e,'Color','k','Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',1.3,'LineStyle','None')
stairs(t,fir_min,'-r','LineWidth',2)
stairs(t,fir_max,'-b','LineWidth',2)
xlabel('Čas [s]')
ylabel('Výstupný signál')
axis([0 20 -0.5 10])
legend('Namerané údaje','FIR_{min}','FIR_{max}','Location','Best')
set(gca,'FontSize',15)
box on

figure
hold on
errorbar(t,y_measured,-e,e,'Color','k','Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',1.3,'LineStyle','None')
stairs(t,arx_min,'-r','LineWidth',2)
stairs(t,arx_max,'-b','LineWidth',2)
xlabel('Čas [s]')
ylabel('Výstupný signál')
axis([0 20 -0.5 10])
legend('Namerané údaje','ARX_{min}','ARX_{max}','Location','Best')
set(gca,'FontSize',15)
box on

%% Pareto fornt FIR
clear all
close all
clc

addpath('e:\skola\semestralny projekt\1')
addpath('funs')

K = 8;
T = 2;
Gs = tf(K,[T 1]);

tf = 20;
Ts = 0.5;
t = linspace(0,tf,tf/Ts);

input = ones(size(t));
output = step(Gs,t)';

model_error = 1;
max_iter = 100;
imax = 100;

for i = 1:1:imax
    k = 1;
    % noise
    rng(i)
    power = 0.9;
    noise = (2*rand(size(t)) - 1);

    y_measured = power.*noise + output;
    y_measured = y_measured';
    
    r_min = gpe_fir_min_order(input,y_measured,model_error,max_iter,0,0);
    
    for r = r_min:1:r_min+4
       C = get_matica_vstupov(input,r);
       P_lsm = (C'*C)\C'*y_measured;
       [P_min, P_max] = gpe_fir_min_max_bound(r, input, y_measured, model_error,0);

       fir_out = C*P_lsm;
       fir_min = C*P_min;
       fir_max = C*P_max;

       acc = sqrt(sum((y_measured-fir_out).^2))/length(y_measured);
       oest = sqrt(sum(max([(fir_out-fir_max).^2, (fir_out-fir_min).^2]')))/length(y_measured);
       
       fir_pareto_data(k,:) = [r acc oest];
       
       k = k + 1;
    end
    pareto_cell{i} = fir_pareto_data;
end

%% data conditioning
addpath('data')
addpath('funs')
load('ex_fir_pareto_data.mat')

fir_pareto_data = [];

for i = 1:1:imax
   fir_pareto_data = [fir_pareto_data; pareto_cell{i}]; 
end
[avg_mat,cnt] = average_if_matrix(fir_pareto_data);

rel_data = fir_pareto_data(:,2:end)./avg_mat(5,2:end);
rel_data = average_if_matrix([fir_pareto_data(:,1), rel_data]);

acc = rel_data(5:9,2);
oest = rel_data(5:9,3);
r = rel_data(5:9,1);

% acc = avg_mat(5:9,2);
% oest = avg_mat(5:9,3);
% r = avg_mat(5:9,1);

figure
hold on
plot(acc,oest,'bo','MarkerFaceColor','b')
for i = 1:1:length(acc)
    text(acc(i)-0.0015,oest(i)-0.005,num2str(r(i)))
end
xlabel('Presnosť odhadu')
ylabel('Maximálny rozptyl odhadu')
xlim([0.92 1.005])
ylim([0.98 1.15])
set(gca,'FontSize',12)
grid on
box on

%% Pareto fornt ARX
clear all
close all
clc

addpath('e:\skola\semestralny projekt\1')
addpath('funs')

K = 8;
T = 2;
Gs = tf(K,[T 1]);

tf = 20;
Ts = 0.5;
t = linspace(0,tf,tf/Ts);

zrs = zeros(1,5);

input = [zrs, ones(size(t))];
output = [zrs, step(Gs,t)'];

model_error = 1;
max_iter = 100;
y0 = 0;
imax = 20;

for i = 1:1:imax
    k = 1;
    % noise
    rng(i)
    power = 0.9;
    noise = [zrs, (2*rand(size(t)) - 1)];

    y_measured = power.*noise + output;
    y_measured = y_measured';
    
    arx_gpe_out = gpe_arx_min_order(input,y_measured,model_error,y0,max_iter,20,0);
    na = arx_gpe_out.na;
    nb = arx_gpe_out.nb;
    nk = 1;
    for r = na:1:na+3
       io_data = iddata(y_measured,input');
       orders = [r nb nk];
       sys = arx(io_data,orders);
       A_lsm = sys.A(2:end)';
       B_lsm = sys.B(2:end)';
       
       [A_min,B_min,A_max,B_max] = gpe_arx_min_max_bound(r,nb,input,y_measured,y0,model_error,5);

       arx_out = generate_ARX_output(A_lsm,B_lsm,input,y0);
       arx_min = generate_ARX_output(A_min,B_max,input,y0);
       arx_max = generate_ARX_output(A_max,B_max,input,y0);

       acc = sqrt(sum((y_measured-arx_out).^2))/length(y_measured);
       oest = sqrt(sum(max([(arx_out-arx_max).^2, (arx_out-arx_min).^2]')))/length(y_measured);
       
       arx_pareto_data(k,:) = [r acc oest];
       
       k = k + 1;
    end
    pareto_cell{i} = arx_pareto_data;
end
cd data
save('ex_arx_pareto_data.mat')
%% data conditioning
addpath('data')
load('ex_arx_pareto_data.mat')

arx_pareto_data = [];

for i = 1:1:imax
   arx_pareto_data = [arx_pareto_data; pareto_cell{i}]; 
end
[avg_mat,cnt] = average_if_matrix(arx_pareto_data);

rel_data = arx_pareto_data(:,2:end)./avg_mat(1,2:end);
rel_data = average_if_matrix([arx_pareto_data(:,1), rel_data]);

acc = rel_data(:,2);
oest = rel_data(:,3);
r = rel_data(:,1);

% acc = avg_mat(5:9,2);
% oest = avg_mat(5:9,3);
% r = avg_mat(5:9,1);

figure
hold on
plot(acc,oest,'bo','MarkerFaceColor','b')
for i = 1:1:length(acc)
    text(acc(i)-0.004,oest(i)-0.2,['1/' num2str(r(i))])
end
xlabel('Relatívna presnosť odhadu')
ylabel('Maximálny relatívny rozptyl odhadu')
xlim([0.85 1.02])
ylim([0 8])
set(gca,'FontSize',15)
grid on
box on

%% kontrola
clear all
close all
clc

addpath('funs')
addpath('data')
load('ex_arx_pareto_data.mat')

rng(103)
power = 0.9;
noise = (2*rand(size(t)) - 1);

noise = [zeros(1,5), noise];
input = [zeros(1,5), input];
output = [zeros(1,5), output];

y_measured = power.*noise + output;
y_measured = y_measured';

na = [1 2 3 4];
nb = 1;
nk = 1;


io_data = iddata(y_measured,input');

figure(1)
hold on
plot(y_measured,'ok','MarkerFaceColor','k')
for i = 1:1:length(na)
    
    orders = [na(i) nb nk];
    sys = arx(io_data,orders);
    A_lsm = sys.A(2:end)';
    B_lsm = sys.B(2:end)';
    
    disp(sys.B)
    disp(sys.A)
    
    arx_out = generate_ARX_output(A_lsm,B_lsm,input,y0);
    
    delta(i,:) = (y_measured-arx_out).^2;
    
    stairs(arx_out,'LineWidth',1.5)
end
legend('','1','2','3','4')
hold off

figure(2)
hold on
stairs(delta(1,:),'LineWidth',1.5)
stairs(delta(2,:),'LineWidth',1.5)
stairs(delta(3,:),'LineWidth',1.5)
stairs(delta(4,:),'LineWidth',1.5)
legend('1','2','3','4')
