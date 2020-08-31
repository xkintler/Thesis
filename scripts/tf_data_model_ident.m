clear all
close all
clc

addpath

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
xlabel('Cas [s]')
ylabel('Vstupny signal')
set(gca,'FontSize',15)
box on 
grid on



figure
hold on
plot(t,output,':','LineWidth',1.5,'Color',[105/255 105/255 105/255])
plot(t,y_measured,'ok','LineWidth',1.3,'MarkerFaceColor','k')
xlabel('Cas [s]')
ylabel('Vystupny signal')
set(gca,'FontSize',15)
box on 
grid on

%% FIR model identification
model_error = 0.8;
max_iter = 100;
[r_min,theta,~] = gpe_fir_min_order(input,y_measured,model_error,max_iter);
[theta_min, theta_max] = gpe_fir_min_max_bound(r_min, input, y_measured, model_error);

fir_min = generate_FIR_output(theta_min,input);
fir_max = generate_FIR_output(theta_max,input);

%% ARX model identification
model_error = 0.8;
y0 = 0;
max_iter = 10;
param_bound = 20;

[na_min,nb_min,A,B,info] = gpe_arx_min_order(input,y_measured,model_error,y0,max_iter,param_bound);
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
xlabel('Cas [s]')
ylabel('Vystupny signal')
axis([0 20 -0.5 10])
%legend('Namerané údaje','FIR_{min}','FIR_{max}','Location','Best')
box on

figure
hold on
errorbar(t,y_measured,-e,e,'Color','k','Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',1.3,'LineStyle','None')
stairs(t,arx_min,'-r','LineWidth',2)
stairs(t,arx_max,'-b','LineWidth',2)
xlabel('Cas [s]')
ylabel('Vystupny signal')
axis([0 20 -0.5 10])
%legend('Namerané údaje','ARX_{min}','ARX_{max}','Location','Best')
box on
