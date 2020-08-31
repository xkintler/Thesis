function [A_min,B_min,A_max,B_max] = gpe_arx_min_max_bound(na,nb,u,y,y0,model_error,bound)
    % tato funkcia hlada kombinaciu parametrov A,B, ktora poskytuje
    % minimalnu a maximalnu hranicu vystupu ARX modelu. Vyber je
    % realizovany na zaklade min a max ustalenej hodnoty modelu. 
    % !!! DOLEZITE !!! MINIMALNA HODNOTA PARAMETRA NEZARUCUJE MINIMALNY
    % USTALENY STAV !!!
    
    if size(u,1) < size(u,2)
        u = u';
    end
    
    if size(y,1) < size(y,2)
        y = y';
    end
    
    for s = 1:1:na+nb
        for j = 1:1:2
            if s == 1
                %Objective
                fun = @(x) ((-1)^j)*x(s);
                %Nonlinear Constraints
                nlcon = @(x) get_gpe_arx_constrains(x,u,y0,na,nb);

                cl = zeros([length(u)-1 1]);
                cu = zeros([length(u)-1 1]);
                %Bounds
                lb = [(-bound)*ones([na+nb 1]); y(2:end)-model_error];
                ub = [bound*ones([na+nb 1]); y(2:end)+model_error];

                %Solve
                x = baron(fun,[],[],[],lb,ub,nlcon,cl,cu);
                
                if j == 2
                    min_P = x(1:na+nb);
                    sum_min_P = sum(min_P(na+1:na+nb,1))/(sum(min_P(1:na,1))+1);
                else
                    max_P = x(1:na+nb);
                    sum_max_P = sum(max_P(na+1:na+nb,1))/(sum(max_P(1:na,1))+1);
                end    
            else
               %Objective
                fun = @(x) ((-1)^j)*x(s);
                %Nonlinear Constraints
                nlcon = @(x) get_gpe_arx_constrains(x,u,y0,na,nb);

                cl = zeros([length(u)-1 1]);
                cu = zeros([length(u)-1 1]);
                %Bounds
                lb = [(-bound)*ones([na+nb 1]); y(2:end)-model_error];
                ub = [bound*ones([na+nb 1]); y(2:end)+model_error];

                %Solve
                x = baron(fun,[],[],[],lb,ub,nlcon,cl,cu);
                new_P = x(1:na+nb);
                new_sum_P = sum(new_P(na+1:na+nb,1))/(sum(new_P(1:na,1))+1);
                
                if j == 2
                    if new_sum_P < sum_min_P
                        min_P = new_P;
                        sum_min_P = new_sum_P;
                    end
                else
                    if new_sum_P > sum_max_P
                        max_P = new_P;
                        sum_max_P = new_sum_P;
                    end
                end
                
            end
        end
    end
    
    if sum_min_P < sum_max_P
        A_min = min_P(1:na,1);
        B_min = min_P(na+1:na+nb,1);
        A_max = max_P(1:na,1);
        B_max = max_P(na+1:na+nb,1);
    else
       A_min = max_P(1:na,1);
       B_min = max_P(na+1:na+nb,1);
       A_max = min_P(1:na,1);
       B_max = min_P(na+1:na+nb,1);
    end
end