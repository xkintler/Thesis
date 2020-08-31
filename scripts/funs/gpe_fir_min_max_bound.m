function [min_P, max_P] = gpe_fir_min_max_bound(r_min,u,y,model_error,id)
    % id = indikator, ci FIR ma byt v case t=0 rovny y=0 v ramci chyby
    % merania; ako ano id = 1 inak id = 0.
    
    if size(u,1) < size(u,2)
        u = u';
    end
    
    if size(y,1) < size(y,2)
        y = y';
    end
    
    if id == 1
        C = napln_maticu(u,r_min);
    else
        C = get_matica_vstupov(u,r_min);
    end
    F = eye(r_min);

    for r = 1:1:r_min
        
        if r == 1
           min_P = linprog(F(r,:),[-C;C],[model_error-y; model_error+y]);
           sum_min_P = sum(min_P);
           
           max_P = linprog(-F(r,:),[-C;C],[model_error-y; model_error+y]);
           sum_max_P = sum(max_P);
        else
            % minimal bound P
            new_P = linprog(F(r,:),[-C;C],[model_error-y; model_error+y]);
            if sum(new_P) < sum_min_P
                min_P = new_P;
                sum_min_P = sum(min_P);
            end
            % maximal bound P
            new_P = linprog(-F(r,:),[-C;C],[model_error-y; model_error+y]);
            if sum(new_P) > sum_max_P
                max_P = new_P;
                sum_max_P = sum(max_P);
            end
        end  
    end    
end