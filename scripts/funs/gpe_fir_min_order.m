function [rad_p,LSM,gpeP,P] = gpe_fir_min_order(u,y,model_error,max_r,id,int_est)
    % P = [Pmin Pmax]
    % LSM - parametre vypocitane metodou najmensich stvorcov
    % id = indikator, ci FIR ma byt v case t=0 rovny y=0 v ramci chyby
    % merania; ako ano id = 1.
    
    if size(u,1) < size(u,2)
        u = u';
    end
    
    if size(y,1) < size(y,2)
        y = y';
    end
    % odhad minimalneho radu modelu
    rad_p = 1;
    C = zeros(rad_p);
    v = 1;
    while v==1
        if id == 1
            C = napln_maticu(u,rad_p);
        else
            C = get_matica_vstupov(u,rad_p);
        end
        s = size(C);
        n = zeros(1,s(2));
        v = isempty(linprog(n,[-C;C],[model_error-y; model_error+y]));
        if v == 1
            rad_p = rad_p + 1;
        else
            rad_p;
        end
        if rad_p == max_r
            v = 0;
            error('Warning! Maximal order reached!')
        end
    end
    LSM = (C'*C)\C'*y;
    gpeP = linprog(n,[-C;C],[model_error-y; model_error+y]);
    if isnan(LSM)
        LSM = gpeP;
    end
    
    % intervalovy odhad parametrov
    if int_est == 1
        n = size(C);
        F = kron(eye(n(2)),[1;-1]);
        P = [];
        i = 1;

        for r = 1:1:(n(2))
            for s = 1:1:2
                temp_p = linprog(F(i,:),[-C;C],[model_error-y; model_error+y]);
                P(r,s) = temp_p(r);
                i=i+1;
            end
        end
    end
end