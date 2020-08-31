function [gpe_output,interval,info] = gpe_arx_min_order(u,y,model_error,y0,max_it,bound,int_est)
    % ohranicenie parametrov bound - treba upravit pri identifikacii
    % na urcenie parametrov metodou najmensich stvorcov treba pouzit ARX
    % toolbox
    
    na = 1;
    nb = 1;
    ef = 0;
    turn = 1;
    
    if size(u,1) < size(u,2)
        u = u';
    end
    
    if size(y,1) < size(y,2)
        y = y';
    end
    
    % minimal solution
    while ef~= 1
        %Objective
        fun = @(x) x(1);
        %Nonlinear Constraints
        nlcon = @(x) get_gpe_arx_constrains(x,u,y0,na,nb);

        cl = zeros([length(u)-1 1]);
        cu = zeros([length(u)-1 1]);
        %Bounds
        lb = [(-bound)*ones([na+nb 1]); y(2:end)-model_error];
        ub = [bound*ones([na+nb 1]); y(2:end)+model_error];

        %Solve
        [x,~,ef,info] = baron(fun,[],[],[],lb,ub,nlcon,cl,cu);

        if ef~=1
            if na > 2
                nb = nb + 1;
                na = 1;
            else
                na = na + 1;
            end
            if turn == max_it
                error('Warning! Maximal number of iterations reached!')
            end
            turn = turn + 1;
        end
    end
    
    %data conditioning
    gpeA = x(1:na);
    gpeB = x(na+1:na+nb);
    arx_out = [y0; x(na+nb+1:end)];
    gpe_output = struct('na',na,'nb',nb,'A',gpeA,'B',gpeB,'ARX_output',arx_out);
    
    % parameters interval estimation
    if int_est == 1
        P = zeros([na+nb 2]);

        for s = 1:1:na+nb
            for j = 1:1:2
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
                [x,~,ef] = baron(fun,[],[],[],lb,ub,nlcon,cl,cu);
                if ef ~= 1
                    error('Warning! Problem is infeasible!')
                end

                P(s,j) = x(s);
            end
        end
        A = P(1:na,1:2);
        B = P(na+1:na+nb,1:2);
        interval = struct('A',A,'B',B);
    end
end