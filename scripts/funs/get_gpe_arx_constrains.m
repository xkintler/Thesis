function fk = get_gpe_arx_constrains(x,u,y0,na,nb)
    if isnumeric(x)
        A = zeros([length(u)-1 na]);
    else
        A = barvec([length(u)-1 na]);
    end
    for i = 1:1:length(u)-1
        if i <= na
            for j = 1:1:na
                if i == j
                    A(i,j) = -y0;
                elseif j < i
                    A(i,j) = -x(na+nb+i-j);
                else
                    A(i,j) = 0;
                end
            end
        else
            for j = 1:1:na
                A(i,j) = -x(na+nb+i-j);
            end
        end

        if(i<nb)
            for j = 1:1:i
                B(i,j) = u(i-j+1,1);
            end
        else
            for j = 1:1:nb
                B(i,j) = u(i-j+1,1);
            end
        end
    end
    C = (-1)*eye([length(u)-1 length(u)-1]);
    for i = 1:1:(length(u)-1+na+nb)
        theta(i,1) = x(i);
    end
    fk = [A B C]*theta;
end
