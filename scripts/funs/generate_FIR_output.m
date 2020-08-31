function y = generate_FIR_output(theta,u)

    if size(u,1) < size(u,2)
        u = u';
    end
    
    r = length(theta);
    C = napln_maticu(u,r);
    
    y = C*theta;
end