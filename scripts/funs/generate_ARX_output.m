function yk = generate_ARX_output(A,B,u,y0)
    
    if size(u,1) < size(u,2)
        u = u';
    end
    
    na = length(A);
    nb = length(B);
    y_mat = [y0 zeros([1 na-1])];
    u_mat = [u(1) zeros([1 nb-1])];
    for i = 2:1:length(u)
        y_mat(i,:) = [(u_mat(i-1,:)*B - y_mat(i-1,:)*A) y_mat(i-1,1:end-1)];
        u_mat(i,:) = [u(i) u_mat(i-1,1:end-1)];
    end
    yk = y_mat(:,1);
end