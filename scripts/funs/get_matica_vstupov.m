function C = get_matica_vstupov(u,rad_p)
    a = zeros([1 rad_p]);
    C = zeros([length(u) rad_p]);

    for i = 1:1:length(u)
       a =  [u(i) a(1:end-1)];
       C(i,:) = a;
    end
end
