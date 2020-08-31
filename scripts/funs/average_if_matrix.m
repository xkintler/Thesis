function [avg_m,ele_count] = average_if_matrix(A)
size_m = size(A);
init_n = unique(A(:,1))';

avg_m = zeros(length(init_n),size_m(2)-1);
ele_count = zeros([length(init_n) 1]);
for i = 1:length(init_n)
    elmt_count = 0;
    for j = 1:size_m(1)
        if A(j,1) == init_n(i)
            elmt_count = elmt_count + 1;
            avg_m(i,:) = avg_m(i,:) + A(j,2:end);
        end
    end
    avg_m(i,:) = avg_m(i,:)./elmt_count;
    ele_count(i,1) = elmt_count;
end
avg_m = [init_n' avg_m];
ele_count = [init_n' ele_count];
end