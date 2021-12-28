clear;
% close all;
clc;
format long
tic;

% Definition of parameters
N = 100; %size
p = 3;
num = 1e6;

ind = zeros(num,1);
ol = zeros(num,1);

for n = 1:num
    mem_con = cell(p,1);
    coeff_save = zeros(1,N^2);
    coeff_save2 = cell(1,N);
    mem_con{1} = ones(1,N);
    for k = 2:p
        mem_con{k} = round(rand(1,N))*2-1;
    end
    
    N1 = 0;
    N2 = 0;
    N3 = 0;
    N4 = 0;
    for i = 1:N
        if mem_con{2}(i)*mem_con{3}(i) > 0
            if mem_con{2}(i) == 1
                N1 = N1 + 1;
            else
                N4 = N4 + 1;
            end
        else
            if mem_con{2}(i) == 1
                N2 = N2 + 1;
            else
                N3 = N3 + 1;
            end
        end
    end
    
    M = [N1 N2 N3 N4];
    M = sort(M);
    N_diff = diff(M);
    ind(n) = sum(N_diff);
    
    ol12 = mean(mem_con{1}.*mem_con{2});
    ol13 = mean(mem_con{1}.*mem_con{3});
    ol23 = mean(mem_con{2}.*mem_con{3});
    ol(n) = mean(abs([ol12 ol13 ol23]));
end

figure;
histogram(ind,100,'DisplayStyle','stairs');

figure;
histogram(ol,100,'DisplayStyle','stairs');

toc;