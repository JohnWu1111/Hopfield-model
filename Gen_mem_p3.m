clear;
% close all;
clc;
format long
tic;

% Definition of parameters
N = 50;%size
p = 3;
num = 1;


for n = 1:num
    % parameters of memeory
    mem_con = cell(p,1);
    Jij = zeros(N,N);
    coeff_save = zeros(1,N^2);
    coeff_save2 = cell(1,N);
    mem_con{1} = ones(1,N);
    for k = 2:p
        mem_con{k} = round(rand(1,N))*2-1;
        coeff_save = coeff_save + kron(mem_con{k},mem_con{k});
        temp = mem_con{k};
        for i = 1:N
            for j = i+1:N
                Jij(i,j) = Jij(i,j) + temp(i)*temp(j);
            end
        end
    end
    Jij = Jij/N;
    
    for i = 1:N
        coeff_save2{i} = coeff_save((i-1)*N+1:i*N);
    end
    
    overlap12 = mean(mem_con{1}.*mem_con{2})
    overlap13 = mean(mem_con{1}.*mem_con{3})
    overlap23 = mean(mem_con{2}.*mem_con{3})
%     if overlap12 == 0 && overlap13 == 0 && overlap23 == 0
%         break;
%     end
end

fname = ['mem_N',num2str(N),'_p',num2str(p),'_No1.mat'];
save(fname,'mem_con','coeff_save2','Jij','-v7.3');

toc;