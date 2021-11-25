clear;
% close all;
clc;
format long
tic;

% Definition of parameters 
N = 20;%size
p = 2;

% parameters of memeory
mem_con = cell(p,1);
Jij = zeros(N,N);
coeff_save = zeros(1,N^2);
coeff_save2 = cell(1,N);
for k = 1:p
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

fname = ['mem_N',num2str(N),'_p',num2str(p),'_No1.mat'];
save(fname,'mem_con','coeff_save2','Jij','-v7.3');

toc;