clear;
% close all;
clc;
format long
tic;

% Definition of parameters
n = 4; %size
N = n^2;

% parameters of memeory
p = 10;
mem_con = cell(p,1);
coeff_save = zeros(N,N);
coeff_save2 = cell(n,n);
for i = 1:p
    mem_con{i} = round(rand(n))*2-1;
    coeff_save = coeff_save + kron(mem_con{i},mem_con{i});
end

for i = 1:n
    for j = 1:n
        coeff_save2{i,j} = coeff_save((i-1)*n+1:i*n,(j-1)*n+1:j*n);
    end
end

fname = ['mem_N',num2str(N),'_p',num2str(p),'_No1.mat'];
save(fname,'mem_con','coeff_save2','-v7.3');

toc;