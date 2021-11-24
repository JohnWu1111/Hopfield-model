clear;
% close all;
clc;
format long
tic;

% Definition of parameters
N = 1000;
num = 1e4;

p = 2:100;

count = zeros(num,length(p));

% parameters of memeory
parfor m = 1:length(p)
       
    for k = 1:num
        mem_con = cell(p(m),1);
        coeff_save = zeros(1,N^2);
        coeff_save2 = cell(1,N);
    
        for i = 1:p(m)
            mem_con{i} = round(rand(1,N))*2-1;
            coeff_save = coeff_save + kron(mem_con{i},mem_con{i});
        end
        
        for i = 1:N
            coeff_save2{i} = coeff_save((i-1)*N+1:i*N);
        end
        
        spin = mem_con{1};    
%         spintotal = zeros(p,N);
        
        
        for i = 1:N
            coeff_it = coeff_save2{i};
            dEt = (mean(coeff_it.*spin,'all')-coeff_it(i)/N)*spin(i);
            
            if dEt < 0
                count(k,m) = count(k,m) +1;
            end
            
            %         Etotal(i) = Ett;
            %     for j = 1:p
            %         temp = mem_con{j};
            %         spintotal(j,i) = mean(spin.*temp,'all')^2;
            %     end
        end
    end    
end
per = mean(count);
per_std = std(count);

figure;
errorbar(p,per,per_std);
% m2 = sqrt(mean(spintotal,2))


toc;