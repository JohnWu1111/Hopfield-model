clear;
% close all;
clc;
format long
tic;

% Definition of parameters
N0 = 10:200;
num = 1e4;

stop = zeros(1,length(N0));

% parameters of memeory
parfor n = 1:length(N0)
    N = N0(n);
    
    p = 2;
    judge = 0;
    while judge == 0
        
        for k = 1:num
            mem_con = cell(p,1);
            coeff_save = zeros(1,N^2);
            coeff_save2 = cell(1,N);
            
            for i = 1:p
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
                    judge = 1;
                    break;
                end
                
                %         Etotal(i) = Ett;
                %     for j = 1:p
                %         temp = mem_con{j};
                %         spintotal(j,i) = mean(spin.*temp,'all')^2;
                %     end
            end
            
            if judge == 1
                break;
            end
        end        
        p = p + 1;
    end
    stop(n) = p-2;
end

figure;
plot(N0,stop);
% m2 = sqrt(mean(spintotal,2))


toc;