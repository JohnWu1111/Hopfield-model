% calculation of some thermaldynamical observables

clear;
% close all;
format long
tic;

% Definition of parameters
p = 4;
J = 1;
num = 100;

N_a = 18:2:26;
lN = length(N_a);

h = 0.1:0.1:1;
lh = length(h);

% sus_all = zeros(lN,lh);
mm = zeros(lN,lh);

% rng(3e5 * No);

% initialization
ol = zeros(p, p);
len_p = 2^(p - 1);
N_div = zeros(p, len_p);

for i = p:-1:1
    a = 2^(i - 1);
    
    for j = 1:len_p
        temp1 = floor((j - 1) / a);
        temp2 = mod(temp1, 2);
        N_div(p - i + 1, j) = -2 * temp2 + 1;
    end
    
end

Jij = zeros(len_p, len_p);

for i = 1:len_p
    
    for j = i + 1:len_p
        Jij(i, j) = 2 * sum(N_div(:, i) .* N_div(:, j));
        Jij(j, i) = Jij(i, j);
    end
    
end

Jij = Jij + eye(len_p) * p;

for t = 1:lN
    N = N_a(t);
    for No = 1:num
        
        % generation of memeory
        N_s = zeros(len_p, 1);
        
        while degen(N_s, Jij) == 0
            N_s = zeros(len_p, 1);
            mem_con = zeros(p, N);
            mem_con(1, :) = ones(1, N);
            mem_con(2:p, :) = round(rand(p - 1, N)) * 2 - 1;
            
            for i = 1:N
                
                for j = 1:len_p
                    temp = mean(mem_con(:, i) .* N_div(:, j));
                    
                    if temp == 1
                        N_s(j) = N_s(j) + 1;
                        break
                    end
                    
                end
                
            end
            
        end
        
        NN = sum(N_s);
        
        NN = 1;
        
        for i = 1:len_p
            NN = NN * (N_s(i) + 1);
        end
        
        S = N_s / 2;
        
        % construction of matrice
        S_z = cell(len_p, 1);
        
        for i = 1:len_p
            S_z{i} = zeros(N_s(i) + 1, 1);
            
            for m = 1:N_s(i) + 1
                S_z{i}(m) = N_s(i) / 2 - (m - 1);
            end
            
        end
        
        S_p = cell(len_p, 1);
        S_m = cell(len_p, 1);
        S_x = cell(len_p, 1);
        
        for i = 1:len_p
            S_p{i} = zeros(N_s(i) + 1);
            S_m{i} = zeros(N_s(i) + 1);
            
            for m = 1:N_s(i)
                S_p{i}(m, m + 1) = sqrt(S(i) * (S(i) + 1) - S_z{i}(m + 1) * (S_z{i}(m + 1) + 1));
                S_m{i}(m + 1, m) = sqrt(S(i) * (S(i) + 1) - S_z{i}(m) * (S_z{i}(m) - 1));
            end
            
            S_x{i} = (S_p{i} + S_m{i}) / 2;
        end
        
        % construction of Hamiltonian
        
        for n = 1:lh
            H1 = zeros(NN, 1);
            
            for i = 1:len_p
                temp = 1;
                
                for k = 1:i - 1
                    temp = kron(temp, ones(N_s(k) + 1, 1));
                end
                
                temp = kron(temp, S_z{i}.^2);
                
                for k = i + 1:len_p
                    temp = kron(temp, ones(N_s(k) + 1, 1));
                end
                
                H1 = H1 + 4 * temp;
                
                for j = i + 1:len_p
                    
                    if Jij(i, j) == 0
                        continue
                    end
                    
                    temp = ones(1);
                    
                    for k = 1:i - 1
                        temp = kron(temp, ones(N_s(k) + 1, 1));
                    end
                    
                    temp = kron(temp, S_z{i});
                    
                    for k = i + 1:j - 1
                        temp = kron(temp, ones(N_s(k) + 1, 1));
                    end
                    
                    temp = kron(temp, S_z{j});
                    
                    for k = j + 1:len_p
                        temp = kron(temp, ones(N_s(k) + 1, 1));
                    end
                    
                    H1 = H1 + Jij(i, j) * temp;
                end
                
            end
            
            H1 = -J * (H1) / (2 * N);
            H1 = sparse(H1);
            H1 = diag(H1);
            
            H2 = sparse(NN, NN);
            
            for i = 1:len_p
                temp = sparse(1);
                
                for k = 1:i - 1
                    temp = kron(temp, diag(sparse(ones(N_s(k) + 1, 1))));
                end
                
                temp = kron(temp, S_x{i});
                
                for k = i + 1:len_p
                    temp = kron(temp, diag(sparse(ones(N_s(k) + 1, 1))));
                end
                
                H2 = H2 + temp;
            end
            
            H = H1 + h(n) * H2;
            
            clear H1 H2;
            
            [V,D] = eigs(H,1,'smallestreal');
            e = diag(D);
            
            clear H D;
            VV = V.^2;
            clear V
            
            % energy distribution
%             sus = 0;
%             
%             for i = 1:len_p
%                 
%                 if N_s(i) == 0
%                     continue
%                 end
%                 
%                 temp = 1;
%                 
%                 for k = 1:i - 1
%                     temp = kron(temp, ones(N_s(k) + 1, 1));
%                 end
%                 
%                 temp = kron(temp, 4 * S_z{i}.^2);
%                 
%                 for k = i + 1:len_p
%                     temp = kron(temp, ones(N_s(k) + 1, 1));
%                 end
%                 
%                 temp1 = temp' * VV;
%                 
%                 if N_s(i) > 1
%                     temp2 = temp1 - N_s(i);
%                     temp1 = N_s(i) + temp2.^2 / (2 * nchoosek(N_s(i), 2));
%                 else
%                     temp1 = 1;
%                 end
%                 
%                 sus = sus + temp1;
%                 %         sus(n) = sus(n) + mean(temp.^2);
%                 
%                 for j = i + 1:len_p
%                     
%                     if N_s(j) == 0
%                         continue
%                     end
%                     
%                     temp = 1;
%                     
%                     for k = 1:i - 1
%                         temp = kron(temp, ones(N_s(k) + 1, 1));
%                     end
%                     
%                     temp = kron(temp, 2 * S_z{i});
%                     
%                     for k = i + 1:j - 1
%                         temp = kron(temp, ones(N_s(k) + 1, 1));
%                     end
%                     
%                     temp = kron(temp, 2 * S_z{j});
%                     
%                     for k = j + 1:len_p
%                         temp = kron(temp, ones(N_s(k) + 1, 1));
%                     end
%                     
%                     temp1 = temp' * VV;
%                     %             sus(n) = sus(n) + mean(temp.^2);
%                     sus = sus + 2 * temp1.^2 / (N_s(i) * N_s(j));
%                 end
%                 
%             end
%             
%             sus_all(t,n) = sus_all(t,n) + sus / N;
            M2 = cell(p, 1);

            for i = 1:p
                temp = 1;
                
                for j = 1:len_p
                    temp = kron_p(temp, N_div(i, j) * S_z{j});
                end
                
                temp = temp.^2;
                M2{i} = temp;
            end
            mmm = zeros(p,1);
            for i = 1:p
                mmm(i) = sqrt(M2{i}' * VV) / N;
            end
            mm(t,n) = mm(t,n) + max(mmm);
        end
    end
end

% sus_all = sus_all/num;
mm = mm/num;

le = cell(1, lN);

for i = 1:lN
    le{i} = strcat('N', num2str(N_a(i)));
end

figure;
plot(h, mm);
xlabel('h')
ylabel('order parameter')
legend(le);

toc;

function y = kron4(a, b, c, d)
y = kron(kron(kron(a, b), c), d);
end

function y = kron3(a, b, c)
y = (kron(kron(a, b), c));
end

function y = kron_p(a, b)
la = length(a);
lb = length(b);
y = zeros(la * lb, 1);

for i = 1:la
    
    for j = 1:lb
        y((i - 1) * lb + j) = a(i) + b(j);
    end
    
end

end

function y = kron_p4(a, b, c, d)
y = kron_p(kron_p(kron_p(a, b), c), d);
end

function y = degen(A, J)
y = 1;
len = length(A);
judge = A .* J .* A';
judge = abs(judge);

for i = 1:len
    
    for j = i + 1:len
        a = judge(i, :);
        b = judge(j, :);
        a = sort(a);
        b = sort(b);
        
        if isequal(a, b)
            y = 0;
            return;
        end
        
    end
    
end

end
