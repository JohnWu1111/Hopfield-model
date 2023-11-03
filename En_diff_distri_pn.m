clear;
% close all;
clc;
format long
tic;

N = 26; %size
p = 4;
J = 1;
num = 1;
hx = 1;
hy = 1;

q_save = cell(num,1);
e_save = cell(num,1);
s_save = cell(num,1);
count1 = 0;

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

for m = 1:num
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
    S_y = cell(len_p, 1);
    
    for i = 1:len_p
        S_p{i} = zeros(N_s(i) + 1);
        S_m{i} = zeros(N_s(i) + 1);
        
        for m = 1:N_s(i)
            S_p{i}(m, m + 1) = sqrt(S(i) * (S(i) + 1) - S_z{i}(m + 1) * (S_z{i}(m + 1) + 1));
            S_m{i}(m + 1, m) = sqrt(S(i) * (S(i) + 1) - S_z{i}(m) * (S_z{i}(m) - 1));
        end
        
        S_x{i} = (S_p{i} + S_m{i}) / 2;
        S_y{i} = (S_p{i} - S_m{i}) / (2*1i);
    end
    
    % construction of Hamiltonian
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
        temp1 = sparse(1);  
        temp2 = sparse(1);
        for k = 1:i - 1
            temp1 = kron(temp1, diag(sparse(ones(N_s(k) + 1, 1))));
            temp2 = kron(temp2, diag(sparse(ones(N_s(k) + 1, 1))));
        end
        temp1 = kron(temp1, S_x{i});
        temp2 = kron(temp2, S_y{i});   
        for k = i + 1:len_p
            temp1 = kron(temp1, diag(sparse(ones(N_s(k) + 1, 1))));
            temp2 = kron(temp2, diag(sparse(ones(N_s(k) + 1, 1))));
        end        
        H2 = H2 + hx * temp1;
%         H2 = H2 + hy * temp2;
    end
    
    H = H1 + H2;
    H = full(H);
    
    [V, D] = eig(H);
    e = diag(D);
    
    len = length(e);
    e = sort(e);
    s = zeros(len-1,1);
    r = zeros(len-2,1);
    for i = 1:len-1
        s(i) = e(i+1) - e(i);
    end
    
    for i = 1:len-2
        r(i) = s(i+1)/s(i);
    end
    
    q = min([r 1./r],[],2);
    
    count1 = count1 + 1;
    q_save{count1,1} = q;
    e_save{count1,1} = e;
    s_save{count1,1} = s;
    
    % s = sort(s);
end

count2 = 0;
for i = 1:count1
    len = length(q_save{i});
    qq(count2+1:count2+len,1) = q_save{i};
    count2 = count2 + len;
end

figure;
histogram(qq,200,'Normalization','pdf','DisplayStyle','stairs');
% axis([0 0.02 0 inf])

mean(qq)

toc;

function y = kron4(a,b,c,d)
y = kron(kron(kron(a,b),c),d);
end

function y = kron3(a,b,c)
y = (kron(kron(a,b),c));
end

function y = degen(A,J)
    y = 1;
    len = length(A);
    judge = A.*J.*A';
    judge = abs(judge);
    for i = 1:len
        for j = i+1:len
            a = judge(i,:);
            b = judge(j,:);
            a = sort(a);
            b = sort(b);
            if isequal(a,b)
                y = 0;
                return;
            end
        end
    end
end