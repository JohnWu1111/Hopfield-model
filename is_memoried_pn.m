clear;
% close all;
clc;
format long
tic;

% Definition of parameters
N = 250; %size
p = 2;
dt = 10;
T = 1000;
t = 0:dt:T;
nt = length(t);
J = 1;
h = 0.1;
No = 3;
cut = 100;

%rng(3e5 * No);
% rng(30)

% load(strcat('mem_N',num2str(N),'_p',num2str(p),'_No',num2str(No),'.mat'));
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

% generation of memeory
N_s = zeros(len_p, 1);
N_pos = zeros(N,1);
while degen(N_s, Jij) == 0
    S_mem = zeros(p,len_p);
    N_s = zeros(len_p, 1);
    mem_con = zeros(p, N);
    mem_con(1, :) = ones(1, N);
    mem_con(2:p, :) = round(rand(p - 1, N)) * 2 - 1;
    for i = 1:N       
        for j = 1:len_p
            temp = mean(mem_con(:, i) .* N_div(:, j));
            
            if temp == 1
                N_s(j) = N_s(j) + 1;
                N_pos(i) = j;
                S_mem(:,j) = S_mem(:,j) + mem_con(:,i);
                break
            end
            
        end       
    end
    
end

NN = sum(N_s);

NN = 1;
for i = 1:len_p
    NN = NN*(N_s(i)+1);
end

mem_rank = zeros(p,1);
for i = 1:p
    temp = sparse(1);
    for j = 1:len_p
        temp1 = zeros(N_s(j)+1,1);
        pos = (N_s(j)-S_mem(i,j))/2+1;
        temp1(pos) = 1;
        temp = kron(temp,temp1);
    end
    [mem_rank(i),~] = find(temp);
%     if mem_rank(i) > NN/2
%         mem_rank(i) = NN - mem_rank(i) +1;
%     end
end

S = N_s/2;

% construction of matrice
S_z = cell(len_p,1);
for i = 1:len_p
    S_z{i} = zeros(N_s(i)+1,1);
    for m = 1:N_s(i)+1
        S_z{i}(m) = N_s(i)/2 - (m-1);
    end
end

S_p = cell(len_p,1);
S_m = cell(len_p,1);
S_x = cell(len_p,1);
for i = 1:len_p
    S_p{i} = zeros(N_s(i)+1);
    S_m{i} = zeros(N_s(i)+1);
    for m = 1:N_s(i)
        S_p{i}(m,m+1) = sqrt(S(i)*(S(i)+1)-S_z{i}(m+1)*(S_z{i}(m+1)+1));
        S_m{i}(m+1,m) = sqrt(S(i)*(S(i)+1)-S_z{i}(m)*(S_z{i}(m)-1));
    end
    
    S_p{i} = sparse(S_p{i});
    S_m{i} = sparse(S_m{i});
    
    S_x{i} = (S_p{i} + S_m{i})/2;
end

% construction of Hamiltonian
H1 = zeros(NN,1);
for i = 1:len_p
    temp = 1;
    for k = 1:i-1
        temp = kron(temp,ones(N_s(k)+1,1));
    end
    temp = kron(temp,S_z{i}.^2);
    for k = i+1:len_p
        temp = kron(temp,ones(N_s(k)+1,1));
    end
    H1 = H1 + 4*temp;
    
    for j = i+1:len_p
        if Jij(i,j) == 0
            continue
        end
        temp = ones(1);
        for k = 1:i-1
            temp = kron(temp,ones(N_s(k)+1,1));
        end
        temp = kron(temp,S_z{i});
        for k = i+1:j-1
            temp = kron(temp,ones(N_s(k)+1,1));
        end
        temp = kron(temp,S_z{j});
        for k = j+1:len_p
            temp = kron(temp,ones(N_s(k)+1,1));
        end
        H1 = H1 + Jij(i,j)*temp;
    end
end
H1 = -J*(H1)/(2*N);
En = H1;
H1 = sparse(H1);
H1 = diag(H1);

H2 = sparse(NN,NN);
for i = 1:len_p
    temp = sparse(1);
    for k = 1:i-1
        temp = kron(temp,diag(sparse(ones(N_s(k)+1,1))));
    end
    temp = kron(temp,S_x{i});
    for k = i+1:len_p
        temp = kron(temp,diag(sparse(ones(N_s(k)+1,1))));
    end
    H2 = H2 + h*temp;
end
H = H1 + H2;
H = full(H);

clear H1 H2;

fprintf("diagonalization started!\n No = %d, N = %d\n", No, NN)

[V, D] = eig(H);
e = diag(D);

clear H D;

% VV = V.^2;
% clear V
V_inv = inv(V);

[En_sort, rank] = sort(En);
[En_re, rank_ori, rank_re] = myreduce(En_sort,rank);
NE = length(En_re);

for i = 1:p
    temp = find(rank == mem_rank(i));
    for j = 1:NE-1
        if rank_re(j) <= temp && rank_re(j+1) > temp
            mem_rank(i) = rank_re(j);
            break
        end
    end
end

% target = [rank_re(1:10)
%     rank_re(floor(NE/3)-4:floor(NE/3)+5)
%     rank_re(floor(NE/2)-4:floor(NE/2)+5)
%     rank_re(end-9:end)];
% target_E = [En_re(1:10)
%     En_re(floor(NE/3)-4:floor(NE/3)+5)
%     En_re(floor(NE/2)-4:floor(NE/2)+5)
%     En_re(end-9:end)];

% if NE < 30
%     warning('En_re < 30')
% end

pro = zeros(NN,nt);
% pro(:,1) = 1;

eiet = exp(-1i*e*t);
V = V(rank,:);
V_inv = V_inv(:,rank);
trans = V'.*V_inv;

for i = 1:NN
%     trans = V(rank(i),:)'.*V_inv(:,rank(i));
%     temp = V*(exp(-1i*e*t).*trans);
%     pro(i,:) = temp(target(i),:);
    pro(i,:) = sum(trans(:,i).*eiet);
end

pro = abs(pro).^2;

pro1 = zeros(NE,nt);

for i = 1:NE-1
    pro1(i,:) = mean(pro(rank_re(i):rank_re(i+1),:));
end
pro1(NE,:) = mean(pro(rank_re(NE):end,:));
mean_pro1 = mean(pro1(:,floor(nt/2):end),2);
e_re = (En_re - En_re(1))/(En_re(end) - En_re(1));

pro1_mean = mean(pro1(:,floor(nt/2):end),2);

figure;
plot(e_re,mean_pro1)

toc;

function y = kron_p(a,b)
la = length(a);
lb = length(b);
y = zeros(la*lb,1);
for i = 1:la
    for j = 1:lb
        y((i-1)*lb+j) = a(i) + b(j);
    end
end
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

function [out, r_ori, r_re] = myreduce(A,rank)
    len = length(A);
    temp = A(1);
    out(1,1) = A(1);
    r_ori(1,1) = rank(1);
    r_re(1,1) = 1;
    count = 2;
    for i = 2:len
        if A(i) ~= temp
            out(count,1) = A(i);
            r_ori(count,1) = rank(i);
            r_re(count,1) = i;
            temp = A(i);
            count = count + 1;
        end
    end
end