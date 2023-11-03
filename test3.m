% calculation of some thermaldynamical observables

clear;
% close all;
format long
tic;

% Definition of parameters
p = 4;
J = 1;
%num = 10;

N = 16;
No = 2;
h = 2;

rng(3e5 * No);

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

% N_s = [16 0 0 0 0 0 0 0];

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

    H2 = H2 + h * temp;
end

H = H1 + H2;
H = full(H);

clear H1 H2;

fprintf("diagonalization started!\n No = %d, N = %d\n", No, NN)

[V, D] = eig(H);
e = diag(D);

clear H D;

VV = V.^2;
clear V

% % order parameter
% M2 = cell(p, 1);
% 
% for i = 1:p
%     temp = 1;
% 
%     for j = 1:len_p
%         temp = kron_p(temp, N_div(i, j) * S_z{j});
%     end
% 
%     temp = temp.^2;
%     M2{i} = temp;
% end

% mm = zeros(p, NN);
% 
% for i = 1:p
%     mm(i, :) = sqrt(M2{i}' * VV) / N;
% end

% SG chi
corre_in = zeros(1, NN);
corre_ext = zeros(1, NN);
weight1 = 0;
weight2 = 0;
core_in_store = zeros(len_p,NN);
core_ext_store = zeros(len_p*(len_p-1)/2,NN);
count = 0;

for i = 1:len_p
    if N_s(i) == 0
        count = count + len_p - i;
        continue
    end
    temp = 1;
    for k = 1:i - 1
        temp = kron(temp, ones(N_s(k) + 1, 1));
    end
    temp = kron(temp, 4 * S_z{i}.^2);
    for k = i + 1:len_p
        temp = kron(temp, ones(N_s(k) + 1, 1));
    end
    temp1 = temp' * VV;

    if N_s(i) > 1
        temp2 = temp1 - N_s(i);
        temp1 = temp2;
        weight1 = weight1 + 2 * nchoosek(N_s(i), 2);
        %/ (2 * nchoosek(N_s(i), 2))
    else
        temp1 = 0;
    end
    temp1 = abs(temp1);

    corre_in = corre_in + temp1;
    core_in_store(i,:) = temp1/N_s(i);

    for j = i + 1:len_p
        count = count +1;
        if N_s(j) == 0
            continue
        end
        temp = 1;
        for k = 1:i - 1
            temp = kron(temp, ones(N_s(k) + 1, 1));
        end
        temp = kron(temp, 2 * S_z{i});
        for k = i + 1:j - 1
            temp = kron(temp, ones(N_s(k) + 1, 1));
        end
        temp = kron(temp, 2 * S_z{j});
        for k = j + 1:len_p
            temp = kron(temp, ones(N_s(k) + 1, 1));
        end
        temp1 = abs(temp' * VV);
        corre_ext = corre_ext + temp1;
        weight2 = weight2 + N_s(i) * N_s(j);
        
        core_ext_store(count,:) = temp1/(N_s(i) * N_s(j));
    end

end

corre_in = corre_in/weight1;
corre_ext = corre_ext / weight2;

% % phi4
% phi4 = sum(VV.^2);
% 
% % corre of the largest two spin
% [B, I] = sort(N_s);
% a = min([I(end), I(end - 1)]);
% b = max([I(end), I(end - 1)]);
% 
% temp = 1;
% 
% for i = 1:a - 1
%     temp = kron(temp, ones(N_s(i) + 1, 1));
% end
% 
% temp = kron(temp, 2 * S_z{a} ./ N_s(a));
% 
% for i = a + 1:b - 1
%     temp = kron(temp, ones(N_s(i) + 1, 1));
% end
% 
% temp = kron(temp, 2 * S_z{b} ./ N_s(b));
% 
% for i = b + 1:len_p
%     temp = kron(temp, ones(N_s(i) + 1, 1));
% end
% 
% corre = temp' * VV;
% 
% % energy distribution
% len = length(e);
% e = sort(e);
% s = zeros(len - 1, 1);
% r = zeros(len - 2, 1);
% 
% for i = 1:len - 1
%     s(i) = e(i + 1) - e(i);
% end
% 
% for i = 1:len - 2
%     r(i) = s(i + 1) / s(i);
% end
% 
% q = min([r 1 ./ r], [], 2);

%save(strcat('main_N', num2str(N), '_p', num2str(p), '_h', num2str(h), '_No', num2str(No), '.mat'), 'N_s', 'sus_save', 'mm_save', 'phi4_save', 'corre_save', 'h', 'e_save', 'r_save', 's_save', 'q_save', '-v7.3');
% save(strcat('main_N', num2str(N), '_p', num2str(p), '_h', num2str(h), '_No', num2str(No), '.mat'), 'N_s', 'NN', 'sus', 'mm', 'phi4', 'corre', 'h', 'e', 'r', 's', 'q', '-v7.3');
% 
figure;
plot(1:NN, corre_in);
xlabel('\epsilon')
ylabel('internal correlation')
% saveas(gcf, strcat('order_N', num2str(N), '_p', num2str(p), '_No', num2str(No), '.fig'));
% print(strcat('order_N', num2str(N), '_p', num2str(p), '_No', num2str(No)), '-dpng', '-r0');

figure;
plot(1:NN, corre_ext);
xlabel('\epsilon')
ylabel('external correlation')

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