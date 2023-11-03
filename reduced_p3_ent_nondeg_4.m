% calculation of some thermaldynamical observables

clear;
% close all;
format long
tic;

% Definition of parameters
p = 3;
J = 1;

N = 20;
No = 1;
h = 0.1;

rng(1e4 * No);

N1 = 0;
N2 = 0;
N3 = 0;
N4 = 0;

while isdeg(N1, N2, N3, N4)

    N1 = 0;
    N2 = 0;
    N3 = 0;
    N4 = 0;

    mem_con = cell(p, 1);
    mem_con{1} = ones(1, N);

    for k = 2:p
        mem_con{k} = round(rand(1, N)) * 2 - 1;
    end

    for i = 1:N

        if mem_con{2}(i) * mem_con{3}(i) > 0

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

end

temp = sort([N1,N2,N3,N4],'descend');
N1 = temp(1);
N2 = temp(2);
N3 = temp(3);
N4 = temp(4);

% N1 = 1;
% N2 = 1;
% N3 = 1;
% N4 = 0;

NN = (N1 + 1) * (N2 + 1) * (N3 + 1) * (N4 + 1);

S1 = N1 / 2;
S2 = N2 / 2;
S3 = N3 / 2;
S4 = N4 / 2;

S1_z = zeros(N1 + 1, 1);
S1_p = zeros(N1 + 1);
S1_m = zeros(N1 + 1);
S2_z = zeros(N2 + 1, 1);
S2_p = zeros(N2 + 1);
S2_m = zeros(N2 + 1);
S3_z = zeros(N3 + 1, 1);
S3_p = zeros(N3 + 1);
S3_m = zeros(N3 + 1);
S4_z = zeros(N4 + 1, 1);
S4_p = zeros(N4 + 1);
S4_m = zeros(N4 + 1);

% construction of matrice
for m = 1:N1 + 1
    S1_z(m) = N1 / 2 - (m - 1);
end

for m = 1:N2 + 1
    S2_z(m) = N2 / 2 - (m - 1);
end

for m = 1:N3 + 1
    S3_z(m) = N3 / 2 - (m - 1);
end

for m = 1:N4 + 1
    S4_z(m) = N4 / 2 - (m - 1);
end

for m = 1:N1
    S1_p(m, m + 1) = sqrt(S1 * (S1 + 1) - S1_z(m + 1) * (S1_z(m + 1) + 1));
    S1_m(m + 1, m) = sqrt(S1 * (S1 + 1) - S1_z(m) * (S1_z(m) - 1));
end

for m = 1:N2
    S2_p(m, m + 1) = sqrt(S2 * (S2 + 1) - S2_z(m + 1) * (S2_z(m + 1) + 1));
    S2_m(m + 1, m) = sqrt(S2 * (S2 + 1) - S2_z(m) * (S2_z(m) - 1));
end

for m = 1:N3
    S3_p(m, m + 1) = sqrt(S3 * (S3 + 1) - S3_z(m + 1) * (S3_z(m + 1) + 1));
    S3_m(m + 1, m) = sqrt(S3 * (S3 + 1) - S3_z(m) * (S3_z(m) - 1));
end

for m = 1:N4
    S4_p(m, m + 1) = sqrt(S4 * (S4 + 1) - S4_z(m + 1) * (S4_z(m + 1) + 1));
    S4_m(m + 1, m) = sqrt(S4 * (S4 + 1) - S4_z(m) * (S4_z(m) - 1));
end

S1_x = (S1_p + S1_m) / 2;
S2_x = (S2_p + S2_m) / 2;
S3_x = (S3_p + S3_m) / 2;
S4_x = (S4_p + S4_m) / 2;

len_p = 4;
N_s = [N1 N2 N3 N4];
S_z = {S1_z, S2_z, S3_z, S4_z};
% Jij = [3 2 2 -2;
%        0 3 -2 2;
%        0 0 3 2;
%        0 0 0 3];
Jij = [3 2 2 2;
    0 3 2 2;
    0 0 3 2;
    0 0 0 3];

% construction of Hamiltonian
H1 = -J * (3 * kron(S1_z.^2, ones((N2 + 1) * (N3 + 1) * (N4 + 1), 1)) ...
    + 3 * kron3(ones(N1 + 1, 1), S2_z.^2, ones((N3 + 1) * (N4 + 1), 1)) ...
    +3 * kron3(ones((N1 + 1) * (N2 + 1), 1), S3_z.^2, ones(N4 + 1, 1)) ...
    +3 * kron(ones((N1 + 1) * (N2 + 1) * (N3 + 1), 1), S4_z.^2) ... %inner term
    +2 * kron3(S1_z, kron(S2_z, ones(N3 + 1, 1)) + kron(ones(N2 + 1, 1), S3_z), ones(N4 + 1, 1)) ...
    +2 * kron3(ones(N1 + 1, 1), kron(S2_z, ones(N3 + 1, 1)) + kron(ones(N2 + 1, 1), S3_z), S4_z) ...
    -2 * kron4(ones(N1 + 1, 1), S2_z, S3_z, ones(N4 + 1, 1)) ...
    -2 * kron3(S1_z, ones((N2 + 1) * (N3 + 1), 1), S4_z)) / (2 * N);
H1 = diag(H1);
H2 = h * (kron(S1_x, eye((N2 + 1) * (N3 + 1) * (N4 + 1))) ...
    +kron3(eye(N1 + 1), S2_x, eye((N3 + 1) * (N4 + 1))) ...
    +kron3(eye((N1 + 1) * (N2 + 1)), S3_x, eye(N4 + 1)) ...
    +kron(eye((N1 + 1) * (N2 + 1) * (N3 + 1)), S4_x));
H = H1 + H2;

clear H1 H2;

[V, D] = eig(H);
e = diag(D);

clear H D;

% spin #1
num_r = (N2 + 1) * (N3 + 1) * (N4 + 1);
E1= zeros((N1+1),num_r);
ent1 = zeros(NN,1);
for n = 1:NN
        for i = 1:N1+1
            for j = 1:num_r
                pos = (i-1)*num_r+j;
                E1(i,j) = V(pos,n);
            end
        end
        R1 = zeros(N1+1);
        for j = 1:num_r
            temp = E1(:,j);
            temp1 = temp*temp';
            R1 = R1 + temp1;
        end
    
    [RV,RD] = eig(R1);
    for k = 1:length(RD)
        if abs(RD(k,k))<1e-14
            RD(k,k) = 1;
        end
    end
    RE = diag(RD);
    temp3 = - trace(RD*diag(log(RE)));
    ent1(n) = real(temp3);
end

% spin #2
num_r = (N1 + 1) * (N3 + 1) * (N4 + 1);
E2= zeros((N2+1),num_r);
ent2 = zeros(NN,1);
for n = 1:NN
        for i = 1:N1+1
            for j = 1:N2+1
                for k = 1:(N3 + 1) * (N4 + 1)
                   pos = (j-1)*(N3+1)*(N4+1)+k;
                   pos = (i-1)*(N2+1)+pos;
                   E2(j,(i-1)*(N3 + 1)*(N4 + 1)+k) = V(pos,n);
                end
            end
        end
        R2 = zeros(N2+1);
        for j = 1:num_r
            temp = E2(:,j);
            temp1 = temp*temp';
            R2 = R2 + temp1;
        end
    
    [RV,RD] = eig(R2);
    for k = 1:length(RD)
        if abs(RD(k,k))<1e-14
            RD(k,k) = 1;
        end
    end
    RE = diag(RD);
    temp3 = - trace(RD*diag(log(RE)));
    ent2(n) = real(temp3);
end

% spin #3
num_r = (N1 + 1) * (N2 + 1) * (N4 + 1);
E3= zeros((N3+1),num_r);
ent3 = zeros(NN,1);
for n = 1:NN
        for i = 1:(N1 + 1) * (N2 + 1)
            for j = 1:N3+1
                for k = 1:N4+1
                   pos = (j-1)*(N4+1)+k;
                   pos = (i-1)*(N3+1)+pos;
                   E3(j,(i-1)*(N4+1)+k) = V(pos,n);
                end
            end
        end
        R3 = zeros(N3+1);
        for j = 1:num_r
            temp = E3(:,j);
            temp1 = temp*temp';
            R3 = R3 + temp1;
        end
    
    [RV,RD] = eig(R3);
    for k = 1:length(RD)
        if abs(RD(k,k))<1e-14
            RD(k,k) = 1;
        end
    end
    RE = diag(RD);
    temp3 = - trace(RD*diag(log(RE)));
    ent3(n) = real(temp3);
end

% spin #1
num_r = (N2 + 1) * (N3 + 1) * (N1 + 1);
E4= zeros((N4+1),num_r);
ent4 = zeros(NN,1);
for n = 1:NN
        for i = 1:N4+1
            for j = 1:num_r
                pos = (j-1)*(N4+1)+i;
                E4(i,j) = V(pos,n);
            end
        end
        R4 = zeros(N4+1);
        for j = 1:num_r
            temp = E4(:,j);
            temp1 = temp*temp';
            R4 = R4 + temp1;
        end
    
    [RV,RD] = eig(R4);
    for k = 1:length(RD)
        if abs(RD(k,k))<1e-14
            RD(k,k) = 1;
        end
    end
    RE = diag(RD);
    temp3 = - trace(RD*diag(log(RE)));
    ent4(n) = real(temp3);
end

% save(strcat('main_N', num2str(N), '_p', num2str(p), '_h', num2str(h), '_No', num2str(No), '.mat'), 'N_s', 'sus', 'mm', 'phi4', 'h', 'e', 'r', 's', 'q', '-v7.3');
% 
figure;
plot(1:NN,ent4);
xlabel('\epsilon')
ylabel('entanglement')
% saveas(gcf, strcat('ratio_N', num2str(N), '_p', num2str(p), '_No', num2str(No), '.fig'));
% print(strcat('ratio_N', num2str(N), '_p', num2str(p), '_No', num2str(No)), '-dpng', '-r0');

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

function y = isdeg(a,b,c,d)
    y = 0;
    A = [a,b,c,d];
    A = sort(A);
    for i = 1:3
        for j = i+1:4
            t = A(j) - A(i);
            if t == 0
                y = 1;
                return
            end
        end
    end
end

function y = rearrange(A)
    y = A;
    [~,I] = max(A);
    if I == 4
        temp = y(1);
        y(1) = y(4);
        y(4) = temp;
    end
    if I == 2
        temp = y(1);
        y(1) = y(2);
        y(2) = y(4);
        y(4) = y(3);
        y(3) = temp;
    end
    if I == 3
        temp = y(1);
        y(1) = y(3);
        y(3) = y(4);
        y(4) = y(2);
        y(2) = temp;
    end       

end