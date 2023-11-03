% calculation of some thermaldynamical observables

clear;
% close all;
clc;
format long
tic;

% Definition of parameters
N = 20; %size
p = 3;
J = 1;
% h = 0:0.05:1;
h = 100;
lh = length(h);
No = 5;

mem_con = cell(p,1);
mem_con{1} = ones(1,N);
for k = 2:p
    mem_con{k} = round(rand(1,N))*2-1;
end

N1 = 0;
N2 = 0;
N3 = 0;
N4 = 0;
for i = 1:N
    if mem_con{2}(i)*mem_con{3}(i) > 0
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

% if N1 == N2 || N1 == N3 || N1 == N4 || N2 == N3 || N2 == N4 || N3 == N4
%     continue
% end

NN = (N1+1)*(N2+1)*(N3+1)*(N4+1);

S1 = N1/2;
S2 = N2/2;
S3 = N3/2;
S4 = N4/2;

S1_z = zeros(N1+1,1);
S1_p = zeros(N1+1);
S1_m = zeros(N1+1);
S2_z = zeros(N2+1,1);
S2_p = zeros(N2+1);
S2_m = zeros(N2+1);
S3_z = zeros(N3+1,1);
S3_p = zeros(N3+1);
S3_m = zeros(N3+1);
S4_z = zeros(N4+1,1);
S4_p = zeros(N4+1);
S4_m = zeros(N4+1);

% construction of matrice
for m = 1:N1+1
    S1_z(m) = N1/2 - (m-1);
end
for m = 1:N2+1
    S2_z(m) = N2/2 - (m-1);
end
for m = 1:N3+1
    S3_z(m) = N3/2 - (m-1);
end
for m = 1:N4+1
    S4_z(m) = N4/2 - (m-1);
end

for m = 1:N1
    S1_p(m,m+1) = sqrt(S1*(S1+1)-S1_z(m+1)*(S1_z(m+1)+1));
    S1_m(m+1,m) = sqrt(S1*(S1+1)-S1_z(m)*(S1_z(m)-1));
end
for m = 1:N2
    S2_p(m,m+1) = sqrt(S2*(S2+1)-S2_z(m+1)*(S2_z(m+1)+1));
    S2_m(m+1,m) = sqrt(S2*(S2+1)-S2_z(m)*(S2_z(m)-1));
end
for m = 1:N3
    S3_p(m,m+1) = sqrt(S3*(S3+1)-S3_z(m+1)*(S3_z(m+1)+1));
    S3_m(m+1,m) = sqrt(S3*(S3+1)-S3_z(m)*(S3_z(m)-1));
end
for m = 1:N4
    S4_p(m,m+1) = sqrt(S4*(S4+1)-S4_z(m+1)*(S4_z(m+1)+1));
    S4_m(m+1,m) = sqrt(S4*(S4+1)-S4_z(m)*(S4_z(m)-1));
end

S1_x = (S1_p + S1_m)/2;
S2_x = (S2_p + S2_m)/2;
S3_x = (S3_p + S3_m)/2;
S4_x = (S4_p + S4_m)/2;

sus = zeros(1,lh);
len_p = 4;
N_s = [N1 N2 N3 N4];
S_z = {S1_z,S2_z,S3_z,S4_z};
% Jij = [3 2 2 -2;
%        0 3 -2 2;
%        0 0 3 2;
%        0 0 0 3];
Jij = [3 2 2 2;
    0 3 2 2;
    0 0 3 2;
    0 0 0 3];
for n = 1:lh
    % construction of Hamiltonian
    H1 = -J*(3*kron(S1_z.^2,ones((N2+1)*(N3+1)*(N4+1),1))...
        +3*kron3(ones(N1+1,1),S2_z.^2,ones((N3+1)*(N4+1),1))...
        +3*kron3(ones((N1+1)*(N2+1),1),S3_z.^2,ones(N4+1,1))...
        +3*kron(ones((N1+1)*(N2+1)*(N3+1),1),S4_z.^2)...%inner term
        +2*kron3(S1_z,kron(S2_z,ones(N3+1,1))+kron(ones(N2+1,1),S3_z),ones(N4+1,1))...
        +2*kron3(ones(N1+1,1),kron(S2_z,ones(N3+1,1))+kron(ones(N2+1,1),S3_z),S4_z)...
        -2*kron4(ones(N1+1,1),S2_z,S3_z,ones(N4+1,1))...
        -2*kron3(S1_z,ones((N2+1)*(N3+1),1),S4_z))/(2*N);
    H1 = diag(H1);
    H2 = h(n)*(kron(S1_x,eye((N2+1)*(N3+1)*(N4+1)))...
        +kron3(eye(N1+1),S2_x,eye((N3+1)*(N4+1)))...
        +kron3(eye((N1+1)*(N2+1)),S3_x,eye(N4+1))...
        +kron(eye((N1+1)*(N2+1)*(N3+1)),S4_x));
    H = H1 + H2;
    
    
    [V,D] = eig(H);
    e = diag(D);
    
    sus_it = zeros(1,NN);
    
    for i = 1:len_p
        if N_s(i) == 0
            continue
        end
        temp = 1;
        for k = 1:i-1
            temp = kron(temp,ones(N_s(k)+1,1));
        end
        temp = kron(temp,4*S_z{i}.^2);
        for k = i+1:len_p
            temp = kron(temp,ones(N_s(k)+1,1));
        end
        temp = sparse(diag(temp));
%         temp = diag(temp);
        temp1 = sum(V.*(temp*V));
        if N_s(i) > 1
            temp2 = temp1 - N_s(i);
            temp1 = N_s(i) + temp2.^2/(2*nchoosek(N_s(i),2));
        else
            temp1 = 1;
        end
        sus_it = sus_it + temp1;
%         sus(n) = sus(n) + mean(temp.^2);
        
        for j = i+1:len_p
            if N_s(j) == 0
                continue
            end
            temp = 1;
            for k = 1:i-1
                temp = kron(temp,ones(N_s(k)+1,1));
            end
            temp = kron(temp,2*S_z{i});
            for k = i+1:j-1
                temp = kron(temp,ones(N_s(k)+1,1));
            end
            temp = kron(temp,2*S_z{j});
            for k = j+1:len_p
                temp = kron(temp,ones(N_s(k)+1,1));
            end
            temp = sparse(diag(temp));
%             temp = diag(temp);
            temp = sum(V.*(temp*V));
%             sus(n) = sus(n) + mean(temp.^2);
            sus_it = sus_it + 2*temp.^2/(N_s(i)*N_s(j));
        end
    end
    sus_it = sus_it/N;
    sus(n) = mean(sus_it);
end

figure;
plot(h,sus);
sus = sus';

toc;

function y = kron4(a,b,c,d)
y = kron(kron(kron(a,b),c),d);
end

function y = kron3(a,b,c)
y = (kron(kron(a,b),c));
end

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

function y = kron_p4(a,b,c,d)
y = kron_p(kron_p(kron_p(a,b),c),d);
end