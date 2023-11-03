clear;
% close all;
clc;
format long
tic;

% Definition of parameters
% N = 56; %size
J = 1;
h = 0.1;
num = 1;

pro = zeros(num,3);

for n = 1:1
    
    N1 = 20;
    N2 = 10;
    N3 = 8;
    N4 = 6;
    
    N = N1+N2+N3+N4;
    
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
    H2 = h*(kron(S1_x,eye((N2+1)*(N3+1)*(N4+1)))...
        +kron3(eye(N1+1),S2_x,eye((N3+1)*(N4+1)))...
        +kron3(eye((N1+1)*(N2+1)),S3_x,eye(N4+1))...
        +kron(eye((N1+1)*(N2+1)*(N3+1)),S4_x));
    H = H1 + H2;
%     H = H1;
    
    spin1 = zeros(NN,1);
    spin1a = zeros(NN,1);
    spin1(1) = 1;
    spin1a(end) = 1;
    
    temp1 = zeros((N1+1)*(N2+1),1);
    temp2 = zeros((N3+1)*(N4+1),1);
    temp1(1) = 1;
    temp2(end) = 1;
    spin2 = kron(temp1,temp2);
    temp1 = zeros((N1+1)*(N2+1),1);
    temp2 = zeros((N3+1)*(N4+1),1);
    temp1(end) = 1;
    temp2(1) = 1;
    spin2a = kron(temp1,temp2);
    
    temp1 = zeros(N1+1,1);
    temp2 = zeros(N2+1,1);
    temp3 = zeros(N3+1,1);
    temp4 = zeros(N4+1,1);
    temp1(1) = 1;
    temp2(end) = 1;
    temp3(1) = 1;
    temp4(end) = 1;
    spin3 = kron4(temp1,temp2,temp3,temp4);
    temp1 = zeros(N1+1,1);
    temp2 = zeros(N2+1,1);
    temp3 = zeros(N3+1,1);
    temp4 = zeros(N4+1,1);
    temp1(end) = 1;
    temp2(1) = 1;
    temp3(end) = 1;
    temp4(1) = 1;
    spin3a = kron4(temp1,temp2,temp3,temp4);
    
%     H = sparse(H);
%     [V,D] = eigs(H,1,'smallestreal');
%     pro1 = spin1'*V(:,1);
%     pro2 = spin2'*V(:,1);
%     pro3 = spin3'*V(:,1);
%     pro(n,:) = abs([pro1 pro2 pro3]);
    
%     e = H;
%     E1 = H(1);
%     E2 = spin2'*e;
%     E3 = spin3'*e;
%     e = sort(e); 
    [V,D] = eigs(H,16,'smallestreal');
    
end

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