% calculation of some thermaldynamical observables

clear;
% close all;
clc;
format long
tic;

% Definition of parameters
N = 20; %size
dt = 0.1;
T = 1000;
t = 0:dt:T;
nt = length(t);
J = 1;
% h = 0:0.05:1;
h = 0.1;
lh = length(h);
No = 5;
cut = 1000;



% N1 = 2;
% N2 = 2;
% N3 = 2;
% N4 = 2;
% N = N1 + N2 + N3 + N4;
% N2ol = [1 1 -1 -1;
%     1 -1 1 -1;
%     1 -1 -1 1;
%     1 1 1 1];
% temp = N2ol*[N1,N2,N3,N4]';
% ol_12 = temp(1)/2;
% ol_13 = temp(2)/2;
% ol_23 = temp(3)/2;

% ol_12 = 1;
% ol_13 = 1;
% ol_23 = 0;
% N2ol = [1 1 -1 -1;
%     1 -1 1 -1;
%     1 -1 -1 1;
%     1 1 1 1];
% temp = N2ol\(2*[ol_12,ol_13,ol_23,N]');
% % temp = inv(N2ol)*[ol_12,ol_13,ol_23,N]';
% N1 = temp(1);
% N2 = temp(2);
% N3 = temp(3);
% N4 = temp(4);

load(strcat('p3_pattern\p3_nondeg_N',num2str(N),'.mat'));
ol = sortrows(ol,[5 6 7],'ComparisonMethod','abs');
N1 = ol(No,1);
N2 = ol(No,2);
N3 = ol(No,3);
N4 = ol(No,4);
% N1 = ol(end-No+1,1);
% N2 = ol(end-No+1,2);
% N3 = ol(end-No+1,3);
% N4 = ol(end-No+1,4);
% ol_12 = ol(end-No+1,5);
% ol_13 = ol(end-No+1,6);
% ol_23 = ol(end-No+1,7);
% ol_12 = ol(No,5);
% ol_13 = ol(No,6);
% ol_23 = ol(No,7);

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
    
    % spin1 = zeros(NN,1);
    % spin1a = zeros(NN,1);
    % spin1(1) = 1;
    % spin1a(end) = 1;
    %
    % temp1 = zeros((N1+1)*(N2+1),1);
    % temp2 = zeros((N3+1)*(N4+1),1);
    % temp1(1) = 1;
    % temp2(end) = 1;
    % spin2 = kron(temp1,temp2);
    % temp1 = zeros((N1+1)*(N2+1),1);
    % temp2 = zeros((N3+1)*(N4+1),1);
    % temp1(end) = 1;
    % temp2(1) = 1;
    % spin2a = kron(temp1,temp2);
    %
    % temp1 = zeros(N1+1,1);
    % temp2 = zeros(N2+1,1);
    % temp3 = zeros(N3+1,1);
    % temp4 = zeros(N4+1,1);
    % temp1(1) = 1;
    % temp2(end) = 1;
    % temp3(1) = 1;
    % temp4(end) = 1;
    % spin3 = kron4(temp1,temp2,temp3,temp4);
    % temp1 = zeros(N1+1,1);
    % temp2 = zeros(N2+1,1);
    % temp3 = zeros(N3+1,1);
    % temp4 = zeros(N4+1,1);
    % temp1(end) = 1;
    % temp2(1) = 1;
    % temp3(end) = 1;
    % temp4(1) = 1;
    % spin3a = kron4(temp1,temp2,temp3,temp4);
    
    % time revolution
    [V,D] = eig(H);
    e = diag(D);
    
    % spin4 = V(:,2500);
    %
    % spin0 = V'*spin4;
    % % spin0 = gpuArray(spin0);
    % % e = gpuArray(e);
    % % t = gpuArray(t);
    % trans = exp(-1i*e*t);
    % spin = spin0.*trans;
    % spin = gather(spin);
    % spint = V*spin;
    
    % ob = kron(ones((N1+1)*(N2+1)*(N3+1),1),S4_z.^2);
    % ob = diag(ob);
    % % ob = kron(eye((N1+1)*(N2+1)*(N3+1)),S4_x);
    % obb = sum(V.*(ob*V));
    % obb = obb./(N4^2);
    
%     for i = 1:len_p
%         if N_s(i) == 0
%             continue
%         end
%         temp = 1;
%         for k = 1:i-1
%             temp = kron(temp,ones(N_s(k)+1,1));
%         end
%         temp = kron(temp,4*S_z{i}.^2);
%         for k = i+1:len_p
%             temp = kron(temp,ones(N_s(k)+1,1));
%         end
% %         temp = sparse(diag(temp));
%         temp = diag(temp);
%         temp1 = sum(V(:,1:10).*(temp*V(:,1:10)));
%         temp2 = temp1 - N_s(i);
%         temp1 = N_s(i) + mean(temp2.^2)/(2*nchoosek(N_s(i),2));
%         sus(n) = sus(n) + temp1;
% %         sus(n) = sus(n) + mean(temp.^2);
%         
%         for j = i+1:len_p
%             if N_s(j) == 0
%                 continue
%             end
%             temp = 1;
%             for k = 1:i-1
%                 temp = kron(temp,ones(N_s(k)+1,1));
%             end
%             temp = kron(temp,2*S_z{i});
%             for k = i+1:j-1
%                 temp = kron(temp,ones(N_s(k)+1,1));
%             end
%             temp = kron(temp,2*S_z{j});
%             for k = j+1:len_p
%                 temp = kron(temp,ones(N_s(k)+1,1));
%             end
% %             temp = sparse(diag(temp));
%             temp = diag(temp);
%             temp = sum(V(:,1:10).*(temp*V(:,1:10)));
% %             sus(n) = sus(n) + mean(temp.^2);
%             sus(n) = sus(n) + 2*mean(temp.^2)/(N_s(i)*N_s(j));
%         end
%     end
%     
%     sus(n) = sus(n)/N;

    VV = V.^2;
    [B,I] = sort(N_s);
    a = min([I(1),I(2)]);
    b = max([I(1),I(2)]);
end

% figure;
% plot(1:NN,obb);
% sus = sus';

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