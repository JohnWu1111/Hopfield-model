clear;
% close all;
clc;
format long
tic;

% Definition of parameters
N = 5000; %size
lN = length(N);
J = 1;
h = 0.1;

N1 = N;
S1 = N1/2;

S1_z = zeros(N1+1,1);
S1_p = zeros(N1+1);
S1_m = zeros(N1+1);

% construction of matrice
for m = 1:N1+1
    S1_z(m) = N1/2 - (m-1);
end

for m = 1:N1
    S1_p(m,m+1) = sqrt(S1*(S1+1)-S1_z(m+1)*(S1_z(m+1)+1));
    S1_m(m+1,m) = sqrt(S1*(S1+1)-S1_z(m)*(S1_z(m)-1));
end

S1_x = (S1_p + S1_m)/2;
S1_y = (S1_p - S1_m)/2i;

% construction of Hamiltonian
H1 = -J*S1_z.^2/(2*N1);
H1 = diag(H1);
H2 = h*S1_x;
H = H1 + H2;

% time revolution
[V,D] = eig(H);
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
%
% figure;
% histogram(q,100,'Normalization','pdf','DisplayStyle','stairs');
%
% mean(q)

VV = V.^2;
Mz2 = S1_z.^2/N;
Mx = S1_x/N;
mz = Mz2'*VV;
mx = zeros(1,N+1);
for i = 1:N+1
    temp = V(:,i);
    mx(i) = temp'*Mx*temp;
end
% figure
% x = (e-e(1))/(e(end)-e(1));
% plot(x,mm);

VV = V.^2;


figure;
plot(1:N-1,q)
xlabel('index of eigenstates')
ylabel('ratio')

figure;
plot(1:N+1,mz)
xlabel('index of eigenstates')
ylabel('S_z^2/N')

figure;
plot(1:N+1,mx)
xlabel('index of eigenstates')
ylabel('S_x/N')

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