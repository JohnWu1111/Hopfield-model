clear;
% close all;
clc;
format long
tic;

% Definition of parameters
N = 50:5:100; %size
lN = length(N);
J = 1;
J2 = 0.5;
h = 0.1;
res = 100;

corre_in_N = zeros(1,lN);
corre_ext_N = zeros(1,lN);
ent_N = zeros(1,lN);

corre_in_res = zeros(lN,res);
corre_ext_res = zeros(lN,res);
ent_res = zeros(lN,res);
count = zeros(lN,res);

for t = 1:lN
    
    N1 = N(t);
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
    
    % construction of Hamiltonian
    H1 = -(J*kron(S1_z.^2, ones(N1+1,1))...
        +J*kron(ones(N1+1,1), S1_z.^2)...
        +J2*kron(S1_z, S1_z))/(2*N1);
    H1 = diag(H1);
    H2 = h*(kron(S1_x, eye(N1+1))+kron(eye(N1+1), S1_x));
    H = H1 + H2;
    
    % time revolution
    [V,D] = eig(H);
    e = diag(D);
    NN = length(e);
    
    VV = V.^2;
    temp1 = kron(4*S1_z.^2, ones(N1+1,1))/(N1^2);
    corre_in = abs(temp1'*VV);
    temp2 = 4*kron(S1_z, S1_z)/(N1^2);
    corre_ext = abs(temp2'*VV);
    
    E1= zeros(N1+1,N1+1);
    ent = zeros(N1+1,1);
    for n = 1:NN
        for i = 1:N1+1
            for j = 1:N1+1
                pos = (i-1)*(N1+1)+j;
                E1(i,j) = V(pos,n);
            end
        end
        R1 = zeros(N1+1);
        for j = 1:N1+1
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
        ent(n) = real(temp3);
    end
    
    for i = 1:NN
        range = e(end) - e(1);
        temp = ceil((e(i) - e(1)) * res / range);
        if temp == 0
            temp = 1;
        end
        if temp > res
            temp = res;
        end
        
        corre_in_res(t,temp) = corre_in_res(t,temp) + corre_in(i);
        corre_ext_res(t,temp) = corre_ext_res(t,temp) + corre_ext(i);
        ent_res(t,temp) = ent_res(t,temp) + ent(i);
        count(t,temp) = count(t,temp) + 1;
    end
end

corre_in_res = corre_in_res./count;
corre_ext_res = corre_ext_res./count;
ent_res = ent_res./count;

corre_in_N = mean(corre_in_res(:,50:60),2);
corre_ext_N = mean(corre_ext_res(:,50:60),2);
ent_N = mean(ent_res(:,50:60),2);


figure;
plot(N,corre_in_N)
xlabel('size')
ylabel('internal correlation')

figure;
plot(N,corre_ext_N)
xlabel('size')
ylabel('external correlation')

figure;
plot(N,ent_N)
xlabel('size')
ylabel('entanglement')

% figure;
% plot(1/res:1/res:1,corre_in_res)
% xlabel('\epsilon')
% ylabel('internal correlation')
% 
% figure;
% plot(1/res:1/res:1,corre_ext_res)
% xlabel('\epsilon')
% ylabel('external correlation')
% 
% figure;
% plot(1/res:1/res:1,ent_res)
% xlabel('\epsilon')
% ylabel('entanglement')

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