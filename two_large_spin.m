% calculation of some thermaldynamical observables

clear;
% close all;
format long
tic;

% Definition of parameters
J = 1;
N_a = 20:4:100;
h = 2;
lN = length(N_a);
ent = zeros(1,lN);

for t = 1:lN
    N = N_a(t);
    N1 = N/2;
    N2 = N/2;
    
    NN = (N1 + 1) * (N2 + 1);
    
    S1 = N1 / 2;
    S2 = N2 / 2;
    
    S1_z = zeros(N1 + 1, 1);
    S1_p = zeros(N1 + 1);
    S1_m = zeros(N1 + 1);
    S2_z = zeros(N2 + 1, 1);
    S2_p = zeros(N2 + 1);
    S2_m = zeros(N2 + 1);
    
    % construction of matrice
    for m = 1:N1 + 1
        S1_z(m) = N1 / 2 - (m - 1);
    end
    
    for m = 1:N2 + 1
        S2_z(m) = N2 / 2 - (m - 1);
    end
    
    for m = 1:N1
        S1_p(m, m + 1) = sqrt(S1 * (S1 + 1) - S1_z(m + 1) * (S1_z(m + 1) + 1));
        S1_m(m + 1, m) = sqrt(S1 * (S1 + 1) - S1_z(m) * (S1_z(m) - 1));
    end
    
    for m = 1:N2
        S2_p(m, m + 1) = sqrt(S2 * (S2 + 1) - S2_z(m + 1) * (S2_z(m + 1) + 1));
        S2_m(m + 1, m) = sqrt(S2 * (S2 + 1) - S2_z(m) * (S2_z(m) - 1));
    end
    
    S1_x = (S1_p + S1_m) / 2;
    S2_x = (S2_p + S2_m) / 2;
    
    N_s = [N1 N2];
    S_z = {S1_z, S2_z};
    
    % construction of Hamiltonian
    H1 = -J * (2*kron(S1_z.^2, ones((N2 + 1), 1)) ...
        + 2*kron(ones(N1 + 1, 1), S2_z.^2) ...%inner term
        + kron(S1_z, S2_z)) / (2 * N);
    H1 = diag(H1);
    H2 = h * (kron(S1_x, eye((N2 + 1))) ...
        +kron(eye(N1 + 1), S2_x));
    H = H1 + H2;
    
    clear H1 H2;
    
    [V, D] = eig(H);
    e = diag(D);
    
    clear H D;
    
    VV = V.^2;
    
    num_r = (N2 + 1);
    E12 = zeros((N1+1),num_r);
    S12 = zeros(NN,1);
    for n = 1:NN
        for i = 1:N1+1
            for j = 1:num_r
                pos = (i-1)*num_r+j;
                E12(i,j) = V(pos,n);
            end
        end
        R12 = zeros(N1+1);
        for j = 1:num_r
            temp = E12(:,j);
            temp1 = temp*temp';
            R12 = R12 + temp1;
        end
        
        [RV,RD] = eig(R12);
        for k = 1:length(RD)
            if abs(RD(k,k))<1e-14
                RD(k,k) = 1;
            end
        end
        RE = diag(RD);
        temp3 = - trace(RD*diag(log(RE)));
        S12(n) = real(temp3);
    end
    
    count = 0;
    range = e(end)-e(1);
    for i = 1:NN
        pos = (e(i)-e(1))/range;
        if pos > 0.45 && pos < 0.55
            ent(t) = ent(t) + S12(i);
            count = count +1;
        end
    end
    ent(t) = ent(t)/count;
end

% figure;
% plot(1:NN,S12);
% xlabel('\epsilon')
% ylabel('entanglement')

% figure;
% plot(N_a,ent);
% xlabel('N')
% ylabel('entanglement')

figure;
loglog(N_a,ent);
xlabel('N')
ylabel('entanglement')

x = log(N_a);
y = log(ent);

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