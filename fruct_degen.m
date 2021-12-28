clear;
% close all;
clc;
format long
tic;

% Definition of parameters
N = 40;% 56

p = 3;
J = 1;
h = 0;
de_ori = 8;
count = 0;

for i = N:-1:ceil(N/4)
    temp1 = min(N,2*i);
    for j = temp1:-1:i
        temp2 = min(N,2*j-i);
        temp3 = max(ceil((N+j)/2),j);
        for k = temp2:-1:temp3
%             if (i == k-i && N-k>0)|| (j-i == N-k+i && k-j>0)|| (k-j == N-k+i && j-i>0)|| (N-k == k-i && j-i>0)
            if 2*(i+k)== 3*N && N-k>0
                count = count +1;
            end            
        end
    end
end
GS = zeros(count,9);
ol = zeros(count,7);
% GS = zeros(floor((N/3+1)^3),1);
N1_s = zeros(count,1);
N2_s = zeros(count,1);
N3_s = zeros(count,1);
N4_s = zeros(count,1);

count = 0;
for i = N:-1:ceil(N/4)
    temp1 = min(N,2*i);
    for j = temp1:-1:i
        temp2 = min(N,2*j-i);
        temp3 = max(ceil((N+j)/2),j);
        for k = temp2:-1:temp3
%             if (i == k-i && N-k>0)|| (j-i == N-k+i && k-j>0)|| (k-j == N-k+i && j-i>0)|| (N-k == k-i && j-i>0)
            if 2*(i+k)== 3*N && N-k>0
                count = count +1;
                N1_s(count,1) = i;
                N2_s(count,1) = j-i;
                N3_s(count,1) = k-j;
                N4_s(count,1) = N-k;
            end
        end
    end
end

degen = zeros(count,1);
diff = zeros(count,de_ori-1);

for i = 1:count
    N1 = N1_s(i,1);
    N2 = N2_s(i,1);
    N3 = N3_s(i,1);
    N4 = N4_s(i,1);
    
    NN = (N1+1)*(N2+1)*(N3+1)*(N4+1);
    
    S1 = N1/2;
    S2 = N2/2;
    S3 = N3/2;
    S4 = N4/2;
    
    S1_z = zeros(N1+1,1);
    S2_z = zeros(N2+1,1);
    S3_z = zeros(N3+1,1);
    S4_z = zeros(N4+1,1);
    
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

    S1_p = zeros(N1+1);
    S1_m = zeros(N1+1);
    S2_p = zeros(N2+1);
    S2_m = zeros(N2+1);    
    S3_p = zeros(N3+1);
    S3_m = zeros(N3+1);
    S4_p = zeros(N4+1);
    S4_m = zeros(N4+1);
    
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
    
    H1 = -J*(3*kron(S1_z.^2,ones((N2+1)*(N3+1)*(N4+1),1))...
        +3*kron3(ones(N1+1,1),S2_z.^2,ones((N3+1)*(N4+1),1))...
        +3*kron3(ones((N1+1)*(N2+1),1),S3_z.^2,ones(N4+1,1))...
        +3*kron(ones((N1+1)*(N2+1)*(N3+1),1),S4_z.^2)...%inner term
        +2*kron3(S1_z,kron(S2_z,ones(N3+1,1))+kron(ones(N2+1,1),S3_z),ones(N4+1,1))...
        +2*kron3(ones(N1+1,1),kron(S2_z,ones(N3+1,1))+kron(ones(N2+1,1),S3_z),S4_z)...
        -2*kron4(ones(N1+1,1),S2_z,S3_z,ones(N4+1,1))...
        -2*kron3(S1_z,ones((N2+1)*(N3+1),1),S4_z))/(2*N);    
%     H2 = h*(kron(S1_x,eye((N2+1)*(N3+1)*(N4+1)))...
%         +kron3(eye(N1+1),S2_x,eye((N3+1)*(N4+1)))...
%         +kron3(eye((N1+1)*(N2+1)),S3_x,eye(N4+1))...
%         +kron(eye((N1+1)*(N2+1)*(N3+1)),S4_x));
%     
%     H1 = diag(H1);
%     H = H1 + H2;
    H = H1;

%     H = sparse(H);
%     [V,D] = eigs(H,de_ori,'smallestreal');
%     e = diag(D);
%     M = e(1);
%     ind = 0;
%     for j = 1:length(e)
%         if e(j) == M
%             ind = ind +1;
%         end
%     end
%     ind = ind/2;
%     for j = 1:length(e)-1
%         diff(i,j) = e(j+1) - e(j);
%     end
    
    ee = H(1:round(length(H)/2));
    [M,I] = min(ee);    
    ind = 0;
    for j = 1:length(ee)
        if ee(j) == M
            ind = ind +1;
        end
    end
    degen(i) = ind;
    
end
GS(:,1) = N1_s;
GS(:,2) = N2_s;
GS(:,3) = N3_s;
GS(:,4) = N4_s;

toc;

function y = kron4(a,b,c,d)
y = kron(kron(kron(a,b),c),d);
end

function y = kron3(a,b,c)
y = (kron(kron(a,b),c));
end
