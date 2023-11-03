clear;
% close all;
clc;
format long
tic;

% Definition of parameters
N = 22;% 56

p = 3;
J = 1;
h = 0;
count = 0;

for i = N:-1:ceil(N/4)
    temp1 = min(N,2*i);
    for j = temp1:-1:i
        temp2 = min(N,2*j-i);
        temp3 = max(ceil((N+j)/2),j);
        for k = temp2:-1:temp3
            if i == j-i || i == k-j || i == N-k || j-i == k-j || j-i == N-k || k-j == N-k
                continue
            end
            count = count +1;
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
            if i == j-i || i == k-j || i == N-k || j-i == k-j || j-i == N-k || k-j == N-k
                continue
            end
            count = count +1;
            N1_s(count,1) = i;
            N2_s(count,1) = j-i;
            N3_s(count,1) = k-j;
            N4_s(count,1) = N-k;
        end
    end
end

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
%     S1_p = zeros(N1+1,1);
%     S1_m = zeros(N1+1,1);
    S2_z = zeros(N2+1,1);
%     S2_p = zeros(N2+1,1);
%     S2_m = zeros(N2+1,1);
    S3_z = zeros(N3+1,1);
%     S3_p = zeros(N3+1,1);
%     S3_m = zeros(N3+1,1);
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
    
%     M1 = diag(S1_z);
%     M2 = diag(S2_z);
%     M3 = diag(S3_z);
    
%     for m = 1:N1
%         S1_p(m,m+1) = sqrt(S1*(S1+1)-M1(m+1)*(M1(m+1)+1));
%         S1_m(m+1,m) = sqrt(S1*(S1+1)-M1(m)*(M1(m)-1));
%     end
%     for m = 1:N2
%         S2_p(m,m+1) = sqrt(S2*(S2+1)-M2(m+1)*(M2(m+1)+1));
%         S2_m(m+1,m) = sqrt(S2*(S2+1)-M2(m)*(M2(m)-1));
%     end
%     for m = 1:N3
%         S3_p(m,m+1) = sqrt(S3*(S3+1)-M3(m+1)*(M3(m+1)+1));
%         S3_m(m+1,m) = sqrt(S3*(S3+1)-M3(m)*(M3(m)-1));
%     end
    
%     S1_x = (S1_p + S1_m)/2;
%     S2_x = (S2_p + S2_m)/2;
%     S3_x = (S3_p + S3_m)/2;
    H = -J*(3*kron(S1_z.^2,ones((N2+1)*(N3+1)*(N4+1),1))...
        +3*kron3(ones(N1+1,1),S2_z.^2,ones((N3+1)*(N4+1),1))...
        +3*kron3(ones((N1+1)*(N2+1),1),S3_z.^2,ones(N4+1,1))...
        +3*kron(ones((N1+1)*(N2+1)*(N3+1),1),S4_z.^2)...%inner term
        +2*kron3(S1_z,kron(S2_z,ones(N3+1,1))+kron(ones(N2+1,1),S3_z),ones(N4+1,1))...
        +2*kron3(ones(N1+1,1),kron(S2_z,ones(N3+1,1))+kron(ones(N2+1,1),S3_z),S4_z)...
        -2*kron4(ones(N1+1,1),S2_z,S3_z,ones(N4+1,1))...
        -2*kron3(S1_z,ones((N2+1)*(N3+1),1),S4_z))/(2*N);%exchange term
    %                 +h*(kron(kron(S1_x,eye(N2+1)),eye(N3+1))...
    %                 +kron(kron(eye(N1+1),S2_x),eye(N3+1))...
    %                 +kron(kron(eye(N1+1),eye(N2+1)),S3_x));
    
%     [V,D] = eig(H);
%     H = sparse(H);
%     [V,D] = eigs(H,1,'smallestabs');
    %             e = diag(D);
    e = H;
    ee = e(1:round(length(e)/2));
    [M,I] = min(ee);
    I4 = mod(I-1,N4+1)+1;
    II = (I-I4)/(N4+1) + 1;
    I3 = mod(II-1,N3+1)+1;
    II = (II-I3)/(N3+1) + 1;
    I2 = mod(II-1,N2+1)+1;
    I1 = (II-I2)/(N2+1) + 1;

    GS(i,5) = I1;
    GS(i,6) = I2;
    GS(i,7) = I3;
    GS(i,8) = I4;
    GS(i,9) = NN;
%     ol(i,2) = N1;
%     ol(i,3) = N2;
%     ol(i,4) = N3;
%     GS(1:NN,i) = V(:,1);
    
end
GS(:,1) = N1_s;
GS(:,2) = N2_s;
GS(:,3) = N3_s;
GS(:,4) = N4_s;

temp1 = GS(:,1:4);
N2ol = [1 1 -1 -1;
    1 -1 1 -1;
    1 -1 -1 1;
    1 1 1 1];
temp2 = N2ol*temp1'/2;
ol(:,1:4) = GS(:,1:4);
ol(:,5:7) = temp2(1:3,:)';

ol = sortrows(ol,[1 2 3 4],'ComparisonMethod','abs');

save(strcat('p3_pattern\p3_nondeg_N',num2str(N),'.mat'),'GS','ol','count','-v7.3');

toc;

function y = kron4(a,b,c,d)
y = kron(kron(kron(a,b),c),d);
end

function y = kron3(a,b,c)
y = (kron(kron(a,b),c));
end
