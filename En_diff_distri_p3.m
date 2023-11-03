clear;
% close all;
clc;
format long
tic;

N = 30; %size
p = 3;
num = 10;

q_save = cell(num,1);
e_save = cell(num,1);
s_save = cell(num,1);
count1 = 0;

N1 = 0;
N2 = 0;
N3 = 0;
N4 = 0;

while isdeg(N1, N2, N3, N4)
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
    
    %     if N1 == N2 && N2 == N3
    %         continue
    %     end
    %     if N1 == N2 && N2 == N4
    %         continue
    %     end
    %     if N1 == N3 && N3 == N4
    %         continue
    %     end
    
end

count1 = count1 + 1;
H = construt(N1,N2,N3,N4,N);
% H = gpuArray(H);
e = eig(H);
% e = gather(e);

len = length(e);
%     e = e(1:20);
%     e = e(round(len*0.1):round(len*0.3));
%     e = e(round(len/3):round(2*len/3));

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

%     qq(count+1:count+len-2,1) = q;
%     count = count + len-2;
q_save{count1,1} = q;
e_save{count1,1} = e;
s_save{count1,1} = s;

% s = sort(s);

count2 = 0;
for i = 1:count1
    len = length(q_save{i});
    qq(count2+1:count2+len,1) = q_save{i};
    count2 = count2 + len;
end

figure;
histogram(qq,200,'Normalization','pdf','DisplayStyle','stairs');
% axis([0 0.02 0 inf])

mean(qq)

%save(strcat('En_diff_N', num2str(N), '_p', num2str(p), '_h0.1', '.mat'), '-v7.3');

% saveas(gcf,strcat('Energy_diff_N',num2str(N),'_p',num2str(p),'_No',num2str(No),'.fig'));
% print(strcat('Energy_diff_N',num2str(N),'_p',num2str(p),'_No',num2str(No)),'-dpng','-r0');
% save(strcat('Energy_diff_N',num2str(N),'_p',num2str(p),'_No',num2str(No),'.mat'),'q','s','e','-v7.3');

toc;

function y = construt(N1,N2,N3,N4,N)

hx = 2;
hy = 0;
J = 1;

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

S1_y = (S1_p - S1_m)/(2*1i);
S2_y = (S2_p - S2_m)/(2*1i);
S3_y = (S3_p - S3_m)/(2*1i);
S4_y = (S4_p - S4_m)/(2*1i);

% construction of Hamiltonian
H1 = -J*(3*kron(S1_z.^2,ones((N2+1)*(N3+1)*(N4+1),1))...
    +3*kron3(ones(N1+1,1),S2_z.^2,ones((N3+1)*(N4+1),1))...
    +3*kron3(ones((N1+1)*(N2+1),1),S3_z.^2,ones(N4+1,1))...
    +3*kron(ones((N1+1)*(N2+1)*(N3+1),1),S4_z.^2)...%inner term
    +2*kron3(S1_z,kron(S2_z,ones(N3+1,1))+kron(ones(N3+1,1),S2_z),ones(N4+1,1))...
    +2*kron3(ones(N1+1,1),kron(S2_z,ones(N3+1,1))+kron(ones(N3+1,1),S2_z),S4_z)...
    -2*kron4(ones(N1+1,1),S2_z,S3_z,ones(N4+1,1))...
    -2*kron3(S1_z,ones((N2+1)*(N3+1),1),S4_z))/(2*N);
H1 = diag(H1);
H2 = hx*(kron(S1_x,eye((N2+1)*(N3+1)*(N4+1)))...
    +kron3(eye(N1+1),S2_x,eye((N3+1)*(N4+1)))...
    +kron3(eye((N1+1)*(N2+1)),S3_x,eye(N4+1))...
    +kron(eye((N1+1)*(N2+1)*(N3+1)),S4_x));
H3 = hy*(kron(S1_y,eye((N2+1)*(N3+1)*(N4+1)))...
    +kron3(eye(N1+1),S2_y,eye((N3+1)*(N4+1)))...
    +kron3(eye((N1+1)*(N2+1)),S3_y,eye(N4+1))...
    +kron(eye((N1+1)*(N2+1)*(N3+1)),S4_y));
H = H1 + H2 + H3;
y = H;

end

function y = kron4(a,b,c,d)
y = kron(kron(kron(a,b),c),d);
end

function y = kron3(a,b,c)
y = (kron(kron(a,b),c));
end

function y = isdeg(a,b,c,d)
y = 0;
if a == d || b == c
    y = 1;
    return
end
if a == b && c == d
    y = 1;
    return
end
if a == c && b == d
    y = 1;
    return
end
end