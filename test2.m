load('H\H_N14_p2_No1.mat');
[row,col,v] = find(H);

N = 14;

temp1 = mem_con{1};
temp1a = round((-temp1 + 1) / 2) + 1;
temp1 = round((temp1 + 1) / 2) + 1;
spin1 = sigma_z(:, temp1(1));
spin1a = sigma_z(:, temp1a(1));

temp2 = mem_con{2};
temp2a = round((-temp2 + 1) / 2) + 1;
temp2 = round((temp2 + 1) / 2) + 1;
spin2 = sigma_z(:, temp2(1));
spin2a = sigma_z(:, temp2a(1));

for i = 2:N
    spin1 = kron(spin1, sigma_z(:, temp1(i)));
    spin1a = kron(spin1a, sigma_z(:, temp1a(i)));
    spin2 = kron(spin2, sigma_z(:, temp2(i)));
    spin2a = kron(spin2a, sigma_z(:, temp2a(i)));
end

[r1,c1] = find(spin1);
[r1a,c1a] = find(spin1);
[r2,c2] = find(spin2);
[r2a,c2a] = find(spin2a);

save('H_N14_p2_No1.txt','r1','r1a','r2','r2a','row','col','v','-ascii');