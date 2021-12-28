I = 4;

I4 = mod(I-1,N4+1)+1;
II = (I-I4)/(N4+1) + 1;
I3 = mod(II-1,N3+1)+1;
II = (II-I3)/(N3+1) + 1;
I2 = mod(II-1,N2+1)+1;
I1 = (II-I2)/(N2+1) + 1;

s1 = zeros(N1+1,1);
s2 = zeros(N2+1,1);
s3 = zeros(N3+1,1);
s4 = zeros(N4+1,1);
s1(1) = 1;
s2(1) = 1;
s3(1) = 1;
s4(end) = 1;

ss = kron4(s1,s2,s3,s4);

function y = kron4(a,b,c,d)
y = kron(kron(kron(a,b),c),d);
end