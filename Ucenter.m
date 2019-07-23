function A = Ucenter(a,n)
m = sum(a,2);
M = sum(m)/((n - 1) * (n - 2));
m = m./(n-2);
A = a - repmat(m,[1 n]);
A = A - repmat(m,[1 n])';
A = A+M;
A(eye(size(A))==1)=0;
end