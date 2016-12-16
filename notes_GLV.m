% R. May, Will al Large Complex System be Stable? Nature 1972 setup
C=0.5;
alpha=0.1;
n=100;
d=1;
M=alpha.*randn(n);
M=M.*(rand(n)<C);
for i=1:n
    M(i,i)=-d;
end;
% Eigs can only calculate up to n-2
% eigenvalues/eigenvectors of a non-symmetric matrix
[V,D,FLAG] = eigs(M,n-2);
Dr=real(diag(D));
Di=imag(diag(D));
figure; plot(Dr,Di,'ko');
hold on;
radius=alpha.*sqrt(n*C);
x=-radius:0.001:radius;
y=sqrt(radius.^2-x.^2);
plot(x-d,y,'k-');
plot(x-d,-y,'k-');
%%
% GLV setup
% dx_i/dt=x_i(g_i+sum_j M_ij*x_j);
% count how many out of 2^n states is feasible
C=0.5;
alpha=0.1;
n=100;
d=1;
M=alpha.*randn(n);
M=M.*(rand(n)<C);
for i=1:n
    M(i,i)=-d;
end;
g=ones(n,1);
pa_array=0.1:0.1:0.9;
k_array=length(pa_array);
Stats=100;
n_act=zeros(k_array, 1);
sum_neg=zeros(k_array, 1);
f_feas=zeros(k_array, 1);
for k=1:k_array;
 pa=pa_array(k);
for t=1:Stats;
active=ones(n,1).*(rand(n,1)<pa);
n_act(k)=n_act(k)+sum(active);
ia=find(active);
Ma=M(ia,ia);
ga=g(ia);
xa=inv(Ma)*(-ga);
neg_a=sum(xa<=0);
sum_neg(k)=sum_neg(k)+neg_a;
if neg_a==0;
    f_feas(k)=f_feas(k)+1;
end;
end;
end;
n_act=n_act./Stats;
f_feas=f_feas./Stats;
sum_neg=sum_neg./Stats;
hold on;
%plot(n_act, f_feas,'ko-')
figure; plot(n_act, f_feas,'ko-')
figure; plot(n_act, sum_neg,'ko-')
%%