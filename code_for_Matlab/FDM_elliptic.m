function solution=FDM_elliptic(Diff,deg,f)
% 二维的possion方程求解
% -D*Delta u = f-d*u   边界(a,b)*(c,d)
%   Neumann boundary condition

% Diff=D; deg=d

global mLap hx hy

a = 0;
b = 100;
c = 0;
d = 100;
M = 100;
N = 100;
hx = (b-a)/M;
hy = (d-c)/N;
x = (a+hx:hx:b)';
y = (c+hy:hy:d)';
[X Y] = meshgrid(x,y);
% x_in = (a+hx:hx:b-hx);
% y_in = (c+hy:hy:d-hy);
% % [X_in Y_in] = meshgrid(x_in,y_in);

%生成Dxx
d_0 = -2*ones(M-1,1);
d_1 = ones(M-1,1);
d_m1 = d_1;
d_xx = spdiags([d_m1 d_0 d_1],[-1,0,1],M-1,M-1);
%% incprperating boundary conditions
d_xx(1,2)=2; d_xx(M-2,M-1)=2; 
%%
d_xx = d_xx/hx^2;
I_N = speye(N-1,N-1);
D_xx = kron(I_N,d_xx);

%D_yy
d_0 = -2*ones(N-1,1);
d_1 = ones(N-1,1);
d_m1 = d_1;
d_yy = spdiags([d_m1 d_0 d_1],[-1,0,1],N-1,N-1);
%% incorporating boundary conditions
d_yy(1,2)=2; d_yy(N-2,N-1)=2; 
%%
d_yy = d_yy/hy^2;
I_M = speye(M-1,M-1);
D_yy = kron(d_yy,I_M);

%左端矩阵
mLap = -(D_xx+D_yy);

%右端矩阵
F = f(2:N,2:M);
F = F';

%% First class boundary condition
% g_d = g(x_in,0);
% g_u = g(x_in,1);
% g_l = g(0,y_in);
% g_r = g(1,y_in);
% F(:,1) = F(:,1)+g_d'/hy^2;
% F(:,N-1) = F(:,N-1)+g_u'/hy^2;
% F(1,:) = F(1,:)+g_l/hx^2;
% F(M-1,:) = F(M-1,:)+g_r/hx^2;
%%
%求解u
ff = F(:);

u = (Diff*mLap+deg*speye(size(Diff*mLap)))\ff;
U = reshape(u,M-1,N-1);
U = U';

solution=zeros(M,N);
solution(2:N,2:M) = U;
solution(1,:)=solution(2,:);
solution(:,1)=solution(:,2);

%绘制结果
% 
% figure,
% mesh(X,Y,solution);

end

