clear all
clc
%% set parameters
xl=0; xr=1; yl=0; yr=1; %[xl,xr]x[yl,yr] is the the computational domain
N=2; %4M=2*N*(N+1) is the size of quadrature set 
I=32; J=I;
hx=(xr-xl)/I; hy=(yr-yl)/J;
tol=2*log(10);% delta=10^(-n), then tol=n*log(10)
parameters={N I J xl xr yl yr tol};

%% cross section functions, source term
% Example 1
g=0;
f_sigma_T=@(x,y)(1).*(x<=xr).*(y<=yr);
f_sigma_a=@(x,y)(0.5).*(x<=xr).*(y<=yr);
f_q=@(x,y) (0).*(x<=xr).*(y<=yr);
f_varepsilon=@(x,y) 0.001+0.999.*(x>=0.25).*(x<0.5).*(y>=0).*(y<0.25)+0.999.*(x>=0.75).*(x<1).*(y>=0).*(y<0.25)+...
    0.999.*(x>=0).*(x<0.25).*(y>=0.25).*(y<0.5)+0.999.*(x>=0.5).*(x<0.75).*(y>=0.25).*(y<0.5)+...
    0.999.*(x>=0.25).*(x<0.5).*(y>=0.5).*(y<0.75)+0.999.*(x>=0.75).*(x<1).*(y>=0.5).*(y<0.75)+...
    0.999.*(x>=0).*(x<0.25).*(y>=0.75).*(y<1)+0.999.*(x>=0.5).*(x<0.75).*(y>=0.75).*(y<1);

% Example 2
% g=0.2;
% % f_varepsilon=@(x,y)(0.02.*(x./xr)+0.001).*(x<=xr).*(y<=yr);
% f_sigma_T=@(x,y)(1+x.^2+y.^2)./(0.02.*x+0.001).*(y<=yr);
% f_sigma_a=@(x,y)(0.5+x.^2+y.^2).*(0.02.*x+0.001).*(x<=xr).*(y<=yr);
% f_q=@(x,y) sin(x.*y).*(0.02.*x+0.001);
% f_varepsilon=@(x,y) (1).*(x<=xr).*(y<=yr);

% Example 3
% g=0;
% f_sigma_T=@(x,y)(-10+10.*(x.^2+y.^2)).*(x<=0.5).*(y<=0.5)+11;
% f_sigma_a=@(x,y)(-4+10.*(x.^2+y.^2)).*(x<=0.5).*(y<=0.5)+5;
% f_q=@(x,y)(5).*(x<=0.5).*(y<=0.5);
% f_varepsilon=@(x,y) (1).*(x<=xr).*(y<=yr);


% Example 4
% g=0;
% f_sigma_T=@(x,y)(0.999*(x.^2+y.^2).^2.*(x.^2+y.^2-2).^2+0.001).*((x.^2+y.^2)<1)+1.*((x.^2+y.^2)>=1)+20;
% f_sigma_a=@(x,y) (0).*(x<=xr).*(y<=yr);
% f_q=@(x,y) 1/(4*pi*0.01).*exp(-(x.^2+y.^2)/(4*0.01)).*20.*(x<=xr).*(y<=yr);
% f_varepsilon=@(x,y) (0.01).*(x<=xr).*(y<=yr);

[omega,ct,st,M,the]=qnwlege2(N); % generate quadrature set
[ Kappa, val ] = HG( N, g );% only for HG kernel

Sigma_T=zeros(I,J); Sigma_a=zeros(I,J); Q=zeros(I,J); Varepsilon=zeros(I,J);
for i=1:I
    for j=1:J
        Sigma_T(i,j)=integral2(f_sigma_T,(i-1)*hx+xl,i*hx+xl,(j-1)*hy+yl,j*hy+yl)/hx/hy;
        Sigma_a(i,j)=integral2(f_sigma_a,(i-1)*hx+xl,i*hx+xl,(j-1)*hy+yl,j*hy+yl)/hx/hy;
        Q(i,j)=integral2(f_q,(i-1)*hx+xl,i*hx+xl,(j-1)*hy+yl,j*hy+yl)/hx/hy;
        Varepsilon(i,j)=integral2(f_varepsilon,(i-1)*hx+xl,i*hx+xl,(j-1)*hy+yl,j*hy+yl)/hx/hy;
    end
end

sigma_T=Sigma_T./Varepsilon; sigma_a=Sigma_a.*Varepsilon; q=Q.*Varepsilon; 

coefficients={Kappa,sigma_T,sigma_a,q};

%% boundary condition
% 直接给边界条件 <--每个边界点上入射值相同

% psiL0=[zeros(M,1);zeros(M,1)];
% psiR0=[zeros(M,1);zeros(M,1)];
% psiB0=[zeros(M,1);zeros(M,1)];
% psiT0=[zeros(M,1);zeros(M,1)];

% psiL0=[ones(M,1);zeros(M,1)];
% psiR0=[ones(M,1);zeros(M,1)];
% psiB0=[zeros(M,1);ones(M,1)];
% psiT0=[zeros(M,1);ones(M,1)];


psiL0=[ones(M,1);ones(M,1)];
psiR0=[ones(M,1);ones(M,1)];
psiB0=[ones(M,1);ones(M,1)];
psiT0=[ones(M,1);ones(M,1)];


%get boundary data for level k mesh
psiL=kron(ones(J,1),psiL0);
psiR=kron(ones(J,1),psiR0);
psiB=kron(ones(I,1),psiB0);
psiT=kron(ones(I,1),psiT0);

bc={psiL,psiR,psiB,psiT};

%% offline assemble
% assemble corresponding matrices and show the existence of low rankness
[A,Dc,VecSize,Fsm,Ml,Mr,Mb,Mt] = assemble(parameters,coefficients);
%% online solution
[b,b_full,Q,Q100]=solution(parameters,coefficients,bc,Ml,Mr,Mb,Mt);
%% phi
alpha=A\b; % coefficients of Adaptive TFPS basis 
psi=Dc*alpha+Q; % angular flux at cell centers
phi=zeros(I,J); % scalar flux at cell centers
for j=1:J
    for i=1:I
        phi(i,j)=omega'*psi(((j-1)*I+(i-1))*4*M+(1:4*M));
    end
end

% psi=Dc100*alpha;
% phi=zeros(10*6,1);
% for i=1:60
%     phi(i)=omega'*psi((i-1)*4*M+(1:4*M));
% end
%% plot solution at cell centers
[xmesh,ymesh]=meshgrid(xl+1/2*hx:hx:xr-1/2*hx,yl+1/2*hy:hy:yr-1/2*hy);
figure(1)
surf(xmesh,ymesh,phi');
xlabel('x'); ylabel('y');
view(2)
shading interp
colorbar
set(gca,'fontsize',14)

%% plot the number of basis functions in each cell
figure(4)
% 颜色定义
colormap default
% 画图
% figureHandle = figure;
b = bar3(VecSize,1);
% zdata作为cdata的值，也就是使用高度赋色
for n=1:numel(b)
   cdata=get(b(n),'zdata');
   cdata=repmat(max(cdata,[],2),1,4);
   set(b(n),'cdata',cdata,'facecolor','flat')
end
% 添加标题、坐标轴标签
% hTitle = title('Adaptive selection of basis functions for transport regime case');
hXLabel = xlabel('i');
hYLabel = ylabel('j');
% hZLabel = zlabel('Number of basis functions in C_{i,j}');
xlim([0 32])
ylim([0 32])
% 赋色范围
temp = [4,8*M];
caxis(temp)
% 用colorbar命令查看赋色范围
colorbar
% view([70,29])
% view([56,54])
view(2)
% title('The number of basis functions in each cell with M=10, \delta=exp(-5)')
set(gca,'FontSize',16);