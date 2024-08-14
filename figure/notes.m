%% Ex1_center_N6
figure
set(gcf,'position',[250 300 600 300])
[xmesh,ymesh]=meshgrid(1/64:1/32:1-1/64,1/64:1/32:1-1/64);
subplot(1,2,1);
set(gca,'position', [0.1 0.25 0.4 0.58],'fontsize',6);
load('lattice_center_phiN6dinf.mat')
surf(xmesh,ymesh,phi');
% set(gca,'linewidth',0.1)
xlabel('x'); ylabel('y');
view(2)
shading interp
colorbar
axis square
axis tight
% Expand_axis_fill_figure(gca);

subplot(1,2,2);
set(gca,'position', [0.58 0.25 0.4 0.58],'fontsize',6);
load('lattice_center_phiN6d1.mat')
surf(xmesh,ymesh,phi');
% set(gca,'linewidth',0.1)
xlabel('x'); ylabel('y');
view(2)
shading interp
colorbar
axis square
axis tight
% Expand_axis_fill_figure(gca);

%% Ex1_error_ratio
figure
set(gcf,'position',[250 300 600 300])
delta=10.^[-1,-2,-3,-4,-5,-10];
% M=1
subplot(2,3,1);
%set(gca,'position', [0.05 0.55 0.2 0.4]);
error=[7.38964445190504e-13,7.38964445190504e-13,7.38964445190504e-13,7.38964445190504e-13,7.38964445190504e-13,7.38964445190504e-13];
ratio=[0.750000000000000,0.750000000000000,0.750000000000000,0.750000000000000,0.750000000000000,0.750000000000000];
yyaxis left; % ¼¤»î×ó±ßµÄÖá
loglog(delta,error,'-o','linewidth',1,'MarkerSize',2)
xlabel('\delta');
ylabel('error'); 
yyaxis right;
semilogx(delta,ratio,'-+','linewidth',1,'MarkerSize',2)
ylabel('ratio');
set(gca,'linewidth',0.01)
axis square
axis tight
%Expand_axis_fill_figure(gca);
set(gca,'fontsize',6)

% M=3
subplot(2,3,2);
%set(gca,'position', [0.35 0.55 0.2 0.4]);
error=[3.11453529633354e-09,3.11453529633354e-09,3.11453529633354e-09,3.11453529633354e-09,3.11453529633354e-09,3.29875016191750e-13];
ratio=[0.583333333333333,0.583333333333333,0.583333333333333,0.583333333333333,0.583333333333333,0.666666666666667];
yyaxis left; % ¼¤»î×ó±ßµÄÖá
loglog(delta,error,'-o','linewidth',1,'MarkerSize',2)
xlabel('\delta');
ylabel('error'); 
yyaxis right;
semilogx(delta,ratio,'-+','linewidth',1,'MarkerSize',2)
ylabel('ratio');
set(gca,'linewidth',0.01)
axis square
axis tight
%Expand_axis_fill_figure(gca);
set(gca,'fontsize',6)

% M=6
subplot(2,3,3);
%set(gca,'position', [0.65 0.55 0.2 0.4]);
error=[9.87362669402359e-09,9.87362669402359e-09,9.87362669402359e-09,9.87362669402359e-09,9.87362669402359e-09,4.54747350886464e-13];
ratio=[0.541666666666667,0.541666666666667,0.541666666666667,0.541666666666667,0.541666666666667,0.750000000000000];
yyaxis left; % ¼¤»î×ó±ßµÄÖá
loglog(delta,error,'-o','linewidth',1,'MarkerSize',2)
xlabel('\delta');
ylabel('error'); 
yyaxis right;
semilogx(delta,ratio,'-+','linewidth',1,'MarkerSize',2)
ylabel('ratio');
set(gca,'linewidth',0.01)
axis square
axis tight
%Expand_axis_fill_figure(gca);
set(gca,'fontsize',6)

% M=10
subplot(2,3,4);
%set(gca,'position', [0.05 0.1 0.2 0.4]);
error=[1.43488388948398e-08,1.43488388948398e-08,1.43488388948398e-08,1.43488388948398e-08,1.43488388948398e-08,1.26476606965298e-12];
ratio=[0.525000000000000,0.525000000000000,0.525000000000000,0.525000000000000,0.525000000000000,0.675000000000000];
yyaxis left; % ¼¤»î×ó±ßµÄÖá
loglog(delta,error,'-o','linewidth',1,'MarkerSize',2)
xlabel('\delta');
ylabel('error'); 
yyaxis right;
semilogx(delta,ratio,'-+','linewidth',1,'MarkerSize',2)
ylabel('ratio');
set(gca,'linewidth',0.01)
axis square
axis tight
%Expand_axis_fill_figure(gca);
set(gca,'fontsize',6)

% M=15
subplot(2,3,5);
%set(gca,'position', [0.35 0.1 0.2 0.4]);
error=[1.62733784225821e-08,1.62733784225821e-08,1.62733784225821e-08,1.62733784225821e-08,1.62733784225821e-08,1.22213350550737e-12];
ratio=[0.516666666666667,0.516666666666667,0.516666666666667,0.516666666666667,0.516666666666667,0.700000000000000];
yyaxis left; % ¼¤»î×ó±ßµÄÖá
loglog(delta,error,'-o','linewidth',1,'MarkerSize',2)
xlabel('\delta');
ylabel('error'); 
yyaxis right;
semilogx(delta,ratio,'-+','linewidth',1,'MarkerSize',2)
ylabel('ratio');
set(gca,'linewidth',0.01)
axis square
axis tight
%Expand_axis_fill_figure(gca);
set(gca,'fontsize',6)

% M=21
subplot(2,3,6);
%set(gca,'position', [0.65 0.1 0.2 0.4]);
error=[1.66835193438075e-08,1.66835193438075e-08,1.66835193438075e-08,1.66835193438075e-08,1.66835193438075e-08,1.75111036782027e-11];
ratio=[0.511904761904762,0.511904761904762,0.511904761904762,0.511904761904762,0.511904761904762,0.654761904761905];
yyaxis left; % ¼¤»î×ó±ßµÄÖá
loglog(delta,error,'-o','linewidth',1,'MarkerSize',2)
xlabel('\delta');
ylabel('error'); 
yyaxis right;
semilogx(delta,ratio,'-+','linewidth',1,'MarkerSize',2)
ylabel('ratio');
set(gca,'linewidth',0.01)
axis square
axis tight
%Expand_axis_fill_figure(gca);
set(gca,'fontsize',6)

%% Ex1_layer
figure
set(gcf,'position',[250 300 600 200])
[xmesh,ymesh]=meshgrid(1/64:1/32:1-1/64,1/64:1/32:1-1/64);
xmesh=1/32+1/32/20:1/32/10:7/32-1/32/20;
load('lattice_layer_phiN6dinf.mat')
% delta=10^{-1}
subplot(1,3,1);
load('lattice_layer_phiN6d1.mat')
temp1=phi_inf-phi;
plot(xmesh,temp1,'linewidth',0.5);
set(gca,'linewidth',0.01)
axis square
axis tight
% Expand_axis_fill_figure(gca);
set(gca,'fontsize',6)

% delta=10^{-10}
subplot(1,3,2);
load('lattice_layer_phiN6d10.mat')
temp1=phi_inf-phi;
plot(xmesh,temp1,'linewidth',0.5);
set(gca,'linewidth',0.01)
axis square
axis tight
% Expand_axis_fill_figure(gca);
set(gca,'fontsize',6)

% delta=10^{-20}
subplot(1,3,3);
load('lattice_layer_phiN6d20.mat')
temp1=phi_inf-phi;
plot(xmesh,temp1,'linewidth',0.5);
set(gca,'linewidth',0.01)
axis square
axis tight
% Expand_axis_fill_figure(gca);
set(gca,'fontsize',6)
%% Ex2_vecsize
figure
set(gcf,'position',[250 300 600 300])

subplot(2,3,1);
set(gca,'position', [0.07 0.55 0.3 0.4]);
colormap default
load('buffer_vecsize_N2d1.mat')
b = bar3(VecSize,1);
for n=1:numel(b)
   cdata=get(b(n),'zdata');
   cdata=repmat(max(cdata,[],2),1,4);
   set(b(n),'cdata',cdata,'facecolor','flat')
end
xlabel('i'); ylabel('j');
xlim([0 32])
ylim([0 32])
temp = [4,24];
caxis(temp)
colorbar
view(2)
set(gca,'linewidth',0.01)
axis square
axis tight
% Expand_axis_fill_figure(gca);
set(gca,'fontsize',6)


subplot(2,3,2);
set(gca,'position', [0.39 0.55 0.3 0.4]);
colormap default
load('buffer_vecsize_N2d3.mat')
b = bar3(VecSize,1);
for n=1:numel(b)
   cdata=get(b(n),'zdata');
   cdata=repmat(max(cdata,[],2),1,4);
   set(b(n),'cdata',cdata,'facecolor','flat')
end
xlabel('i'); ylabel('j');
xlim([0 32])
ylim([0 32])
temp = [4,24];
caxis(temp)
colorbar
view(2)
set(gca,'linewidth',0.01)
axis square
axis tight
% Expand_axis_fill_figure(gca);
set(gca,'fontsize',6)

subplot(2,3,3);
set(gca,'position', [0.71 0.55 0.3 0.4]);
colormap default
load('buffer_vecsize_N2d5.mat')
b = bar3(VecSize,1);
for n=1:numel(b)
   cdata=get(b(n),'zdata');
   cdata=repmat(max(cdata,[],2),1,4);
   set(b(n),'cdata',cdata,'facecolor','flat')
end
xlabel('i'); ylabel('j');
xlim([0 32])
ylim([0 32])
temp = [4,24];
caxis(temp)
colorbar
view(2)
set(gca,'linewidth',0.01)
axis square
axis tight
% Expand_axis_fill_figure(gca);
set(gca,'fontsize',6)

subplot(2,3,4);
set(gca,'position', [0.07 0.05 0.3 0.4]);
colormap default
load('buffer_vecsize_N6d1.mat')
b = bar3(VecSize,1);
for n=1:numel(b)
   cdata=get(b(n),'zdata');
   cdata=repmat(max(cdata,[],2),1,4);
   set(b(n),'cdata',cdata,'facecolor','flat')
end
xlabel('i'); ylabel('j');
xlim([0 32])
ylim([0 32])
temp = [4,168];
caxis(temp)
colorbar
view(2)
set(gca,'linewidth',0.01)
axis square
axis tight
% Expand_axis_fill_figure(gca);
set(gca,'fontsize',6)

subplot(2,3,5);
set(gca,'position', [0.39 0.05 0.3 0.4]);
colormap default
load('buffer_vecsize_N6d3.mat')
b = bar3(VecSize,1);
for n=1:numel(b)
   cdata=get(b(n),'zdata');
   cdata=repmat(max(cdata,[],2),1,4);
   set(b(n),'cdata',cdata,'facecolor','flat')
end
xlabel('i'); ylabel('j');
xlim([0 32])
ylim([0 32])
temp = [4,168];
caxis(temp)
colorbar
view(2)
set(gca,'linewidth',0.01)
axis square
axis tight
% Expand_axis_fill_figure(gca);
set(gca,'fontsize',6)

subplot(2,3,6);
set(gca,'position', [0.71 0.05 0.3 0.4]);
colormap default
load('buffer_vecsize_N6d5.mat')
b = bar3(VecSize,1);
for n=1:numel(b)
   cdata=get(b(n),'zdata');
   cdata=repmat(max(cdata,[],2),1,4);
   set(b(n),'cdata',cdata,'facecolor','flat')
end
xlabel('i'); ylabel('j');
xlim([0 32])
ylim([0 32])
temp = [4,168];
caxis(temp)
colorbar
view(2)
set(gca,'linewidth',0.01)
axis square
axis tight
% Expand_axis_fill_figure(gca);
set(gca,'fontsize',6)


%% Ex2_center_N6
figure
set(gcf,'position',[250 300 600 300])
[xmesh,ymesh]=meshgrid(1/64:1/32:1-1/64,1/64:1/32:1-1/64);
subplot(1,2,1);
set(gca,'position', [0.1 0.25 0.4 0.58],'fontsize',6);
load('buffer_center_phiN6dinf.mat')
surf(xmesh,ymesh,phi');
% set(gca,'linewidth',0.1)
xlabel('x'); ylabel('y');
view(2)
shading interp
colorbar
axis square
axis tight
% Expand_axis_fill_figure(gca);

subplot(1,2,2);
set(gca,'position', [0.58 0.25 0.4 0.58],'fontsize',6);
load('buffer_center_phiN6d1.mat')
surf(xmesh,ymesh,phi');
% set(gca,'linewidth',0.1)
xlabel('x'); ylabel('y');
view(2)
shading interp
colorbar
axis square
axis tight
% Expand_axis_fill_figure(gca);


%% Ex2_error_ratio
figure
set(gcf,'position',[250 300 600 300])
delta=10.^[-1,-2,-3,-4,-5,-10];
% M=1
subplot(2,3,1);
%set(gca,'position', [0.05 0.55 0.2 0.4]);
error=[4.24416925874360e-05,1.58569349484328e-06,1.48057726551132e-07,9.19563503121168e-09,1.23443666399936e-09,5.80438475061840e-15];
ratio=[0.500000000000000,0.819824218750000,0.895996093750000,0.930664062500000,0.950683593750000,0.988281250000000];
yyaxis left; % ¼¤»î×ó±ßµÄÖá
loglog(delta,error,'-o','linewidth',1,'MarkerSize',2)
xlabel('\delta');
ylabel('error'); 
yyaxis right;
semilogx(delta,ratio,'-+','linewidth',1,'MarkerSize',2)
ylabel('ratio');
set(gca,'linewidth',0.01)
axis square
axis tight
%Expand_axis_fill_figure(gca);
set(gca,'fontsize',6)

% M=3
subplot(2,3,2);
%set(gca,'position', [0.35 0.55 0.2 0.4]);
error=[7.73671936797338e-05,1.63905538917986e-05,5.05773921188890e-07,2.37236734423885e-08,1.71630110123999e-09,5.54417622922188e-15];
ratio=[0.226888020833333,0.468750000000000,0.721679687500000,0.816406250000000,0.865071614583333,0.954101562500000];
yyaxis left; % ¼¤»î×ó±ßµÄÖá
loglog(delta,error,'-o','linewidth',1,'MarkerSize',2)
xlabel('\delta');
ylabel('error'); 
yyaxis right;
semilogx(delta,ratio,'-+','linewidth',1,'MarkerSize',2)
ylabel('ratio');
set(gca,'linewidth',0.01)
axis square
axis tight
%Expand_axis_fill_figure(gca);
set(gca,'fontsize',6)

% M=6
subplot(2,3,3);
%set(gca,'position', [0.65 0.55 0.2 0.4]);
error=[0.000117254818158996,2.30049695189860e-05,2.21703497027548e-06,1.45528142703188e-07,3.90306773012661e-09,1.74998904256540e-14];
ratio=[0.171549479166667,0.398030598958333,0.591878255208333,0.745442708333333,0.817708333333333,0.936930338541667];
yyaxis left; % ¼¤»î×ó±ßµÄÖá
loglog(delta,error,'-o','linewidth',1,'MarkerSize',2)
xlabel('\delta');
ylabel('error'); 
yyaxis right;
semilogx(delta,ratio,'-+','linewidth',1,'MarkerSize',2)
ylabel('ratio');
set(gca,'linewidth',0.01)
axis square
axis tight
%Expand_axis_fill_figure(gca);
set(gca,'fontsize',6)

% M=10
subplot(2,3,4);
%set(gca,'position', [0.05 0.1 0.2 0.4]);
error=[0.000101900192052251,6.35693729911235e-06,3.05173474668408e-06,2.76114333952471e-07,1.64280729020305e-08,3.28903571045203e-14];
ratio=[0.142626953125000,0.428710937500000,0.529394531250000,0.675683593750000,0.776416015625000,0.924072265625000];
yyaxis left; % ¼¤»î×ó±ßµÄÖá
loglog(delta,error,'-o','linewidth',1,'MarkerSize',2)
xlabel('\delta');
ylabel('error'); 
yyaxis right;
semilogx(delta,ratio,'-+','linewidth',1,'MarkerSize',2)
ylabel('ratio');
set(gca,'linewidth',0.01)
axis square
axis tight
%Expand_axis_fill_figure(gca);
set(gca,'fontsize',6)

% M=15
subplot(2,3,5);
%set(gca,'position', [0.35 0.1 0.2 0.4]);
error=[0.000109605645031596,1.24899206910278e-05,1.53003749236991e-06,2.81363394194489e-07,2.46860158426010e-08,4.32154312335342e-14];
ratio=[0.136393229166667,0.438476562500000,0.540592447916667,0.627766927083333,0.732682291666667,0.912500000000000];
yyaxis left; % ¼¤»î×ó±ßµÄÖá
loglog(delta,error,'-o','linewidth',1,'MarkerSize',2)
xlabel('\delta');
ylabel('error'); 
yyaxis right;
semilogx(delta,ratio,'-+','linewidth',1,'MarkerSize',2)
ylabel('ratio');
set(gca,'linewidth',0.01)
axis square
axis tight
%Expand_axis_fill_figure(gca);
set(gca,'fontsize',6)

% M=21
subplot(2,3,6);
%set(gca,'position', [0.65 0.1 0.2 0.4]);
error=[0.000110258758869886,1.51680994716075e-05,9.51237750179690e-07,3.13269000462579e-07,2.79712495387940e-08,5.27911048209262e-14];
ratio=[0.125348772321429,0.434175037202381,0.563360305059524,0.620047433035714,0.699288504464286,0.903994605654762];
yyaxis left; % ¼¤»î×ó±ßµÄÖá
loglog(delta,error,'-o','linewidth',1,'MarkerSize',2)
xlabel('\delta');
ylabel('error'); 
yyaxis right;
semilogx(delta,ratio,'-+','linewidth',1,'MarkerSize',2)
ylabel('ratio');
set(gca,'linewidth',0.01)
axis square
axis tight
%Expand_axis_fill_figure(gca);
set(gca,'fontsize',6)

%% boundedness of Crho
figure
set(gcf,'position',[200 100 450 600])
x=exp(-25:0.5:-0.5);
y=exp(-25:0.5:-0.5);
I=length(x); 
xmesh=ones(I,1)*x; ymesh=y'*ones(1,I);
% g_{C_{-}}=g_{C_{+}}=0
load('crho00.mat')
for N=2:5
    subplot(4,3,3*(N-2)+1);
    surf(xmesh,ymesh,crho(:,:,N));
    xlim([0,1]); ylim([0,1]);
    set(gca,'xtick',[10^(-10) 10^(-5) 1],'xticklabel',{'1-10^{-10}','1-10^{-5}','0'});
    set(gca,'ytick',[10^(-10) 10^(-5) 1],'yticklabel',{'1-10^{-10}','1-10^{-5}','0'});
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    shading interp 
    view([55,29])
    pos=axis;
    xlabel('\gamma_{C-}')
    ylabel('\gamma_{C+}')
    set(gca,'linewidth',0.01)
    % axis square
    % axis tight
    set(gca,'FontSize',6);
end

% g_{C_{-}}=0.3, g_{C_{+}}=0
load('crho30.mat')
for N=2:5
    subplot(4,3,3*(N-2)+2);
    surf(xmesh,ymesh,crho(:,:,N));
    xlim([0,1]); ylim([0,1]);
    set(gca,'xtick',[10^(-10) 10^(-5) 1],'xticklabel',{'1-10^{-10}','1-10^{-5}','0'});
    set(gca,'ytick',[10^(-10) 10^(-5) 1],'yticklabel',{'1-10^{-10}','1-10^{-5}','0'});
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    shading interp 
    view([55,29])
    xlabel('\gamma_{C-}')
    ylabel('\gamma_{C+}')
    set(gca,'linewidth',0.01)
    % axis square
    % axis tight
    set(gca,'FontSize',6);
end

% g_{C_{-}}=0.2, g_{C_{+}}=-0.3
load('crho23.mat')
for N=2:5
    subplot(4,3,3*(N-2)+3);
    surf(xmesh,ymesh,crho(:,:,N));
    xlim([0,1]); ylim([0,1]);
    set(gca,'xtick',[10^(-10) 10^(-5) 1],'xticklabel',{'1-10^{-10}','1-10^{-5}','0'});
    set(gca,'ytick',[10^(-10) 10^(-5) 1],'yticklabel',{'1-10^{-10}','1-10^{-5}','0'});
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    shading interp 
    view([55,29])
    xlabel('\gamma_{C-}')
    ylabel('\gamma_{C+}')
    set(gca,'linewidth',0.01)
    % axis square
    % axis tight
    set(gca,'FontSize',6);
end

%% linear independence of eigenvectors at interior edge center
% rank v.s. \gamma for different values of g and M
for N=2:5
g1=0.2; g2=-0.3;
% g1=0.3; g2=0;
% g1=0; g2=0;
[omega,ct,st,M,the]=qnwlege2(N);
[ Kappa1, val1 ] = HG2( N, g1 );
val1;
[ Kappa2, val2 ] = HG2( N, g2 );
val2;
x=exp(-25:0.1:-0.1);
y=exp(-25:0.1:-0.1);
I=length(x); 
xmesh=ones(I,1)*x; ymesh=y'*ones(1,I); zmesh=zeros(I,I);
for j=1:I
    for i=1:I
        % left cell
        [vel,val]=eig(diag(ct.^-1)*((1-x(i))*Kappa1*diag(omega)-eye(4*M)));
        val=diag(val);
        [val,orderl]=sort(val); % sort eigenvalues
        vel=vel(:,orderl); % sort eigenvector matrix to be consistent with eigenvalue vector.
        % right cell
        [ver,var]=eig(diag(ct.^-1)*((1-y(j))*Kappa2*diag(omega)-eye(4*M)));
        var=diag(var);
        [var,orderr]=sort(var); % sort eigenvalues
        ver=ver(:,orderr); % sort eigenvector matrix to be consistent with eigenvalue vector.
        E=[vel(:,2*M+1:end),ver(:,1:2*M)];
        zmesh(i,j)=rank(E)/(4*M);
    end
end
figure(N)
set(gcf,'position',[200,200,300,300]);
surf(xmesh,ymesh,zmesh);
shading interp 
set(gca,'xtick',[10^(-10) 10^(-5) 1],'xticklabel',{'1-10^{-10}','1-10^{-5}','0'});
set(gca,'ytick',[10^(-10) 10^(-5) 1],'yticklabel',{'1-10^{-10}','1-10^{-5}','0'});
xlabel('\gamma_{C-}')
ylabel('\gamma_{C+}')
set(gca,'XScale','log') 
set(gca,'YScale','log') 
xlim([0,1]); ylim([0,1]);
title(['rank ratio when g_{C-}\equiv',num2str(g1),' and g_{C+}\equiv',num2str(g2)])
title('rank ratio')
set(gca,'FontSize',12);
end

%% linear independence of eigenvectors at left boundary edge center
% rank v.s. \gamma for different values of g and M
for N=2:5
g=-0.3;
% g=0.3;
[omega,ct,st,M,the]=qnwlege2(N);
[ Kappa, va ] = HG2( N, g );
xmesh=exp(-25:0.1:-0.1);
I=length(xmesh); 
ymesh=zeros(1,I);
for i=1:I
    [ve,va]=eig(diag(ct.^-1)*((1-xmesh(i))*Kappa*diag(omega)-eye(4*M)));
    va=diag(va);
    [va,order]=sort(va); % sort eigenvalues
    ve=ve(:,order); % sort eigenvector matrix to be consistent with eigenvalue vector.
    E=[ve(3*M+1:4*M,1:2*M);ve(1:M,1:2*M)];
    ymesh(i)=rank(E)/(2*M);
end
figure(N)
set(gcf,'position',[200,200,300,300]);
semilogx(xmesh,ymesh)
set(gca,'xtick',[10^(-10) 10^(-7.5) 10^(-5) 10^(-2.5) 1],'xticklabel',{'1-10^{-10}','1-10^{-7.5}','1-10^{-5}','1-10^{-2.5}','0'});
xlabel('\gamma_{C}')
xlim([0,1])
title(['rank ratio when M=',num2str(M),' and g=',num2str(g)])
title('rank ratio')
set(gca,'FontSize',12);
end