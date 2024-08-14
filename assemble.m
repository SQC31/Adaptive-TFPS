function [A,Dc,VecSize,Fsm,Ml,Mr,Mb,Mt]=assemble(parameters,coefficients)
%% nested functions
    function [ve,va]=eigensx(i,j)
        [ve,va]=eig(diag(ct.^-1)*((1-sigma_a(i,j)/sigma_T(i,j))*Kappa*diag(omega)-eye(4*M)));
        va=diag(va);
        %Be careful! In some case, may generate complex eigenvalues
        if sum((abs(imag(va))>0))>0
            if sigma_a(i,j)<10^-12
                if sum((abs(imag(va))>0))==2
                    va(abs(imag(va))>0)=0;
                    % when sigma_a=0, 0 is a double eigenvalue of B with some eigenvectors we do not care about
                else
                    error('complex eigenvalues!');
                end
            else
                error('complex eigenvalues!');
            end
        end
        % 调整特征向量的模和方向
        [va,order]=sort(va); % sort eigenvalues
        ve=ve(:,order); % sort eigenvector matrix to be consistent with eigenvalue vector.
        len=length(va);
        for i=1:len
            ve(:,i)=ve(:,i)/max(abs(ve(:,i)));
            ve(:,i)=ve(:,i)/sign(ve(find(abs(ve(:,i))>1e-8,1),i));
        end
    end
   
    function [ve,va]=eigensy(i,j)
        [ve,va]=eig(diag(st.^-1)*((1-sigma_a(i,j)/sigma_T(i,j))*Kappa*diag(omega)-eye(4*M)));
        va=diag(va);
        if sum((abs(imag(va))>0))>0
            if sigma_a(i,j)<10^-12
                if sum((abs(imag(va))>0))==2
                    va(abs(imag(va))>0)=0;
                    % when sigma_a=0, 0 is a double eigenvalue of B with some eigenvectors we do not care about
                else
                    error('complex eigenvalues!');
                end
            else
                error('complex eigenvalues!');
            end
        end
        %Be careful! In some case, may generate complex eigenvalues
        [va,order]=sort(va); % sort eigenvalues
        ve=ve(:,order); % sort eigenvector matrix to be consistent with eigenvalue vector.
        len=length(va);
        for i=1:len
            ve(:,i)=ve(:,i)/max(abs(ve(:,i)));
            ve(:,i)=ve(:,i)/sign(ve(find(abs(ve(:,i))>1e-8,1),i));
        end
    end
   
    function [mxy]=fsm(x,y,i,j)
        % 只选取其中一部分基函数
        [vx,dx]=eigensx(i,j);
        idn=(abs(1/2*dx.*hx*sigma_T(i,j))<tol);
        dx=dx(idn); vx=vx(:,idn); nx=size(vx,2)/2;
        [vy,dy]=eigensy(i,j);
        idn=(abs(1/2*dy.*hy*sigma_T(i,j))<tol);
        dy=dy(idn); vy=vy(:,idn); ny=size(vy,2)/2;
        % 网格(i,j)中心点(x_0,y_0)
        x_0=(i-0.5)*hx+xl; y_0=(j-0.5)*hy+yl;
        % generate fundamental solution matrix when sigma_a\neq 0
        mx=zeros(4*M,2*nx); my=zeros(4*M,2*ny);
        mx(:,1:nx)=vx(:,1:nx)*diag(exp((dx(1:nx)*(x-x_0+0.5*hx))*sigma_T(i,j)));
        mx(:,nx+1:end)=vx(:,nx+1:end)*diag(exp((dx(nx+1:end)*(x-x_0-0.5*hx))*sigma_T(i,j)));
        my(:,1:ny)=vy(:,1:ny)*diag(exp((dy(1:ny)*(y-y_0+0.5*hy))*sigma_T(i,j)));
        my(:,ny+1:end)=vy(:,ny+1:end)*diag(exp((dy(ny+1:end)*(y-y_0-0.5*hy))*sigma_T(i,j)));
        mxy=[mx,my];
    end

%% Prepare parameters and matrices
N=parameters{1}; I=parameters{2}; J=parameters{3}; xl=parameters{4}; xr=parameters{5}; yl=parameters{6}; yr=parameters{7};
tol=parameters{8};
hx=(xr-xl)/I;hy=(yr-yl)/J; 
Kappa=coefficients{1}; sigma_T=coefficients{2}; sigma_a=coefficients{3};
% generate Gauss quadrature set
[omega,ct,st,M]=qnwlege2(N); % generate quadrature set

%% fsm
% fsml{i,j} maps the coefficients of Adaptive TFPS basis functions
% localized within cell(i,j) to the value of solution at the left edge
% center of cell(i,j)
% fsmr{i,j}, fsmb{i,j}, fsmt{i,j} similar
fsml=cell(I,J); fsmr=cell(I,J); fsmb=cell(I,J); fsmt=cell(I,J); fsmc=cell(I,J);
for j=1:J
    for i=1:I
        fsml{i,j}=fsm(xl+(i-1)*hx,yl+(j-1/2)*hy,i,j);
        fsmr{i,j}=fsm(xl+i*hx,yl+(j-1/2)*hy,i,j);
        fsmb{i,j}=fsm(xl+(i-1/2)*hx,yl+(j-1)*hy,i,j);
        fsmt{i,j}=fsm(xl+(i-1/2)*hx,yl+j*hy,i,j);
        fsmc{i,j}=fsm(xl+(i-1/2)*hx,yl+(j-1/2)*hy,i,j);
    end
end
Fsm={fsml,fsmr,fsmb,fsmt};
%% v2p
% Ml{i,j} maps the value of solution to the contributions to important
% velocity modes corresponding to the left edges of cell(i,j)
% Mr{i,j}, Mb{i,j}, Mt{i,j} similar
Ml=cell(I,J); Mr=Ml; Mb=Ml; Mt=Ml;
condx=zeros(I+1,J); condy=zeros(I,J+1);
for j=1:J
    % left boundary
    [Vx,Dx]=eigensx(1,j);
    idx=(abs(1/2*Dx.*hx*sigma_T(1,j))<tol); num=sum(idx)/2;
    Matrix=[Vx(3*M+1:4*M,1:2*M);Vx(1:M,1:2*M)]; idx=idx(1:2*M);
    [Q1,~] = qr(Matrix(:,idx),0); [Q2,~] = qr(Matrix(:,~idx),0);
    condx(1,j)=cond([Q1,Q2]);
    Vxi=inv([Q1,Q2]);
    Ms=zeros(num,4*M);
    Ms(:,3*M+1:4*M)=Vxi(1:num,1:M);
    Ms(:,1:M)=Vxi(1:num,M+1:2*M);
    Ml{1,j}=Ms;
    for i=1:I-1
        % between cell(i,j) and cell(i+1,j), i\ne I
        % cell(i,j)
        [Vxl,Dxl]=eigensx(i,j);
        idxl=(abs(1/2*Dxl.*hx*sigma_T(i,j))<tol); numl=sum(idxl)/2;
        % cell(i+1,j)
        [Vxr,Dxr]=eigensx(i+1,j);
        idxr=(abs(1/2*Dxr.*hx*sigma_T(i+1,j))<tol); numr=sum(idxr)/2;
        % Ms represents the projection to eigenvectors in cell(i,j)
        % Mo represents the projection to eigenvectors in other cell, cell(i+1,j)
        % Vx=[Vxl(:,2*M+1:4*M),Vxr(:,1:2*M)];
        Matrix=[Vxl(:,2*M+1:4*M),Vxr(:,1:2*M)]; idx=[idxl(2*M+1:4*M);idxr(1:2*M)];
        [Q1,~] = qr(Matrix(:,idx),0); [Q2,~] = qr(Matrix(:,~idx),0);
        condx(i+1,j)=cond([Q1,Q2]);
        Vxi=inv([Q1,Q2]);
        Mr{i,j}=Vxi(1:numl,:);
        Ml{i+1,j}=Vxi(numl+1:numl+numr,:);
    end
    % right boundary
    [Vx,Dx]=eigensx(I,j);
    idx=(abs(1/2*Dx.*hx*sigma_T(I,j))<tol); num=sum(idx)/2;
    Matrix=Vx(M+1:3*M,2*M+1:4*M); idx=idx(2*M+1:4*M);
    [Q1,~] = qr(Matrix(:,idx),0); [Q2,~] = qr(Matrix(:,~idx),0);
    condx(I+1,j)=cond([Q1,Q2]);
    Vxi=inv([Q1,Q2]);
    Ms=zeros(num,4*M);
    Ms(:,M+1:3*M)=Vxi(1:num,:);
    Mr{I,j}=Ms;
end
% assemble Mb and Mt for horizontal points
for i=1:I
    % bottom boundary
    [Vy,Dy]=eigensy(i,1);
    idy=(abs(1/2*Dy.*hy*sigma_T(i,1))<tol); num=sum(idy)/2;
    Matrix=Vy(1:2*M,1:2*M); idy=idy(1:2*M);
    [Q1,~] = qr(Matrix(:,idy),0); [Q2,~] = qr(Matrix(:,~idy),0);
    condy(i,1)=cond([Q1,Q2]);
    Vyi=inv([Q1,Q2]);
    Ms=zeros(num,4*M);
    Ms(:,1:2*M)=Vyi(1:num,:);
    Mb{i,1}=Ms;
    for j=1:J-1
        % top edge of cell(i,j), between cell(i,j) and cell(i,j+1), j\ne J
        % cell(i,j)
        [Vyb,Dyb]=eigensy(i,j);
        idyb=(abs(1/2*Dyb.*hy*sigma_T(i,j))<tol); numb=sum(idyb)/2;
        % cell(i,j+1)
        [Vyt,Dyt]=eigensy(i,j+1);
        idyt=(abs(1/2*Dyt.*hy*sigma_T(i,j+1))<tol); numt=sum(idyt)/2;
        % Ms represents the projection to eigenvectors in cell(i,j)
        % Mo represents the projection to eigenvectors in other cell, cell(i,j+1)
        % Vy=[Vyb(:,2*M+1:4*M),Vyt(:,1:2*M)];
        Matrix=[Vyb(:,2*M+1:4*M),Vyt(:,1:2*M)]; idy=[idyb(2*M+1:4*M);idyt(1:2*M)];
        [Q1,~] = qr(Matrix(:,idy),0); [Q2,~] = qr(Matrix(:,~idy),0);
        condy(i,j+1)=cond([Q1,Q2]);
        Vyi=inv([Q1,Q2]);
        Mt{i,j}=Vyi(1:numb,:);
        Mb{i,j+1}=Vyi(numb+1:numb+numt,:);
    end
    % top boundary
    [Vy,Dy]=eigensy(i,J);
    idy=(abs(1/2*Dy.*hy*sigma_T(i,J))<tol); num=sum(idy)/2;
    Matrix=Vy(2*M+1:4*M,2*M+1:4*M); idy=idy(2*M+1:4*M);
    [Q1,~] = qr(Matrix(:,idy),0); [Q2,~] = qr(Matrix(:,~idy),0);
    condy(i,J+1)=cond([Q1,Q2]);
    Vyi=inv([Q1,Q2]);
    Ms=zeros(num,4*M);
    Ms(:,2*M+1:4*M)=Vyi(1:num,:);
    Mt{i,J}=Ms;
end

VecSizel=zeros(I,J); VecSizer=VecSizel; VecSizeb=VecSizel; VecSizet=VecSizel; 
VecSize=zeros(I,J); VecSumxy=zeros(I*J+1,1);
for j=1:J
    for i=1:I
        VecSizel(i,j)=size(Ml{i,j},1); VecSizer(i,j)=size(Mr{i,j},1);
        VecSizeb(i,j)=size(Mb{i,j},1); VecSizet(i,j)=size(Mt{i,j},1);
        VecSize(i,j)=VecSizel(i,j)+VecSizer(i,j)+VecSizeb(i,j)+VecSizet(i,j);
        VecSumxy((j-1)*I+i+1)=VecSumxy((j-1)*I+i)+VecSize(i,j);
    end
end

%% assemble A matrix
% A maps the coefficients of Adaptive TFPS basis functions to the
% contribution of the basis functions to important velocity modes
list_row=zeros(0,1); list_col=zeros(0,1); list_num=zeros(0,1);
count=0;
for j=1:J
    for i=1:I
        % left edge
        % cell (i,j)
        [ti,tj,ts]=myfind(Ml{i,j}*fsml{i,j});
        l=length(ts);
        list_row(count+1:count+l)=VecSumxy((j-1)*I+i)+ti;
        list_col(count+1:count+l)=VecSumxy((j-1)*I+i)+tj;
        list_num(count+1:count+l)=ts;
        count=count+l;
        if i~=1
            % cell(i-1,j)
            [ti,tj,ts]=myfind(-Ml{i,j}*fsmr{i-1,j});
            l=length(ts);
            list_row(count+1:count+l)=VecSumxy((j-1)*I+i)+ti;
            list_col(count+1:count+l)=VecSumxy((j-1)*I+i-1)+tj;
            list_num(count+1:count+l)=ts;
            count=count+l;
        end
        % right edge
        % cell (i,j)
        [ti,tj,ts]=myfind(Mr{i,j}*fsmr{i,j});
        l=length(ts);
        list_row(count+1:count+l)=VecSumxy((j-1)*I+i)+VecSizel(i,j)+ti;
        list_col(count+1:count+l)=VecSumxy((j-1)*I+i)+tj;
        list_num(count+1:count+l)=ts;
        count=count+l;
        if i~=I
            % cell(i+1,j)
            [ti,tj,ts]=myfind(-Mr{i,j}*fsml{i+1,j});
            l=length(ts);
            list_row(count+1:count+l)=VecSumxy((j-1)*I+i)+VecSizel(i,j)+ti;
            list_col(count+1:count+l)=VecSumxy((j-1)*I+i+1)+tj;
            list_num(count+1:count+l)=ts;
            count=count+l;
        end
        % bottom edge
        % cell(i,j)
        [ti,tj,ts]=myfind(Mb{i,j}*fsmb{i,j});
        l=length(ts);
        list_row(count+1:count+l)=VecSumxy((j-1)*I+i)+VecSizel(i,j)+VecSizer(i,j)+ti;
        list_col(count+1:count+l)=VecSumxy((j-1)*I+i)+tj;
        list_num(count+1:count+l)=ts;
        count=count+l;
        if j~=1
            % cell(i,j-1)
            [ti,tj,ts]=myfind(-Mb{i,j}*fsmt{i,j-1});
            l=length(ts);
            list_row(count+1:count+l)=VecSumxy((j-1)*I+i)+VecSizel(i,j)+VecSizer(i,j)+ti;
            list_col(count+1:count+l)=VecSumxy((j-2)*I+i)+tj;
            list_num(count+1:count+l)=ts;
            count=count+l;
        end
        % top edge
        % cell(i,j)
        [ti,tj,ts]=myfind(Mt{i,j}*fsmt{i,j});
        l=length(ts);
        list_row(count+1:count+l)=VecSumxy((j-1)*I+i)+VecSizel(i,j)+VecSizer(i,j)+VecSizeb(i,j)+ti;
        list_col(count+1:count+l)=VecSumxy((j-1)*I+i)+tj;
        list_num(count+1:count+l)=ts;
        count=count+l;
        if j~=J
            % cell(i,j+1)
            [ti,tj,ts]=myfind(-Mt{i,j}*fsmb{i,j+1});
            l=length(ts);
            list_row(count+1:count+l)=VecSumxy((j-1)*I+i)+VecSizel(i,j)+VecSizer(i,j)+VecSizeb(i,j)+ti;
            list_col(count+1:count+l)=VecSumxy(j*I+i)+tj;
            list_num(count+1:count+l)=ts;
            count=count+l;
        end
    end
end
A=sparse(list_row,list_col,list_num,VecSumxy(end),VecSumxy(end));

%% assemble Dc matrix 
% Dc maps the coeffcients of basis functions to the value of solution at
% cell centers
list_row=zeros(0,1); list_col=zeros(0,1); list_num=zeros(0,1);
count=0;
for j=1:J
    for i=1:I
        block=fsmc{i,j};
        [ti,tj,ts]=myfind(block);
        l=length(ts);
        list_row(count+1:count+l)=((j-1)*I(1)+i-1)*4*M+ti;
        list_col(count+1:count+l)=VecSumxy((j-1)*I+i)+tj;
        list_num(count+1:count+l)=ts;
        count=count+l;
    end
end
Dc=sparse(list_row,list_col,list_num,4*M*I(1)*J(1),VecSumxy(end));

end