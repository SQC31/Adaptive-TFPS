function [b,b_full,Q,Q25]=solution(parameters,coefficients,bc,Msl,Msr,Msb,Mst)
%% nested functions
    function Psi0 = Psi0(x,y,i,j)
        %special solution vector at (i,j)th cell with position (x,y)
        if sigma_a(i,j)>10^-8
            % Psi0=q(i,j)*ones(4*M,1)/sigma_a(i,j);
            Psi0=(sigma_T(i,j)*eye(4*M)-(sigma_T(i,j)-sigma_a(i,j))*Kappa*diag(omega))^(-1).*q(i,j)*ones(4*M,1);
        else
            qt=q(i,j); % suppose source term is isotropic
            Psi0=-3/8*sigma_T(i,j)*qt*(x^2+y^2)+3/4*varepsilon(i,j)*qt*(ct*x+st*y);
            Psi0=Psi0-3*varepsilon(i,j)^2*qt/2/sigma_T(i,j)*(ct.*ct+st.*st)+varepsilon(i,j)^2/sigma_T(i,j)*qt;
        end
    end

%% Prepare parameters and matrices
N=parameters{1}; I=parameters{2}; J=parameters{3}; xl=parameters{4}; xr=parameters{5}; yl=parameters{6}; yr=parameters{7};
hx=(xr-xl)/I; hy=(yr-yl)/J;
[omega,~,~,M]=qnwlege2(N); % generate quadrature set
Kappa=coefficients{1}; sigma_T=coefficients{2}; sigma_a=coefficients{3}; q=coefficients{4};
psiL=bc{1}; psiR=bc{2}; psiB=bc{3}; psiT=bc{4};

% special solution
Psi0l=zeros(4*M,I(1),J(1)); Psi0r=Psi0l; Psi0b=Psi0l; Psi0t=Psi0l; Psi0c=Psi0l;
for i=1:I(1)
    for j=1:J(1)
            %special solution vectors at the left/right/bottom/top edge of each cell
            [Psi0l(:,i,j)]=Psi0((i-1)*hx+xl,(j-1/2)*hy+yl,i,j);
            [Psi0r(:,i,j)]=Psi0(i*hx+xl,(j-1/2)*hy+yl,i,j);
            [Psi0b(:,i,j)]=Psi0((i-1/2)*hx+xl,(j-1)*hy+yl,i,j);
            [Psi0t(:,i,j)]=Psi0((i-1/2)*hx+xl,j*hy+yl,i,j);
            [Psi0c(:,i,j)]=Psi0((i-1/2)*hx+xl,(j-1/2)*hy+yl,i,j);
    end
end

% for i=1:I(1)
%     for j=1:J(1)
%             %special solution vectors at the left/right/bottom/top edge of each cell
%             [Psi0l(:,i,j)]=Psi0((i-1)*hx+xl,(j-1/2)*hy+yl,i,j);
%             Psi0r(:,i,j)=Psi0l(:,i,j); 
%             Psi0b(:,i,j)=Psi0l(:,i,j); 
%             Psi0t(:,i,j)=Psi0l(:,i,j); 
%             Psi0c(:,i,j)=Psi0l(:,i,j); 
%     end
% end


%% assemble vector b
b=zeros(0,1);
count=0;
for j=1:J
    for i=1:I
        % left edge
        % cell (i,j)
        if i==1
            l=size(Msl{i,j},1);
            b(count+1:count+l,1)=Msl{i,j}*([psiL((j-1)*2*M+(M+1:2*M));zeros(2*M,1);psiL((j-1)*2*M+(1:M))]-Psi0l(:,i,j));
            count=count+l;
        else
            l=size(Msl{i,j},1);
            b(count+1:count+l,1)=Msl{i,j}*(Psi0r(:,i-1,j)-Psi0l(:,i,j));
            count=count+l;
        end
        % right edge
        if i==I
            l=size(Msr{i,j},1);
            b(count+1:count+l,1)=Msr{i,j}*([zeros(M,1);psiR((j-1)*2*M+(1:2*M));zeros(M,1)]-Psi0r(:,i,j));
            count=count+l;
        else
            l=size(Msr{i,j},1);
            b(count+1:count+l,1)=Msr{i,j}*(Psi0l(:,i+1,j)-Psi0r(:,i,j));
            count=count+l;
        end
        % bottom edge
        if j==1
            l=size(Msb{i,j},1);
            b(count+1:count+l,1)=Msb{i,j}*([psiB((i-1)*2*M+(1:2*M));zeros(2*M,1)]-Psi0b(:,i,j));
            count=count+l;
        else
            l=size(Msb{i,j},1);
            b(count+1:count+l,1)=Msb{i,j}*(Psi0t(:,i,j-1)-Psi0b(:,i,j));
            count=count+l;
        end
        % top edge
        if j==J
            l=size(Mst{i,j},1);
            b(count+1:count+l,1)=Mst{i,j}*([zeros(2*M,1);psiT((i-1)*2*M+(1:2*M))]-Psi0t(:,i,j));
            count=count+l;
        else
            l=size(Mst{i,j},1);
            b(count+1:count+l,1)=Mst{i,j}*(Psi0b(:,i,j+1)-Psi0t(:,i,j));
            count=count+l;
        end
    end
end

%% assemble vector b_full
b_full=zeros(8*M*I*J,1);
for j=1:J
    for i=1:I
        % left boundary
        if i==1
            b_full(((j-1)*I+i-1)*8*M+(1:2*M))=psiL((j-1)*2*M+(1:2*M))-[Psi0l(3*M+1:4*M,i,j);Psi0l(1:M,i,j)];
        else
            b_full(((j-1)*I+i-1)*8*M+(1:2*M))=[Psi0r(3*M+1:4*M,i-1,j);Psi0r(1:M,i-1,j)]-[Psi0l(3*M+1:4*M,i,j);Psi0l(1:M,i,j)];
        end
        % right boundary
        if i==I(1)
            b_full(((j-1)*I+i-1)*8*M+2*M+(1:2*M))=psiR((j-1)*2*M+(1:2*M))-Psi0r(M+1:3*M,i,j);
        else
            b_full(((j-1)*I+i-1)*8*M+2*M+(1:2*M))=Psi0l(M+1:3*M,i+1,j)-Psi0r(M+1:3*M,i,j);
        end
        % bottom boundary
        if j==1
            b_full(((j-1)*I+i-1)*8*M+4*M+(1:2*M))=psiB((i-1)*2*M+(1:2*M))-Psi0b(1:2*M,i,j);
        else
            b_full(((j-1)*I+i-1)*8*M+4*M+(1:2*M))=Psi0t(1:2*M,i,j-1)-Psi0b(1:2*M,i,j);
        end
        %top boundary
        if j==J
            b_full(((j-1)*I+i-1)*8*M+6*M+(1:2*M))=psiT((i-1)*2*M+(1:2*M))-Psi0t(2*M+1:4*M,i,j);
        else
            b_full(((j-1)*I+i-1)*8*M+6*M+(1:2*M))=Psi0b(2*M+1:4*M,i,j+1)-Psi0t(2*M+1:4*M,i,j);
        end
    end
end

%% assemble Q matrix
Q=zeros(4*M*I*J,1);
for j=1:J
    for i=1:I
        Q(((j-1)*I+i-1)*4*M+(1:4*M))=Psi0c(:,i,j);
    end
end

%% assemble Q25 matrix
Q25=zeros(25*4*M*I*J,1);
for j=1:J
    for i=1:I
        for y=1:5
            for x=1:5
                Q25(((5*(j-1)+y-1)*5*I+5*(i-1)+x-1)*4*M+(1:4*M))=Psi0c(:,i,j);
            end
        end
    end
end

end