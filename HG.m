function [ Kappa, relres ] = HG( N, g )
%HG compute discrete HG kernel
% Input: ordinate parameter N, and anisotropy parameter g
% kmn∞¥––≈≈¡–
[omega,ct,st,zeta,M]=qnwlege3(N);
u=[ct';st';zeta'];
f_kappa=@(the) (1).*(1-g^2)./(1+g^2-2*g.*cos(the)).^(3/2);
kappa=f_kappa(ct*ct'+st*st'+zeta*zeta');
kappa=kappa(1:4*M,:)';
kappa=kappa(:);
H=eye(4*M*8*M); f=kappa;
E1=zeros(4*M,4*M*8*M); E2=zeros(4*M,4*M*8*M); E3=zeros(12*M,4*M*8*M);
% b1=ones(4*M-3,1); b2=ones(4*M,1); b3=g*reshape(u(:,1:4*M),12*M,1);
b1=ones(4*M,1); b2=ones(4*M,1); b3=g*reshape(u(:,1:4*M),12*M,1);
for m=1:4*M
    for n=1:4*M
        E1(m,(n-1)*8*M+m)=omega(n);
        E1(m,(n-1)*8*M+m+4*M)=omega(n);
    end
end

for n=1:4*M
    for m=1:8*M
        E2(n,(n-1)*8*M+m)=omega(m);
    end
end

for n=1:4*M
    for m=1:8*M
        E3((n-1)*3+(1:3),(n-1)*8*M+m)=omega(m)*u(:,m);
    end
end
% E1=[E1(1:5,:);E1(7:15,:);E1(17:23,:);E1(25:end,:)];
Aeq=[E1;E2;E3]; beq=[b1;b2;b3]; lb=zeros(4*M*8*M,1);
% Aeq=[E1;E2;E4;E5]; beq=[b1;b2;b4;b5]; lb=zeros(8*M*8*M,1);
% options = optimoptions('quadprog','MaxIterations',200);
options = optimoptions('lsqlin','Algorithm','interior-point');
% [x,fval] = quadprog(H,f,[],[],Aeq,beq,lb,[],[],options);
% [x,resnorm,residual,exitflag,output,lambda] = lsqlin(H,f,[],[],Aeq,beq,lb,[],[],options);
[x,resnorm,residual,exitflag,output,lambda] = lsqlin(H,f,[],[],Aeq,beq,lb,[],[],options);
% [x,fval] = quadprog(H,f,[],[],Aeq,beq,lb,[],[],options);
% val=sqrt(2*(fval+1/2*(kappa'*kappa))/size(x,1));
relres=sqrt(resnorm/(4*M*8*M));

Kappa=zeros(4*M,4*M);
% for n=1:4*M
%     for m=1:4*M
%         Kappa((n-1)*4*M+m)=1/2*(x((n-1)*8*M+m)+x((n+4*M-1)*8*M+m));
%     end
% end

% Kappa=transpose(reshape(Kappa,4*M,4*M));
for n=1:4*M
    for m=1:4*M
        Kappa(n,m)=1/2*(x((n-1)*8*M+m)+x((n-1)*8*M+m+4*M));
    end
end
omega=2.*omega(1:4*M);
end

