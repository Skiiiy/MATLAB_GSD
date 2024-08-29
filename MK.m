function [K] = MK(M)
%MK 二阶A-M关系K参数计算
gamma=1.4;
rhoc=1650;W=1;
rc=(3*W/(4*pi*rhoc))^(1/3);
Zc=rc/(W^(1/3));
zeta=2*gamma.*(M.^2-1)./(gamma+1);
syms Z
eqn=zeta-(808.*(1+(Z./(4.5)).^2))./sqrt(((1+(Z./0.048).^2).*(1+(Z./0.32).^2).*(1+(Z./1.35).^2)))==0;
Z=vpasolve(eqn,Z);
Z=eval(Z);
% KZ=(808.*(1+(Z./(4.5)).^2))./sqrt(((1+(Z./0.048).^2).*(1+(Z./0.32).^2).*(1+(Z./1.35).^2)));
dKZ_KZ=2./(Z.*(1+(4.5./Z).^2))-1./(Z.*(1+(0.048./Z).^2))-1./(Z.*(1+(0.32./Z).^2))-1./(Z.*(1+(1.35./Z).^2));

mu=sqrt(((gamma-1).*(M.^2)+2)./(2.*gamma.*(M.^2)+1-gamma));
lambda=(1+(2/(gamma+1)).*((1-mu.^2)./mu)).*(1+2.*mu+1./(M.^2));

Exiang=-(2.*Zc)./Z-(lambda./2).*dKZ_KZ.*Zc;
Axiang=2.*Zc./Z;
K=Exiang./Axiang;
end

