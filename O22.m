clear all
clc
format long e
set(gcf,'DefaultLineLineWidth',2.5)
set(gca,'FontSize',16)

dim = 2;
I1 = 200;
I2 = 3*I1;
I3 = 3*I2;

[s1,tau1,nt1] = we_1d_fun1(I1,dim);
[s2,tau2,nt2] = we_1d_fun1(I2,dim);
[s3,tau3,nt3] = we_1d_fun1(I3,dim);

nt = nt1 - 2;
%nt = 1;

d1 = s1(:,1:nt) - s2(2:3:I2,1:3:3*nt);
d2 = s2(2:3:I2,1:3:3*nt) - s3(5:9:I3,1:9:9*nt);
dl1 = zeros(nt,1);
dl2 = dl1;
dc1 = dl1;
dc2 = dc1;
r9 = dc1;
r9(:,1) = 9;

for n = 1:nt
    dl1(n) = norm(d1(:,n),2);
    dl2(n) = norm(d2(:,n),2);
    
    dc1(n) = max(abs(d1(:,n)));
    dc2(n) = max(abs(d2(:,n)));
end

figure (1);
tt = (0:nt-1)*tau1;
plot(tt,[dl1./dl2, dc1./dc2, r9]);
title(['O2 grid convergence for 1D WE in R',num2str(dim)])
xlabel('Time') 
ylabel('Ratio')
% legend('L2-norm','C-norm','Ref=9')
legend('L2-norm','C-norm','Ref=9','Location','northwest')


% get(gcf,'default')



function [sol,tau1,nt1] = we_1d_fun1(I,d)
% solution of WE for given I cells of the grif. Task 2, grig convergence
 
% input parameters

rmin = 0;
a = 0.6;
b = 1.2;
rmax = 1.8;

cs = 1.5;  % sound speed
%d = 3;     % space dimension 

%I = 50;
CFL = 0.9;
Tmax = 1.5;

% mesh 

h = (rmax-rmin)/I;
r_i = linspace(rmin-h/2,rmax+h/2,I+2)'; 
tau = CFL*h/cs;

nt = fix(Tmax/tau);

Tmax= nt*tau;

disp([' 1D Wave equation',' in R',num2str(d)]);
ekr1 = [...
    ' Tmax =', num2str(Tmax),...
    '  tau=', num2str(tau),...
    '  CFL=', num2str(CFL)...
       ];
disp(ekr1);
disp('pause');

tau1 = tau;
nt1 = nt;
sol = zeros(I,nt);

 % pause;


% initial data

u0 = zeros(I+2,1);
u1 = zeros(I+2,1);

u0 = fun_v0(r_i,a,b);
%%%% u1 = fun_v0(r_i,a,b);
% u1 = u0 + tau*u_t + tau^2/2*u_tt + ...
% u_tt is calculated from the wave equation for v0, see below


% data preparation 

u_nm1 = u0;
% u_n = u1;  see below
u_np1 = u0(2:I+1);
d1 = (tau*cs)^2/h;
ord_i = d1./(r_i(2:I+1).^(d-1));
rdp_i = (r_i(2:I+1)+h/2).^(d-1)/h;
rdm_i = (r_i(2:I+1)-h/2).^(d-1)/h;

%% calculation of u1
u_tt = ord_i.*(rdp_i.*(u0(3:I+2)-u0(2:I+1)) - rdm_i.*(u0(2:I+1)-u0(1:I)));
u1(2:I+1) = u0(2:I+1) + 0.5*u_tt;
u_n = u1;
%%

% solution analysis preparation 

umin = 10000;
umax = -10000;

% time marching

for n = 2:nt

sol(:,n-1) = u_nm1(2:I+1); 

% interal points

u_np1 = 2*u_n(2:I+1) - u_nm1(2:I+1) + ord_i.*...
    (rdp_i.*(u_n(3:I+2)-u_n(2:I+1)) - rdm_i.*(u_n(2:I+1)-u_n(1:I))); 
u_nm1 = u_n;
u_n(2:I+1) = u_np1;

% boundary points

u_n(1) = u_n(2);
u_n(I+2) = u_n(I+1);

% sol(:,n) = u_nm1(2:I+1);

% on-line solution analysis 

% umin = min(umin,min(u_n)); 
% umax = max(umax,max(u_n));
% fekr = '%10.4e'; 
% ekran = [' time=', num2str(n*tau,fekr),...
%     ' umin=', num2str(umin,fekr),...
%     ' umax=', num2str(umax,fekr),...
%     '  ', num2str(n)];
% disp(ekran);
% plot (r_i,u_n);
% axis([rmin,rmax, -1.1,1]);  % d=1
% % axis([rmin,rmax, -2.0,4]);  % d=2
% % axis([rmin,rmax, -7.0,8]);  % d=3
% xlabel('r')
% ylabel('u_n(t)')
% pause(0.1);

end

% disp([' Tmax =', num2str(Tmax),'  tau=', num2str(tau)]);

end




function u = fun_v0(r,a,b)
% initial data at the interval [a,b] \in [min(r), max(r)]

u = r;
n = size(r);

if a >= b 
    error('a > b in initial data. Stop')
elseif a <= min(r) || b >= max(r) 
    error('check "[a,b] \in [min(r), max(r)]" in initial data. Stop')
end

for i = 1:n
  if r(i) < a || r(i) > b
      u(i) = 0;
  else
      u(i) = exp(-4*(2*r(i)-(a+b))^2/((b-a)^2-(2*r(i)-(a+b))^2));
  end
  
end
 
end

