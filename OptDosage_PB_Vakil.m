% This file supports paper "Dosage Strategies for Delaying Resistance Emergence in Heterogeneous Tumors"
% submitted to Proceedings of the Royal Society B: Biological Sciences
% Copyright (2020) Vahideh Vakil. All rights reserved.
% Contact: vavakil@winlab.rutgers.edu  WINLAB, Rutgers University, NJ.
% Date: 09/21/2020
%
% Subject: Containtment strategy in management phase
% The case of two cell-types, one sensetive and one resitant
% This file finds the optimal dosage that maximizes the time of resistance.
%% Initial values and model parameters
clear all
lambda0=0.185; 
d0=0.02;
lambda1=0.185; 
d1=0.12;
u01=0.001; 
x0=1e3; x1=1e1; % These are the initial values at the begining of the management phase
X0=[x0;x1];
Smax=sum(X0);

% Nutrient source parameters
alfa=0.9;
R0=1500; 

%global model_param 
model_param=[lambda0 d0 lambda1 d1 u01 alfa R0];  
%% drug parameters
Emax=100; EC50_0=50; EC50_1=800; 
kapa=0.0857;
kl_0=Emax/EC50_0;
kl_1=Emax/EC50_1;

D_min=0.001*EC50_0; 
D_max=0.01*EC50_0;


% treatment period specification
T_final=30; %days length of treatment
N=4; % Number of courses
tau_min=4; %days dosing interval

%global drug_param
drug_param=[N T_final kapa kl_0 kl_1];

D0=D_min*ones(N,1); 
time0=round([1/(N+1):1/(N+1):(1-1/(N+1))]*T_final); %equal timing
dosage0=[D0;time0']; % Initial dosage for optimization 

lb=zeros(2*N,1); ub=zeros(2*N,1);
for n=1:N
    lb(n)=D_min; ub(n)=D_max;
end
for n=N+1:2*N
    lb(n)=(n-N)*tau_min; ub(n)=T_final-(2*N-n+1)*tau_min;
end


 
A=[0 0 0 0 1 -1 0 0;
   0 0 0 0 0 1 -1 0;
   0 0 0 0 0 0 1 -1];
b=-tau_min*ones(3,1);

%% Optimization


[dosage_star,t_star] =  runnested(drug_param,model_param,X0,dosage0,A,b,lb,ub)

function [dosage_star,t_star] =  runnested(drug_param,model_param,initials,dosage0,A,b,lb,ub)
[dosage_star,t_star]=patternsearch(@objfun,dosage0,A,b,[],[],lb,ub,[],[]);


% Objective function
function f=objfun(dosage)

T_final=drug_param(2);
tspan=0:T_final;

sol=ode45(@model2cel,tspan,initials,[],dosage,drug_param,model_param);

if (deval(sol,0,2)-deval(sol,0,1))*(deval(sol,T_final,2)-deval(sol,T_final,1))<0
equaltim=fzero(@(r)(deval(sol,r,2)-deval(sol,r,1)),[0 T_final]);
else
equaltim=0;
end

f = -round(equaltim); % take negative of time for maximization
end
end


%% functions

% Set ODE system
function dXdt=model2cel(t,X,dosage,drug_param,model_param)
%D0=dosage(1:N); time=dosage(N+1:end);
T_final=drug_param(2);
tau=0:T_final;
tau0=tau';tau1=tau';
%model_param=[lambda0 d0 lambda1 d1 u01 alfa R0];
lambda0=model_param(1); d0=model_param(2);
lambda1=model_param(3); d1=model_param(4); u01=model_param(5);
alfa=model_param(6); R0=model_param(7);
Hf=drugRes(dosage,drug_param);
h0f = interp1(tau0,Hf(:,1),t); % Interpolate the data set (tau0,h0f) at time t
h1f = interp1(tau1,Hf(:,2),t); % Interpolate the data set (tau1,h1f) at time t
dXdt=[lambda0*(R0-alfa*(X(1)+X(2)))*X(1)-d0*X(1)-u01*X(1)-h0f.*X(1);
    lambda1*(R0-alfa*(X(1)+X(2)))*X(2)-d1*X(2)+u01*X(1)-h1f.*X(2)];
 
end
 

% Find drug response 
function Hf=drugRes(dosage,drug_param)
%drug_param=[N T_final kapa kl_0 kl_1];

N=drug_param(1); T_final=drug_param(2);
kapa=drug_param(3); kl_0=drug_param(4); kl_1=drug_param(5);
D0=dosage(1:N); time=round(dosage(N+1:end));

y=exp(time);
c0(1)=D0(1); 

for n = 2:N
  c0(n)=D0(n)+c0(n-1)*y(n)^(-kapa)*y(n-1)^(kapa);
  
end
  
cp_t0=[]; cp_t1=[];
cp_t0=0*(0:time(1)-1); 
for n=1:N-1
    tt=[];
    tt=[time(n):time(n+1)-1]-time(n);
    cp_t0=[cp_t0 c0(n)*exp(-kapa*tt)];
    
end
tt=[];
tt=[time(N):T_final]-time(N);
cp_t0=[cp_t0 c0(N)*exp(-kapa*tt)];
cp_t1=cp_t0;

h0=kl_0*cp_t0;
h1=kl_1*cp_t1;

Hf(:,1)=h0';
Hf(:,2)=h1';

end

