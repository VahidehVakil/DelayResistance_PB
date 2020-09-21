% This file supports paper "Dosage Strategies for Delaying Resistance Emergence in Heterogeneous Tumors"
% submitted to Proceedings of the Royal Society B: Biological Sciences
% Copyright (2020) Vahideh Vakil. All rights reserved.
% Contact: vavakil@winlab.rutgers.edu   WINLAB, Rutgers University, NJ.
% Date: 09/21/2020
%
% Subject: Containtment strategy in management phase
% The case of two cell-types, one sensetive and one resitant
% This file plots the drug response and the number of cells vs time. 
%% Initial values
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

%% drug parameters
Emax=100; EC50_0=50; EC50_1=800; 
kapa=0.0857;
kl_0=Emax/EC50_0;
kl_1=Emax/EC50_1;

D_min0=0.0025*EC50_0; %change for other options: 0.01 , 0.005, 0.001 
D_min1=D_min0;


% treatment period specification
t0=0;
T_final=30; %days length of treatment
N=4; % Number of courses
tau=t0:T_final;



% Uncomment D0 and time for optimal results
%%% These are optimal values for D_min=0.05, D_max=0.5 
D0=[0.2277 0.0500 0.0500 0.0500]; 
time=[6 12 18 24];
D1=D0;

% Uncomment D0 and time for any other results than optimal  
% D0=D_min0*ones(N,1); 
% D1=D_min1*ones(N,1); 
% time=round([1/(N+1):1/(N+1):(1-1/(N+1))]*T_final); %equal timing
%time=[4,8,22,26]; Uncomment and edit for other timing options


y=exp(time);

c0(1)=D0(1); 
c1(1)=D1(1); 
for n = 2:N
  c0(n)=D0(n)+c0(n-1)*y(n)^(-kapa)*y(n-1)^(kapa);
  c1(n)=D1(n)+c1(n-1)*y(n)^(-kapa)*y(n-1)^(kapa);
end
  
cp_t0=[]; cp_t1=[];
cp_t0=0*(0:time(1)-1); cp_t1=0*(0:time(1)-1);
for n=1:N-1
    tt=[];
    tt=[time(n):time(n+1)-1]-time(n);
    cp_t0=[cp_t0 c0(n)*exp(-kapa*tt)];
    cp_t1=[cp_t1 c1(n)*exp(-kapa*tt)];
end
tt=[];
tt=[time(N):T_final]-time(N);
cp_t0=[cp_t0 c0(N)*exp(-kapa*tt)];
cp_t1=[cp_t1 c1(N)*exp(-kapa*tt)];


h0=kl_0*cp_t0;
h1=kl_1*cp_t1;
figure(),plot(tau,h0), hold on
plot(tau,h1)
xlabel('Time (day)'); ylabel('Drug response')
h0f=h0';
h1f=h1';
tau0=tau';
tau1=tau';



%% find total number of cells at intermediate points and plot
initials=X0;
tspan=t0:T_final;
[T,celsiz]=ode45(@model2cel,tspan,initials,[],tau0,h0f,tau1,h1f,lambda0,lambda1,R0,alfa,d0,d1,u01);
figure(),plot(T,100*celsiz./sum(celsiz,2)) %percentage of the total size
xlabel('Time (day)'); ylabel('Percentage of tumor size') 
hold on

%%
function dXdt=model2cel(t,X,tau0,h0f,tau1,h1f,lambda0,lambda1,R0,alfa,d0,d1,u01)
h0f = interp1(tau0,h0f,t); % Interpolate the data set (tau0,h0f) at time t
h1f = interp1(tau1,h1f,t); % Interpolate the data set (tau1,h1f) at time t
dXdt=[lambda0*(R0-alfa*(X(1)+X(2)))*X(1)-d0*X(1)-u01*X(1)-h0f.*X(1);
    lambda1*(R0-alfa*(X(1)+X(2)))*X(2)-d1*X(2)+u01*X(1)-h1f.*X(2)];
 
 end