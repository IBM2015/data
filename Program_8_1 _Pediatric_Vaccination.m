function [T,S,I,R] = Program_8_1_Pediatric_Vaccination(beta,gamma,mu,S0,I0,p, tV, MaxTime)
%
% 
%
% RISK_STRUCTURE( beta, gamma, mu, S0, I0, p, tV, MaxTime)
%      This is the MATLAB version of program 8.1 from page 293 of 
% "Modeling Infectious Disease in humans and animals" 
% by Keeling & Rohani.
% 
% It is the SIR epidemic (with equal births and deaths) and prophylactic 
% vaccination. Vaccination starts at time tV, after which a proportion p of
% all newborne individuals are vaccinated
%

% Sets up default parameters if necessary.
if nargin == 0
   beta=520/365.0;
   gamma=1/7.0;
   mu=1/(70*365.0);
   S0=0.1;
   I0=1e-4;
   R0=0;
   p=0.7;
   tV=30*365;
   MaxTime=100*365;
end

% Checks all the parameters are valid
if S0<=0 
    error('Initial level of susceptibles (%g) is less than or equal to zero',S0);
end

if I0<=0 
    error('Initial level of infecteds (%g) is less than or equal to zero',I0);
end

if R0<0 
    error('Initial level of Removable (%g) is less than or equal to zero',R0);
end

if beta<=0 
    error('Transmission rate beta (%g) is less than or equal to zero',beta);
end
    
if gamma<=0 
    error('Recovery rate gamma (%g) is less than or equal to zero',gamma);
end

if mu<0 
    error('Birth / Death rate mu (%g) is less than zero',mu);
end

if p<0 || p>1 
    error('Vaccinated proportion p (%g) is not between zero and one',p);
end
 

if tV<0 || tV>=MaxTime
    error('Onset of vaccination, tV (%g), is not between zero and the end time',tV);
end
    
if MaxTime<=0 
    error('Maximum run time (%g) is less than or equal to zero',MaxTime);
end
    
if S0+I0>1
    warning('Initial level of susceptibles+infecteds (%g+%g=%g) is greater than one',S0,I0,S0+I0);
end

if beta<gamma+mu
    warning('Basic reproductive ratio (R_0=%g) is less than one',beta/(gamma+mu));
end

S=S0; I=I0; R=(S0+I0+R0)-S-I;

% The main iteration 
options = odeset('RelTol', 1e-5, 'AbsTol', 1e-4);
[t, pop]=ode45(@Diff_8_1,[0 tV],[S I R],options,[beta gamma mu 0]);
T=t; S=pop(:,1); I=pop(:,2); R=pop(:,3);

[t, pop]=ode45(@Diff_8_1,[tV MaxTime],pop(end,:),options,[beta gamma mu p]);
T=[T; t]; S=[S; pop(:,1)]; I=[I; pop(:,2)]; R=[R; pop(:,3)];

% plots the graphs with scaled colours
subplot(3,1,1)
h=plot(T,S,'-g'); 
hold on; Y=get(gca,'YLim'); plot(tV+[0 0],Y,'--k'); hold off
xlabel 'Time';
ylabel 'Susceptibles'

subplot(3,1,2) 
h=plot(T,I,'-r');
hold on; Y=get(gca,'YLim'); plot(tV+[0 0],Y,'--k'); hold off
xlabel 'Time';
ylabel 'Infectious'

subplot(3,1,3)
h=plot(T,R,'-k'); 
hold on; Y=get(gca,'YLim'); plot(tV+[0 0],Y,'--k'); hold off
xlabel 'Time';
ylabel 'Recovereds/Vaccinated'


% Calculates the differential rates used in the integration.
function dPop=Diff_8_1(t,pop, parameter)

beta=parameter(1); gamma=parameter(2); mu=parameter(3);  p=parameter(4);
S=pop(1); I=pop(2); R=pop(3);

dPop=zeros(3,1);

dPop(1)= mu*(1-p) -beta*S*I - mu*S;
dPop(2)= beta*S*I - gamma*I - mu*I;
dPop(3)= mu*p + gamma*I - mu*R;

