function [t,S,I] = Program_2_5(beta,gamma,I0,MaxTime)
%
% 
%
% RISK_STRUCTURE( beta, gamma, I0, MaxTime)
%      This is the MATLAB version of program 2.5 from page 39 of 
% "Modeling Infectious Disease in humans and animals" 
% by Keeling & Rohani.
% 
% It is the simple SIS epidemic without births or deaths.
%

% Sets up default parameters if necessary.
if nargin == 0
   beta=0.181;
   gamma=0.1428;
   I0=1e-6;
   MaxTime=1000;
end

S0=1-I0;

% Checks all the parameters are valid
if S0<=0 
    error('Initial level of susceptibles (%g) is less than or equal to zero',S0);
end

if I0<=0 
    error('Initial level of infecteds (%g) is less than or equal to zero',I0);
end

if beta<=0 
    error('Transmission rate beta (%g) is less than or equal to zero',beta);
end
    
if gamma<=0 
    error('Recovery rate gamma (%g) is less than or equal to zero',gamma);
end
    
if MaxTime<=0 
    error('Maximum run time (%g) is less than or equal to zero',MaxTime);
end

if beta<gamma
    warning('Basic reproductive ratio (R_0=%g) is less than one',beta/gamma);
end

S=S0; I=I0;

% The main iteration 
%[~,n] = size(beta);
%for i = 1:n
[t, pop]=ode45(@Diff_2_5,[0 MaxTime],[S I],[],[beta gamma]);

S=pop(:,1); I=pop(:,2);

%hold on
%plot(t,I)
%xlabel('Tempo(Dias)');
%ylabel('Proporção da População');
%hold off
%end

%plots the graphs with scaled colours
subplot(2,1,1)
plot(t,S,'-g'); 
xlabel 'Tempo(Dias)';
ylabel 'Susceptível';

subplot(2,1,2) 
plot(t,I,'-r');
xlabel 'Tempo(Dias)';
ylabel 'Infectado';

%Plot(t,I,'r',t,S,'b');
%title 'Simulação Modelo SIS para Hospitalização por Rotavirose';
%xlabel 'Tempo (Dias)';
%ylabel 'Proporção de Infectados';
%legend('Infectados ','Susceptíveis ')




% Calculates the differential rates used in the integration.
function dPop=Diff_2_5(~,pop, parameter)

beta=parameter(1); gamma=parameter(2);
S=pop(1); I=pop(2);

dPop=zeros(2,1);

dPop(1)= -beta*S*I + gamma*I;
dPop(2)= beta*S*I - gamma*I;


