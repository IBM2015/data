function [t,S,I,R] = SIR_Seazonallity(beta0,beta1,gamma,mu,S0,I0,MaxTime)
%
% 
%
% Program_5_1( beta0, beta1, gamma, mu, S0, I0, MaxTime)
%      This is the MATLAB version of program 5.1 from page 160 of 
% "Modeling Infectious Disease in humans and animals" 
% by Keeling & Rohani.
% 
% It is the simple SIR epidemic with sinusoidal forcing of the transmission rate.
% Note: setting beta1 too high can cause numerical difficulties.
%
% This code can also be used to generate bifurcation diagrams, by setting
% beta1 equal to a vector of seasonality rates. The bifurcation diagram is
% constructed using extrapolated initial conditions. Try:
% Program_5_1(17/13,[0:0.001:0.25],1/13,1/(50*365),1/17,1e-4,20*365);
%


% Sets up default parameters if necessary.
if nargin == 0
   beta0=17/13;
   beta1=0.1;
   gamma=1/13.0;
   mu=1/(50*365.0);
   S0=1/17;
   I0=1e-4;
   MaxTime=60*365;
end

% Checks all the parameters are valid
if S0<=0 
    error('Initial level of susceptibles (%g) is less than or equal to zero',S0);
end

if I0<=0 
    error('Initial level of infecteds (%g) is less than or equal to zero',I0);
end

if beta0<=0 
    error('Transmission rate beta0 (%g) is less than or equal to zero',beta0);
end

if beta1<0 
    error('Seasonality beta1 (%g) is less than zero',beta1);
end

if beta1>=1
    error('Seasonality beta1 (%g) is greater than or equal to one',beta1);
end
    
if gamma<=0 
    error('Recovery rate gamma (%g) is less than or equal to zero',gamma);
end

if mu<0 
    error('Birth / Death rate gamma (%g) is less than zero',mu);
end
    
if MaxTime<=0 
    error('Maximum run time (%g) is less than or equal to zero',MaxTime);
end
    
if S0+I0>1
    warning('Initial level of susceptibles+infecteds (%g+%g=%g) is greater than one',S0,I0,S0+I0);
end

if beta0<gamma+mu
    warning('Basic reproductive ratio (R_0=%g) is less than one',beta0/(gamma+mu));
end

S=S0; I=I0; R=1-S-I;

% The main iteration 

if length(beta1)==1
    
    options = odeset('RelTol', 1e-5);
    [t, pop]=ode45(@Diff_5_1,[0 MaxTime],[S I R],options,[beta0 beta1 gamma mu]);

    T=t/365; S=pop(:,1); I=pop(:,2); R=pop(:,3);

    % plots the graphs with scaled colours
    subplot(3,1,1)
    plot(T,S,'-g');
    xlabel 'Time (years)';
    ylabel 'Susceptibles'

    subplot(3,1,2)
    plot(T,I,'-r');
    xlabel 'Time (years)';
    ylabel 'Infectious'

    subplot(3,1,3)
    plot(T,R,'-k');
    xlabel 'Time (years)';
    ylabel 'Recovereds'

else
    if MaxTime<3650
        MaxTime=3650;
    end
    Last=[S0 I0 R];
    for loop=1:length(beta1)
        B1=beta1(loop);

        options = odeset('RelTol', 1e-5);
        [t, pop]=ode45(@Diff_5_1,[0 MaxTime],Last,options,[beta0 B1 gamma mu]);
    
        Last=[pop(end,1) pop(end,2) pop(end,3)];
        Bifur_I(loop,1:10)=interp1(t,pop(:,2),MaxTime-[0:9]*365,'linear');
        
        semilogy(beta1(1:loop),Bifur_I(1:loop,:),'.k');
        xlabel 'Seasonality, \beta_1'
        ylabel 'Level of Infection'
        set(gca,'XLim',[min(beta1) max(beta1)]); drawnow;
    end
    
    semilogy(beta1,Bifur_I,'.k');
    xlabel 'Seasonality, \beta_1'
    ylabel 'Level of Infection'
end

% Calculates the differential rates used in the integration.
function dPop=Diff_5_1(t,pop, parameter)

beta0=parameter(1); beta1=parameter(2); gamma=parameter(3); mu=parameter(4);
beta=beta0*(1+beta1*sin(2*pi*t/365));

S=pop(1); I=pop(2); R=pop(3);

dPop=zeros(3,1);

dPop(1)= mu -beta*S*I - mu*S;
dPop(2)= beta*S*I - gamma*I - mu*I;
dPop(3)= gamma*I - mu*R;


