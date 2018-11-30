function [t,SC,IC,SA,IA] = Program_3_3_SIS(beta,gamma,~,mu,S0,I0, MaxTime)
%
% 
%
% Program_3_3( beta, gamma, LC, mu, S0, I0, MaxTime)
%      This is the MATLAB version of program 3.3 from page 79 of 
% "Modeling Infectious Disease in humans and animals" 
% by Keeling & Rohani.
% 
% It is the SIR model with two different age-groups.
% Note we do not explicitly model the recovered classes.
%
% mu, S0 and I0 are all vectors.


% Sets up default parameters if necessary.
if nargin == 0
   beta=[100 10; 10 20];
   gamma=10;
   lC=0.066667;
   mu=[0 0.0166667];
   S0=[0.1 0.1];
   I0=[0.0001 0.0001];
   MaxTime=100;
end

nC=mu(2)/(lC+mu(2));
n=[nC 1-nC];

nu=(lC+mu(2))*nC;


% Checks all the parameters are valid

for a=1:2

    if I0(a)<0
        error('Initial level of infection in age-group %d (%g) is less than zero',a,I0(a));
    end
    if S0(a)<0
        error('Initial level of susceptibles in age-group %d (%g) is less than zero',a,S0(a));
    end
    if mu(a)<0
        error('Death rate in age-group %d (%g) is less than zero',a,mu(a));
    end
    if I0(a)>n(a)
        error('Initial level of infection in age-group %d (%g) is greater than calculated group size (%g)',a,I0(a),n(a));
    end
    if S0(a)>n(a)
        error('Initial level of susceptibles in age-group %d (%g) is greater than calculated group size (%g)',a,S0(a),n(a));
    end
    if S0(a)+I0(a)>n(a)
        error('Initial level of infecteds + susceptibles in age-group %d (%g) is greater than calculated group size (%g)',a,I0(a)+S0(a),n(a));
    end
end

if ~isempty(find(beta<0, 1)) 
    error('Transmission rate beta ([%g %g; %g %g]) is less than zero',beta(1,1),beta(1,2),beta(2,1),beta(2,2));
end
    
if gamma<=0 
    error('Recovery rate gamma (%g) is less than or equal to zero',gamma);
end
    
if MaxTime<=0 
    error('Maximum run time (%g) is less than or equal to zero',MaxTime);
end
    
if beta(1,2)~=beta(2,1)
    warning('Transmission rate beta ([%g %g; %g %g]) is not symmetric',beta(1,1),beta(1,2),beta(2,1),beta(2,2));

end

% The main iteration 
options = odeset('RelTol', 1e-5);
[t, pop]=ode45(@Diff_3_3,[0 MaxTime],[S0(1) I0(1) S0(2) I0(2)],options,[reshape(beta,1,4) gamma lC reshape(mu,1,2)]);

SC=pop(:,1); IC=pop(:,2);
SA=pop(:,3); IA=pop(:,4);

% plots the graphs with scaled colours
subplot(2,1,1)
h=plot(t,IC,'-r',t,IA,'--r'); 
legend(h,'Children','Adults');
xlabel 'Time';
ylabel 'Infectious'

subplot(2,1,2) 
h=plot(t,SC,'-g',t,SA,'--g'); 
legend(h,'Children','Adults',2);
xlabel 'Time';
ylabel 'Susceptibles'


% Calculates the differential rates used in the integration.
function dPop=Diff_3_3(~,pop, parameter)

beta=reshape(parameter(1:4),2,2);  gamma=parameter(5); lC=parameter(6); mu=parameter(7:8);
nC=mu(2)/(lC+mu(2));
nu=(lC+mu(2))*nC;

SC=pop(1); IC=pop(2); SA=pop(3); IA=pop(4);

dPop=zeros(4,1);

dPop(1)=nu - (beta(1,1)*IC+beta(1,2)*IA)*SC - mu(1)*SC - lC*SC;
dPop(2)=(beta(1,1)*IC+beta(1,2)*IA)*SC - gamma*IC - mu(1)*IC - lC*IC;
dPop(3)=lC*SC-(beta(2,1)*IC+beta(2,2)*IA)*SA - mu(2)*SA;
dPop(4)=lC*IC+(beta(2,1)*IC+beta(2,2)*IA)*SA - gamma*IA - mu(2)*IA;



