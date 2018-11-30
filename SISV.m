%Programa: Modelo SIS com controle por Vacinação Pediátrica
%Propósito: Explorar um modelo Simples, com controle por vacinação simples,
%para a epidemia de rotavirose no Estado de São Paulo.
%Autor: Aldo Parada Hurtado
%Data: 06/Set/2018
%Observação: Este modelo é baseado nos modelos SIS e SIR das páginas 39 e
%293 do livro "Modeling Infectious Disease in humans and animals" 
% de Keeling & Rohani (2008).

% Função do Programa
function [T,S,I,V] = SISV(beta,gamma,S0,I0,p, tV, MaxTime)

% Sets up default parameters if necessary.
if nargin == 0
   beta=0.181;
   gamma=0.1428;
   S0=1-1e-6;
   I0=1e-6;
   p=0.9148;
   tV=500;
   MaxTime=1000;
end

% Verificação da Validade dos Imputados aos Parâmetros
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

if beta<gamma
    warning('Basic reproductive ratio (R_0=%g) is less than one',beta/(gamma+mu));
end

% Condições Iniciais
S=S0; I=I0; V=1-S-I;

% Integração das Equações Diferenciais Ordinárias pelo ODE45 do Matlab 
options = odeset('RelTol', 1e-5);
% Integração antes do Controle por Vacina (p=0)
[t, pop]=ode45(@Diff_SISV,[0 tV],[S I V],options,[beta gamma 0]);
T=t; S=pop(:,1); I=pop(:,2); V=pop(:,3);
% Integração após o Controle por Vacina (p>0)
[t, pop]=ode45(@Diff_SISV,[tV MaxTime],pop(end,:),options,[beta gamma p]);
T=[T; t]; S=[S; pop(:,1)]; I=[I; pop(:,2)]; V=[V; pop(:,3)];

% Gera o Gráfico da Simulação
plot(T,S,'-g',T,I,'-r',T,V,'-k'); 
hold on; Y=get(gca,'YLim'); plot(tV+[0 0],Y,'--k'); hold off
title 'Modelo SIR com Vacinação para Rotavirose';
xlabel 'Tempo(Dias)';
ylabel 'Proporção da População';


% Calcula as taxas diferenciais usadas na integração
function dPop=Diff_SISV(~,pop, parameter)

beta=parameter(1); gamma=parameter(2);  p=parameter(3);
S=pop(1); I=pop(2); V=pop(3);

dPop=zeros(3,1);

dPop(1)= - (1 - p)*beta*S*I + gamma*I;
dPop(2)= beta*S*I - gamma*I;
dPop(3)= pS;
