clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Configuration
%%
%% Environment speficications
  headspaceVol = 300; % headspace volume [ml]
  temperature  =  23; % temperature [degC]
  tend = 080;         % simulation time in hours
%%
%% Coffee specifications
  totalCoffee  =1000; % total mass of coffee [grams]
  densi= 561E3;   % Density of coffee bean [g/m^3]
  rbean= 0.057;   % [cm]
%%
%% Valve specifications
%% https://www.wojinvalve.com/coffee-valve/wj1902.html
%% Note: We assume flow rate of valve is linear
%%       see implementation in valveFlux formula
  number = 01;    % number of valves
  diam   = 20E-3; % [m]
  closeP = 2E-3;  % [bar]
  openMP = 8E-3;  % [bar]
  Qref   = 1.3;   % measured air flow rate [L/min]
  Qp     = 6E-3;  % pressure flow rate is measured at [bar]
%%
%% Read in information from Excel spreadsheet
  [~,species] = xlsread('Workbook1.xlsx',1,'A:A'); 
  species=species(2:end);
  Deff = xlsread('Workbook1.xlsx',1,sprintf('C2:C%d',1+length(species)));
  Cinf = xlsread('Workbook1.xlsx',1,sprintf('D2:D%d',1+length(species)));
  Pstar = xlsread('Workbook1.xlsx',1,sprintf('E2:E%d',1+length(species)));
  Henry = xlsread('Workbook1.xlsx',1,sprintf('F2:F%d',1+length(species)));
  mw = xlsread('Workbook1.xlsx',1,sprintf('G2:G%d',1+length(species)));
  Pl = xlsread('Workbook1.xlsx',1,sprintf('H2:H%d',1+length(species)));
  yhs_0 = xlsread('Workbook1.xlsx',1,sprintf('I2:I%d',1+length(species)));
%% Thermodynamic Constants and Properties
  RcmHg= 8.205E-5;% ideal gas constant [atm m^3/(mol K)]
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = temperature+273.15;
R = RcmHg; Vhead = headspaceVol;
tspan = 0:tend*3600;

% convert pressue initial conditions to mass
mi0 = mw .* yhs_0 .* Pl .* Vhead .*(10^-6)./( R .* (T+273.15) );

% define plotting placeholders
timeatm = cell(1,length(species));
Patm = cell(1,length(species));
mgHS = cell(1,length(species));

% Helper functions
  % Pressure in atm
  P      = @(miHS) (miHS./mw) .*(R.*(T+273.15)./Vhead) .* (1/10^-6);
  % valve position given the headspaces composition
  ValveP = @(miHS) valvePosition( closeP, openMP, sum(P(miHS)) );
  valveref = valvePosition(closeP, openMP, Qp);
  % flow rate through the valve in mol/s
  Q      = @(miHS) flowrate( T, Qref, valveref, ValveP(miHS) );

% define flux through valve (area is not included since it cancels out)
J = @(miHS) number * Q(miHS) .* ones(size(species)) ;

% get equilibrium pressure
C_bean = @(t,miHS) Cinf - miHS/totalCoffee;
C_mol = @(t,miHS) C_bean(t,miHS) * densi ./ mw;
Peqbrm = @(t,miHS) Henry .* C_mol(t,miHS);

% define ode system
dRel_fun= @(t,miHS) dCrank(Deff,rbean,t,Cinf) .* totalCoffee .* (P(miHS)<Peqbrm(t,miHS));
sys = @(t,miHS) dRel_fun(t,miHS) - P(miHS)/(sum(P(miHS))) .* J(miHS);

% solve ode system and store in species values
sol = ode45(sys, tspan, mi0);

% generate valve plot
y= [];
for i=1:length(sol.y)
  y(i) = ValveP(sol.y(:,i));
end

% get release profile from bean
retained = zeros(length(species),length(sol.x));
for i=1:length(species)
  miHS = sol.y(i,:);
  total_released = zeros(size(miHS));
  in_headspace = zeros(size(miHS));

  for j=2:length(miHS)
    total_released(j) = trapz(sol.x(1:j),sol.y(1:j));
  end

  in_bean = Cinf(i) * totalCoffee - total_released;
  retained(i,:) = total_released;%100-(1 - in_bean/(Cinf(i)*totalCoffee)) .* 100;
end

% plot out profiles of all species
figure;
  subplot(3,1,1);
  hold on
  Phead = P(sol.y);
  for i = 1:length(species)
    if i == 1
      yyaxis left
    else
      yyaxis right
    end
    plot(sol.x./3600,Phead(i,:),'DisplayName',species{i});
  end
  xlabel('Time [hours]');
  yyaxis left
  ylabel('Presusre of species [bar]');
  title('Pressure profile in headspace for one-way valve');
  legend show;
  hold off

  subplot(3,1,2);
  hold on
  for i = 1:length(species)
    plot(sol.x./3600,retained(i,:),'DisplayName',species{i});
  end
  xlabel('Time [hours]');
  ylabel('Percent Retained');
  title('Retention profile of coffee bean');
  legend show
  hold off

  subplot(3,1,3);
  hold on
  plot(sol.x./3600,y);
  xlabel('Time [hours]');
  ylabel('Valve position');
  title('Valve position profile');
  hold off

  y0 = y;


%% Function to model valve behaviour
function Q = flowrate(T, Qref, valveRef, valvePosition)
  % convert from L/min to mol/s
  Q = Qref./valveRef .* valvePosition;
  Q = Q./60 .* (T/273.15);
end
function y = valvePosition(closeP,openMP,Psys)
  if Psys <= closeP
    y = 0;
  elseif Psys >= openMP
    y = 1;
  else
    y = (Psys-closeP)./(openMP-closeP);
  end
end