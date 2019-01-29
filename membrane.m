clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Configuration
%%
%% Environment speficications
  headspaceVol = 300; % headspace volume [ml]
  temperature  =  23; % temperature [degC]
  tend = 001;         % simulation time in hours
%%
%% Membrane specifications
  memThick = 5;   % thickness of membrane [microns]
  Asurf = 5;      % membrane surface area [cm^2]
%%
%% Coffee specifications
  totalCoffee  =1000; % total mass of coffee [grams]
  densi= 561E3;   % Density of coffee bean [g/m^3]
  rbean= 0.057;   % [cm]
%%
%% Read in information from Excel spreadsheet
  % Read in all the species
  [~,species] = xlsread('Workbook1.xlsx',1,'A:A'); species=species(2:end);
  Pig = xlsread('Workbook1.xlsx',1,sprintf('B2:B%d',1+length(species)));
  Deff = xlsread('Workbook1.xlsx',1,sprintf('C2:C%d',1+length(species)));
  Cinf = xlsread('Workbook1.xlsx',1,sprintf('D2:D%d',1+length(species)));
  Pstar = xlsread('Workbook1.xlsx',1,sprintf('E2:E%d',1+length(species)));
  Henry = xlsread('Workbook1.xlsx',1,sprintf('F2:F%d',1+length(species)));
  mw = xlsread('Workbook1.xlsx',1,sprintf('G2:G%d',1+length(species)));
  Pl = xlsread('Workbook1.xlsx',1,sprintf('H2:H%d',1+length(species)));
  yhs_0 = xlsread('Workbook1.xlsx',1,sprintf('I2:I%d',1+length(species)));
%% Thermodynamic Constants and Properties
  RcmHg= 8.205E-5;% ideal gas constant [atm m^3/(mol K)]
  GMV  = 22414;   % molar volume of ideal gas [cm^3 (STP)/mol]
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create useful shorthands
T = temperature; R = RcmHg; Vhead = headspaceVol;
tspan = 0:tend*3600;

% convert pressue initial conditions to mass
mi0 = mw .* yhs_0 .* Pl .* Vhead .*(10^-6)./( R .* (T+273.15) );

% get equilibrium partial pressure using 
% henries law
C_mol  = Cinf * densi ./ mw;
Peqbrm = Henry.*C_mol; 

% define plotting placeholders
timeatm = cell(1,length(species));
Patm = cell(1,length(species));
mgHS = cell(1,length(species));

% for all 3 species loop through and solve system
for k = 1:length(species)
  % Define helper functions
  % J - given mass of species, calculates flux [mg/cm^2 per s]
  % P - calculates pressure of species in headspace [atm]
  J = @(miHS) membraneFlux(memThick,Vhead,T,R,Pig(k),GMV,mw(k),miHS,Pl(k));
  P = @(miHS) (miHS./mw(k)) .*(R.*(T+273.15)./Vhead) .* (1/10^-6);

  % define ode system
  dRel_fun= @(t,miHS) dCrank(Deff(k),rbean,t,Cinf(k)) .* totalCoffee .* (P(miHS)<Peqbrm(k));
  Rel_fun = @(t,miHS) totalCrank(Deff(k),rbean,t,Cinf(k)) .* totalCoffee;
  sys = @(t,miHS) dRel_fun(t,miHS) - Asurf .* J(miHS);

  % solve ode system for mass of species in headspace.
  % store as pressure in headspace
  sol = ode45(sys, tspan, mi0(k));
  timeatm{k} = sol.x./3600; 
  Patm{k} = P(sol.y);

  % build profile of release from bean of species
  % or consumption in the case of oxygen (to consider oxidation)
  release_profile = zeros(size(sol.x));
  for i = 1:length(sol.x)
    release_profile(i) = Rel_fun(sol.x(i)) .* (Peqbrm(k)-P(sol.y(i))) .* ((Peqbrm(k)-P(sol.y(i)))>0);
  end
  mgHS{k} = release_profile;
end

% plot out the profiles of all species
subplot(2,1,1);
hold on
for i=1:length(species)
  plot(timeatm{i},Patm{i},'DisplayName',species{i});
end
xlabel('Time [hours]');
ylabel('Pressure of Species [bar]');
title('Pressure profile in headspace for membrane');
legend show
hold off

subplot(2,1,2);
hold on
for i=1:length(species)
  plot(timeatm{i},mgHS{i},'DisplayName',species{i});
end
xlabel('Time [hours]');
ylabel('Total Mass Released [mg]');
title('Release profile of coffee bean');
legend show
hold off

function J = membraneFlux(memThick, Vhead, T, R, Pig, GMV, mwi, miHS, PiL)
  % membraneFlux (returns flux [mg/cm^2 per s] )
  %  Calculates the flux through the membrane
  % Input:
  %   memThick  membrane thickness [um]
  %   Vhead     headspace volume [ml]
  %   T         temperature [degC]
  %   R         ideal gas constant to use
  %   Pig       permiability through membrane
  %   GMV       molar volume of ideal gas
  %   mwi       molecular weight of species [mg/mol]
  %   miHS      mass of species in headspace [mg]
  %   Pil       pressure of species outside of membrane [atm]

  % convert given mass miHS to pressure using ideal gas law [atm]
  PiHS = (miHS/mwi)*(R*(T+273.15)/Vhead) * (1/10^-6);
  
  % calculate flux of species [mg/cm^2 per s]
  J = (Pig/memThick)*(1/10^-4)*(mwi/GMV)*(PiHS-PiL);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coffee Co2 release functions
%
% These functions describe the dynamics of co2 release in coffee using 
% the crank method. Primary assumption is that the coffee is spherical 
% in shape and that the release does not depend on the concentration
% of CO2 in the headspace. 
% The derivative function is analytical.  

function C_bean = crank(Deff,r,t,Co)
  % CO2_bean
  % Gives the amount of mg of co2 in the bean per
  % gram of coffee
  % 
  % Deff is the effective diffusivity [m^2/s]
  % r is the radius of the coffee  [microns]
  % t is the time to evaluate at [seconds]
  % Co is the co2 in a bean in [mg of Co2/g coffe]
 
  % convert radius to m
  r = r*1E-2;

  % determine the Co2 release in [mg of Co2/g coffee]
  C_bean = 0;
  for i = 1:200
    exponent = -Deff * i^2 * 3.14^2 /(r^2);
    C_bean = C_bean + 1/i^2 * exp(exponent*t);
  end
  C_bean = 6*Co/(3.14)^2 * C_bean;
end

function C_bean = totalCrank(Deff,r,t,Co)
  % CO2_bean
  % Gives the total amount of mg of co2 released per mg of co2 per
  % gram of coffee
  % 
  % Deff is the effective diffusivity [m^2/s]
  % r is the radius of the coffee  [microns]
  % t is the time to evaluate at [seconds]
  % Co is the co2 in a bean in [mg of Co2/g coffe]
 
  % convert radius to m
  r = r*1E-2;

  % determine the Co2 release in [mg of Co2/g coffee]
  C_bean = 0;
  for i = 1:200
    exponent = -Deff * i^2 * 3.14^2 /(r^2);
    C_bean = C_bean + 1/i^2 * exp(exponent*t);
  end
  C_bean = Co - 6*Co/(3.14)^2 * C_bean;
end

function dC_bean = dCrank(Deff,r,t,Co)
  % computes the derivative of crank
  % same arguments as above function
  % the negative is to give rate of release
  r=r*1E-2;
  dC_bean = 0;
  for i = 1:200
    exponent = -Deff * i^2 * 3.14^2 /(r^2);
    dC_bean = dC_bean + exponent/i^2 * exp(exponent*t);
  end
  dC_bean = -6*Co/(3.14)^2 * dC_bean;
end