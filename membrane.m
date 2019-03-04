clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Configuration
%%
%% Environment speficications
  headspaceVol = 300; % headspace volume [ml]
  temperature  =  23; % temperature [degC]
  tend = 100;         % simulation time in hours
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
spacing=tend*36001E-3;
tspan = 0:spacing:tend*3600;

% convert pressue initial conditions to mass
mi0 = mw .* yhs_0 .* Pl .* Vhead .*(10^-6)./( R .* (T+273.15) );

% define plotting placeholders
timeatm = cell(1,length(species));
Patm = zeros(size(tspan));
mgHS = cell(1,length(species));

% for all 3 species loop through and solve system
for k = 1:length(species)
  % Define helper functions
  % J - given mass of species, calculates flux [mg/cm^2 per s]
  % P - calculates pressure of species in headspace [atm]
  J = @(miHS) membraneFlux(memThick,Vhead,T,R,Pig(k),GMV,mw(k),miHS,Pl(k));
  P = @(miHS) (miHS./mw(k)) .*(R.*(T+273.15)./Vhead) .* (1/10^-6);

  % get equilibrium pressure
  C_bean = @(t,miHS) Cinf(k) - miHS/totalCoffee;
  C_mol = @(t,miHS) C_bean(t,miHS) * densi ./ mw(k);
  Peqbrm = @(t,miHS) Henry(k) .* C_mol(t,miHS);

  % define ode system
  dRel_fun= @(t,miHS) dCrank(Deff(k),rbean,t,Cinf(k)) .* totalCoffee .* (P(miHS)<Peqbrm(t,miHS));
  Rel_fun = @(t,miHS) totalCrank(Deff(k),rbean,t,Cinf(k)) .* totalCoffee;
  sys = @(t,miHS) dRel_fun(t,miHS) - Asurf .* J(miHS);

  % solve ode system for mass of species in headspace.
  % store as pressure in headspace
  sol = ode45(sys, tspan, mi0(k));
  timeatm{k} = sol.x./3600; 

  % build profile of release from bean of species
  in_headspace = zeros(size(sol.x));
  in_bean = zeros(size(sol.x));
  total_released = zeros(size(sol.x));
  for i = 1:length(sol.x)
    t = sol.x(i);
    miHS = sol.y(i);
    in_headspace(i) = sol.y(i);
    total_released(i) = sum(in_headspace);
    in_bean(i) = Cinf(k) * totalCoffee - total_released(i);
  end
  mgHS{k} = 100-(1 - in_bean/(Cinf(k)*totalCoffee)) .* 100;

  % getting total pressure
  Patm = Patm + P(interp1(sol.x,sol.y,tspan));
end

% plot out the profiles of all species
figure;
subplot(2,1,1);
hold on
plot(tspan./3600, Patm,'LineWidth',2 );
xlabel('Time [hours]');
ylabel('Pressure in Bag [barg]');
title('Pressure profile in headspace for membrane');
hold off

subplot(2,1,2);
%figure('Position',[0 0 1920 1080]); set(gca,'FontSize',16);
hold on
for i=1:length(species)
  plot(timeatm{i},mgHS{i},'DisplayName',species{i},'LineWidth',4);
end
xlabel('Time [hours]','fontsize',30);
ylabel('Percent Retained','fontsize',30);
title('Release profile of coffee bean','fontsize',30);
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