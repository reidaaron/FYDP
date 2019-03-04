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
  % P - calculates pressure of species in headspace [atm]
  P = @(miHS) (miHS./mw(k)) .*(R.*(T+273.15)./Vhead) .* (1/10^-6);

  % define ode system
  dRel_fun= @(t,miHS) dCrank(Deff(k),rbean,t,Cinf(k)) .* totalCoffee;
  Rel_fun = @(t,miHS) totalCrank(Deff(k),rbean,t,Cinf(k)) .* totalCoffee;
  sys = @(t,miHS) dRel_fun(t,miHS);

  % solve ode system for mass of species in headspace.
  % store as pressure in headspace
  sol = ode45(sys, tspan, mi0(k));
  timeatm{k} = sol.x./3600; 
  Patm{k} = P(sol.y);
end

% plot out the profiles of all species
figure;
%figure('Position',[0 0 1500 1125]); set(gca,'FontSize',16);
hold on
for i=1
  plot(timeatm{i},Patm{i},'LineWidth',4);
end
xlabel('Time [hours]','fontsize',30);
ylabel('Pressure of Carbon Dioxide [bar]','fontsize',30);
title('Pressure profile in headspace without release mechanism','fontsize',30);
hold off
