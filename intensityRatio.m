clear all
clc

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
%% Coffee specifications
  totalCoffee  =1000; % total mass of coffee [grams]
  densi= 561E3;   % Density of coffee bean [g/m^3]
  rbean= 0.057;   % [cm]
%% Environment speficications
  headspaceVol = 300; % headspace volume [ml]
  temperature  =  23; % temperature [degC]
  tend = 480;         % simulation time in hours
%% Membrane specifications
  memThick = 5;   % thickness of membrane [microns]
  Asurf = 5;      % membrane surface area [cm^2]

% plot ratio of pyridie/2-methylbutanol and hexanal/2-methylbutanol
hexanal  = find(strcmp(species,'Hexanal'));
pyridine = find(strcmp(species,'Pyridine'));
methylbutanol = find(strcmp(species,'2-methylpropanal'));

% solve systems
[timesol,P,release_percent] = membrane( temperature, headspaceVol, species,Pig,Deff,Cinf,Pstar,Henry,mw,Pl,yhs_0,tend, totalCoffee, densi, rbean, memThick, Asurf );
%[valve_T,valve_P,valve_R] = valve( temperature, headspaceVol, species,Pig,Deff,Cinf,Pstar,Henry,mw,Pl,yhs_0,tend, totalCoffee, densi, rbean );

tspan = 0:tend*3600;
y_hexanal = interp1(timesol{hexanal},release_percent{hexanal},tspan);
y_pyridine = interp1(timesol{pyridine},release_percent{pyridine},tspan);
y_methbutan = interp1(timesol{methylbutanol},release_percent{methylbutanol},tspan);
tspan = tspan./(3600*24);

% membrane intensity ratio
plot(tspan,y_methbutan./y_pyridine,tspan, y_methbutan./y_hexanal, tspan, y_hexanal./y_pyridine);
title('Membrane'); xlabel('Time [Days]'); ylabel('Intensity Ratio');
legend('Methbu/pyridine','methybut/hexanal','hexanal/pyridein');

% % valve intensity ratio
% valve_T = valve_T./3600;
% v_hexanal = valve_R(hexanal,:);
% v_pyridin = valve_R(pyridine,:);
% v_methbut = valve_R(methylbutanol,:);
% h_p = v_methbut./v_pyridin; m_h = v_methbut./v_hexanal; m_p = v_methbut./v_pyridin;
% plot(m_p(1:5));%,valve_T,m_h,valve_T,h_p);