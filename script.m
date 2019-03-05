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

  % plot out the profiles of all species
[timesol,P,release_percent] = membrane( temperature, headspaceVol, species,Pig,Deff,Cinf,Pstar,Henry,mw,Pl,yhs_0,tend, totalCoffee, densi, rbean, memThick, Asurf );

hold on
for i=1:length(species)
  plot(timesol{i},release_percent{i},'DisplayName',species{i});
end
legend show
hold off