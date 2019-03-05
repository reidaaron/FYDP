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
  Asurf = 4;      % cm^2

% vary surface area and thickness
Asurf_Range = 1:18;
mem_range = 1:10;
dims = [length(mem_range),length(Asurf_Range)];
cost_area = zeros(dims);
ruptured = zeros(dims);
degradation = zeros(dims);
for i=1:length(Asurf_Range)
  Asurf = Asurf_Range(i);
  for j=1:length(mem_range)
    memThick = mem_range(j);

    [timesol,P,release_percent] = membrane( temperature, headspaceVol, species,Pig,Deff,Cinf,Pstar,Henry,mw,Pl,yhs_0,tend, totalCoffee, densi, rbean, memThick, Asurf );

    % determine if ruptured
    ruptured(j,i) = 1-rupture(P,species,timesol,tend,0.8);

    % get degradation of intensity ratio
    degradation(j,i) = degraded(species,timesol,release_percent);

    % determine cost ($125/ft^2);
    cost_area(j,i) = 0.1345 * Asurf;
  end
end

figure;
Z = degradation.*ruptured;
Z(Z==0) = nan;
surf(Asurf_Range,mem_range,Z);
xlabel('Surface Area [cm^2]');
ylabel('Thickness [um]');
zlabel('% Intensity Degredation');
title('Intensity Degradation after 20 days');

% determine if degraded
function max_degrade = degraded(species, time, release)
  hexanal = find(strcmp(species,'Hexanal'));
  pyridin = find(strcmp(species,'Pyridine'));
  methylb = find(strcmp(species,'2-methylpropanal'));

  tspan = [0 480*3600];
  y_hexanal = interp1(time{hexanal}, release{hexanal}, tspan);
  y_pyridin = interp1(time{pyridin}, release{pyridin}, tspan);
  y_methylb = interp1(time{methylb}, release{methylb}, tspan);

  ratio = [y_methylb./y_pyridin; y_methylb./y_hexanal];
  degrade = zeros(1,length(ratio));
  for i = 1:length(degrade)
    ratio_i = ratio(i,:);
    degrade(i) = ( ratio_i(2) - ratio_i(1) )/ratio_i(1);
  end

  max_degrade = max(abs(degrade)) * 100;
end

% determine if exceeds rupture pressure
function y = rupture(P, species, time, tend, P_rupt)
  tspan = 0:tend*3600;
  P_atm = zeros(size(tspan));

  for i=1:length(species)
    P_i = interp1(time{i},P{i},tspan);
    P_atm = P_atm + P_i;
  end

  y = max(P_atm) > P_rupt;
end