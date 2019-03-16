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

% vary surface area and thickness
Asurf_Range = 1:5;
mem_range = [5:9 10:100:1000 1001:1000:5000 5000];
dims = [length(mem_range),length(Asurf_Range)];
cost_area = zeros(dims);
ruptured = zeros(dims);
degradation = zeros(dims);
Z = zeros(size(dims));

for i=1:length(Asurf_Range)
  Asurf = Asurf_Range(i);
  for j=1:length(mem_range)
    memThick = mem_range(j);

    [timesol,P,release_percent] = membrane( temperature, headspaceVol, species,Pig,Deff,Cinf,Pstar,Henry,mw,Pl,yhs_0,tend, totalCoffee, densi, rbean, memThick, Asurf );

    % determine if ruptured
    ruptured(j,i) = rupture(P,species,timesol,tend,1);

    % get degradation of intensity ratio
    degradation(j,i) = degraded(species,timesol,release_percent);

    % determine cost ($125/ft^2);
    cost_area(j,i) = 0.1345 * Asurf;

    % check to make sure we are wihin reasonable limits
    if degradation(j,i) > 10
      Z(j,i) = NaN;
    elseif ruptured(j,i) == 1
      Z(j,i) = NaN;
    else
      Z(j,i) = cost_area(j,i);
    end
  end
end

% find surface area ( cost =f(Area) )
degradation_const = Z./Z .* degradation;
min_cost = min(min(Z));
col = find(min(Z) == min_cost);
membrane_area = Asurf_Range(col)

% find membrane thickness
min_degradation = min(degradation_const(:,col))
row = find(degradation_const(:,col) == min_degradation);
membrane_thickness = mem_range(row)

% plot Area at fixed thickness
figure;
plot(Asurf_Range,degradation(row,:));
xlabel('Surface Area [cm^{2}]');
ylabel('% Degradation after 20 days');
title(sprintf('Percent Degradation after 20 days at %d um thickness',membrane_thickness));

% plot thickness at fixed Area
figure;
plot(mem_range,degradation(:,col));
xlabel('Membrane Thickness [um]');
ylabel('% Degradation after 20 days');
title(sprintf('Percent Degradation after 20 days at %d cm^{2} area',membrane_area));

% simulate membrane
[timesol,P,release_percent] = membrane( temperature, headspaceVol, species,Pig,Deff,Cinf,Pstar,Henry,mw,Pl,yhs_0,10, totalCoffee, densi, rbean, membrane_thickness, membrane_area );

% plot pressure profile of membrane
figure;
hold on
for i=1:length(species)
  plot(timesol{i}./3600,P{i},'DisplayName', species{i});
  break;
end
hold off
xlabel('Time [hrs]');
ylabel('Pressure [bar]');
title(sprintf('Pressure profile in headspace for membrane with area %d cm^{2} and %d um thickness',membrane_area,membrane_thickness));
legend show

% plot percent release from bean
[timesol,P,release_percent] = membrane( temperature, headspaceVol, species,Pig,Deff,Cinf,Pstar,Henry,mw,Pl,yhs_0,20, totalCoffee, densi, rbean, membrane_thickness, membrane_area );
figure;
hold on
for i=2:length(species)
  plot(timesol{i}./3600,release_percent{i},'DisplayName', species{i});
end
hold off
xlabel('Time [hrs]');
ylabel('Pressure [bar]');
title(sprintf('Percent released for membrane with area %d cm^{2} and %d um thickness',membrane_area,membrane_thickness));
legend show