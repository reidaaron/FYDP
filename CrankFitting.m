clear all
clc

[~, species] = xlsread('CranksModelInfo.xlsx',1,'A:A'); 
species = species(2:end);
species_name = unique(species);

data = xlsread('CranksModelInfo.xlsx',1,sprintf('B2:G%d',1+length(species)));

totalmass = data(:, 1);
species_id = data(:, 2);
time = 24 * 3600 * data(:, 3);
M_bean = data(:, 4);
M_released = data(:, 5);
r = data(:,6);

compounds = unique(species_id);

Deff0 = 10E-12 * ones(size(compounds)); % initial guess 
Deff = zeros(size(Deff0));

Co = 1; % initial concentration as a percentage

figure;
hold on


for i = 1:length(compounds)
    indx = find(species_id==i);
    ri = max(r(indx));
    ti = time(indx);
    M_releasedi = M_released(indx);
    Co = ones(size(ti));
    
	% anonymous function in order to input r and Co are constants
	totalCrankm = @(Deff, t) totalCrank(Deff, ri, t, Co);

	% solve for Deff using initial guess Deff0
	Deff(i) = nlinfit(ti, M_releasedi, totalCrankm, Deff0(i));
    M_pred = totalCrankm(Deff(i), ti);
    
    ti = ti / 3600 / 24; % convert units back to day
    
    plot(ti, M_releasedi, 'o', 'DisplayName', strcat(species_name{i},' - measured'));
    plot(ti, M_pred, '-', 'DisplayName', strcat(species_name{i}, ' - predicted'));
    xlabel('Time [days]');
    ylabel('Mass released [%]');
    title('Curve fitting for Diffusion Coefficients');
    legend show
    
    strcat(species_name{i}, ': ', num2str(Deff(i)))
    
end

hold off
