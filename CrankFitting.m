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
%figure('Position',[0 0 1920 1360]); set(gca,'FontSize',16);
hold on

colours = ['b', 'r', 'g', 'y', 'c', 'm'];

for i = 1:length(compounds)
    indx = find(species_id==i);
    ri = max(r(indx));
    ti = time(indx);
    M_releasedi = M_released(indx);
    Co = ones(size(ti));
    
	% anonymous function in order to input r and Co are constants
	totalCrankm = @(Deff, t) totalCrank(Deff, ri, t, Co);

	% solve for Deff using initial guess Deff0
	[Deff(i), R, J, COVB, MSE, ERRORMODELINFO] = nlinfit(ti, M_releasedi, totalCrankm, Deff0(i));
    t_pred = linspace(min(ti), max(ti));
    M_pred = totalCrank(Deff(i), ri, t_pred, ones(size(t_pred)));
    
    % convert units back to day
    ti = ti / 3600 / 24;
    t_pred = t_pred / 3600 / 24;
    plot(ti, M_releasedi,'o','Color',colours(i),'DisplayName',...
        strcat(species_name{i},' - measured'),'LineWidth',2);
    plot(t_pred, M_pred,'-','Color',colours(i),'DisplayName',...
        strcat(species_name{i},' - predicted'),'LineWidth',4);
    xlabel('Time [days]','fontsize',30);
    ylabel('Mass released [%]','fontsize',30);
    title('Curve fitting for Diffusion Coefficients','fontsize',30);
    legend('Location', 'se');
    legend show
    
    %strcat(species_name{i}, ': ', num2str(Deff(i)))
end

hold off
