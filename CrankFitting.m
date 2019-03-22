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

% Crank's model parameters
Deff0 = 10E-12 * ones(size(compounds)); % initial guess for Crank's model
Deff = zeros(size(Deff0));
R = cell(size(Deff0));
tR = R;

% Weibull Distribution's model parameters
k0 = 2.2E-7 * ones(size(compounds)); % initial guess for k for Weibull Distribution
k = zeros(size(k0));
n0 = 0.50 * ones(size(compounds));
n = zeros(size(n0));
Weibullp0 = [k0 n0]
Rw = cell(size(compounds));
tRw = Rw;

Co = 1; % initial concentration as a percentage

figure(1);
%figure('Position',[0 0 1500 1125]); set(gca,'FontSize',16);
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
    weibull = @(p, t) (1 - exp(-(p(1)*t).^p(2)));

	% solve for Deff using initial guess Deff0
	[Deff(i), R{i}, J, COVB, MSE, EMI] = nlinfit(ti, M_releasedi, totalCrankm, Deff0(i));
    t_pred = linspace(min(ti), max(ti));
    [Weibullp, Rw{i}, J, COVB, MSE, EMI] = nlinfit(ti, M_releasedi, weibull, Weibullp0(i,:));
    M_pred = totalCrank(Deff(i), ri, t_pred, ones(size(t_pred)));
    M_pred2 = weibull(Weibullp, t_pred);

    k(i) = Weibullp(1);
    n(i) = Weibullp(2);

    % convert units back to day
    ti = ti / 3600 / 24;
    t_pred = t_pred / 3600 / 24;
    plot(ti, M_releasedi,'o','Color',colours(i),'DisplayName',...
        strcat(species_name{i},' - measured'));
    plot(t_pred, M_pred,'-','Color',colours(i),'DisplayName',...
        strcat(species_name{i},' - Crank'));
    plot(t_pred, M_pred2,'--','Color',colours(i),'DisplayName',...
        strcat(species_name{i},' - Weibull'));
    xlabel('Time [days]');
    ylabel('Mass released [%]');
    title('Curve fitting for Diffusion Coefficients');
    legend('Location', 'se');
    legend show
    
    tR{i} = ti;

    %strcat(species_name{i}, ': ', num2str(Deff(i)))

end

for i = 1:length(Deff0)
    figure(i+2);
    plot(tR{i},R{i},'o',tR{i}, zeros(size(tR{i})), '--');
    xlabel('Time [days]');
    ylabel('Residuals');
    title(strcat('Diffusion Coefficient curve fitting for',{' '},species_name{i},...
        ' using Crank''s Model'));
end

figure_num = i;

for i = 1:length(compounds)
    figure(i+figure_num+2);
    plot(tR{i},Rw{i},'o',tR{i}, zeros(size(tR{i})), '--');
    xlabel('Time [days]');
    ylabel('Residuals');
title(strcat('Diffusion Coefficient curve fitting for',{' '},species_name{i},...
    ' using Weibull Distribution'));
end

hold off
