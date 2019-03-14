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