function y = rupture(P, species, time, tend, P_rupt)
  tspan = 0:tend*3600;
  P_atm = zeros(size(tspan));

  for i=1:length(species)
    P_i = interp1(time{i},P{i},tspan);
    P_atm = P_atm + P_i;
  end

  y = max(P_atm) > P_rupt;
end