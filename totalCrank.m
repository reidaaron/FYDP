function C_bean = totalCrank(Deff,r,t,Co)
  % CO2_bean
  % Gives the total amount of mg of co2 released per mg of co2 per
  % gram of coffee
  % 
  % Deff is the effective diffusivity [m^2/s]
  % r is the radius of the coffee  [microns]
  % t is the time to evaluate at [seconds]
  % Co is the co2 in a bean in [mg of Co2/g coffe]
 
  % convert radius to m
  r = r*1E-2;

  % determine the Co2 release in [mg of Co2/g coffee]
  C_bean = 0;
  for i = 1:200
    exponent = -Deff * i^2 * 3.14^2 /(r^2);
    C_bean = C_bean + 1/i^2 * exp(exponent*t);
  end
  C_bean = Co - 6*Co/(3.14)^2 .* C_bean;
end