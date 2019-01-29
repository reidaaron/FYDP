%% Coffee release function
function dC_bean = dCrank(Deffa,r,t,Coa)
  % computes the derivative of crank
  % same arguments as above function
  % the negative is to give rate of release
  r=r*1E-2;
  dC_bean = zeros(size(Deffa));

  for j = 1:length(Deffa)
    Deff = Deffa(j); Co=Coa(j);

    for i = 1:200
      exponent = -Deff * i^2 * 3.14^2 /(r^2);
      dC_bean(j) = dC_bean(j) + exponent/i^2 * exp(exponent*t);
    end
    dC_bean(j) = -6*Co/(3.14)^2 * dC_bean(j);
  end
end