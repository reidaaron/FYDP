function [timespan, P_atm, retained, y] = valve( T, Vhead, species,Pig,Deff,Cinf,Pstar,Henry,mw,Pl,yhs_0, tend, totalCoffee,density,rbean)

  %% Valve specifications
  %% https://www.wojinvalve.com/coffee-valve/wj1902.html
  %% Note: We assume flow rate of valve is linear
  %%       see implementation in valveFlux formula
    number = 04;    % number of valves
    diam   = 20E-3; % [m]
    closeP = 2E-3;  % [bar]
    openMP = 8E-3;  % [bar]
    Qref   = 1.3;   % measured air flow rate [L/min]
    Qp     = 6E-3;  % pressure flow rate is measured at [bar]
  %%
    
  RcmHg= 8.205E-5;% ideal gas constant [atm m^3/(mol K)]
  T = T+273.15; R = RcmHg;

  % convert pressue initial conditions to mass
  mi0 = mw .* yhs_0 .* Pl .* Vhead .*(10^-6)./( R .* (T+273.15) );

  % define plotting placeholders
  timeatm = cell(1,length(species));
  Patm = cell(1,length(species));
  mgHS = cell(1,length(species));

  % Helper functions
    % Pressure in atm
    P      = @(miHS) (miHS./mw) .*(R.*(T+273.15)./Vhead) .* (1/10^-6);
    % valve position given the headspaces composition
    ValveP = @(miHS) valvePosition( closeP, openMP, sum(P(miHS)) );
    valveref = valvePosition(closeP, openMP, Qp);
    % flow rate through the valve in mol/s
    Q      = @(miHS) flowrate( T, Qref, valveref, ValveP(miHS) );

  % define flux through valve (area is not included since it cancels out)
  J = @(miHS) number * Q(miHS) .* ones(size(species)) ;

  % get equilibrium pressure
  C_bean = @(t,miHS) Cinf - miHS/totalCoffee;
  C_mol = @(t,miHS) C_bean(t,miHS) * density ./ mw;
  Peqbrm = @(t,miHS) Henry .* C_mol(t,miHS);

  % define ode system
  dRel_fun= @(t,miHS) dCrank(Deff,rbean,t,Cinf) .* totalCoffee .* (P(miHS)<Peqbrm(t,miHS));
  sys = @(t,miHS) dRel_fun(t,miHS) - P(miHS)/(sum(P(miHS))) .* J(miHS);

  % solve ode system and store in species values
  sol = ode45(sys, [0 tend*3600], mi0);

  % generate valve plot
  y= [];
  for i=1:length(sol.y)
    y(i) = ValveP(sol.y(:,i));
  end

  % get release profile from bean
  retained = zeros(length(species),length(sol.x));
  for i=1:length(species)
    miHS = sol.y(i,:);
    rel = @(t) totalCrank(Deff(i),rbean,t,Cinf(i)) .* totalCoffee;

    for j=1:length(sol.x)
      t = sol.x(j);
      miHS = sol.y(:,j);
      P_eq = P(miHS) < Peqbrm(t,miHS);

      if P_eq(i)
        release_profile(i,j) = rel(sol.x(j));
      else
        if j> 2
          release_profile(i,j) = release_profile(i,j-1);
        else
          release_profile(i,j) = 0;
        end
      end  
    end

    retained(i,:) = (Cinf(i)*totalCoffee -release_profile(i,:) )./(Cinf(i)*totalCoffee) .* 100;
  end

  P_atm = P(sol.y);
  timespan = sol.x;
end


%% Function to model valve behaviour
function Q = flowrate(T, Qref, valveRef, valvePosition)
  % convert from L/min to mol/s
  Q = Qref./valveRef .* valvePosition;
  Q = Q./60 .* (T/273.15);
end
function y = valvePosition(closeP,openMP,Psys)
  if Psys <= closeP
    y = 0;
  elseif Psys >= openMP
    y = 1;
  else
    y = (Psys-closeP)./(openMP-closeP);
  end
end