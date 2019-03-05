function [timespan,P_atm,released_percent] = membrane( T, Vhead, species,Pig,Deff,Cinf,Pstar,Henry,mw,Pl,yhs_0, tend, totalCoffee,density,rbean, memThick, Asurf )

  % create useful shorthands
  RcmHg= 8.205E-5;% ideal gas constant [atm m^3/(mol K)]
  GMV  = 22414;   % molar volume of ideal gas [cm^3 (STP)/mol]
  R = RcmHg;

  % convert pressue initial conditions to mass
  mi0 = mw .* yhs_0 .* Pl .* Vhead .*(10^-6)./( R .* (T+273.15) );

  % define plotting placeholders
  timeatm = cell(1,length(species));
  released_percent = cell(1,length(species));
  timespan = cell(1,length(species));
  P_atm = cell(1,length(species));

  % for all 3 species loop through and solve system
  for k = 1:length(species)
    % Define helper functions
    % J - given mass of species, calculates flux [mg/cm^2 per s]
    % P - calculates pressure of species in headspace [atm]
    J = @(miHS) membraneFlux(memThick,Vhead,T,R,Pig(k),GMV,mw(k),miHS,Pl(k));
    P = @(miHS) (miHS./mw(k)) .*(R.*(T+273.15)./Vhead) .* (1/10^-6);

    % get equilibrium pressure
    C_bean = @(t,miHS) Cinf(k) - miHS/totalCoffee;
    C_mol = @(t,miHS) C_bean(t,miHS) * density ./ mw(k);
    Peqbrm = @(t,miHS) Henry(k) .* C_mol(t,miHS);

    % define ode system
    dRel_fun= @(t,miHS) dCrank(Deff(k),rbean,t,Cinf(k)) .* totalCoffee .* (P(miHS)<Peqbrm(t,miHS));
    Rel_fun = @(t,miHS) totalCrank(Deff(k),rbean,t,Cinf(k)) .* totalCoffee;
    sys = @(t,miHS) dRel_fun(t,miHS) - Asurf .* J(miHS);

    % solve ode system for mass of species in headspace.
    % store as pressure in headspace
    sol = ode45(sys, [0, tend*3600], mi0(k));
    timeatm{k} = sol.x./3600; 

    % build profile of release from bean of species
    in_headspace = zeros(size(sol.x));
    in_bean = zeros(size(sol.x));
    total_released = zeros(size(sol.x));
    for i = 1:length(sol.x)
      t = sol.x(i);
      miHS = sol.y(i);
      in_headspace(i) = sol.y(i);
      total_released(i) = sum(in_headspace);
      in_bean(i) = Cinf(k) * totalCoffee - total_released(i);
    end
    released_percent{k} = 100-(1 - in_bean/(Cinf(k)*totalCoffee)) .* 100;
    P_atm{k} = P(sol.y);
    timespan{k} = sol.x;
  end
end

function J = membraneFlux(memThick, Vhead, T, R, Pig, GMV, mwi, miHS, PiL)
  % membraneFlux (returns flux [mg/cm^2 per s] )
  %  Calculates the flux through the membrane
  % Input:
  %   memThick  membrane thickness [um]
  %   Vhead     headspace volume [ml]
  %   T         temperature [degC]
  %   R         ideal gas constant to use
  %   Pig       permiability through membrane
  %   GMV       molar volume of ideal gas
  %   mwi       molecular weight of species [mg/mol]
  %   miHS      mass of species in headspace [mg]
  %   Pil       pressure of species outside of membrane [atm]

  % convert given mass miHS to pressure using ideal gas law [atm]
  PiHS = (miHS/mwi)*(R*(T+273.15)/Vhead) * (1/10^-6);
  
  % calculate flux of species [mg/cm^2 per s]
  J = (Pig/memThick)*(1/10^-4)*(mwi/GMV)*(PiHS-PiL);
end