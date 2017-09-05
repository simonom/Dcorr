% This file is part of Dcorr.

% Dcorr is free software: you can redistribute it and/or modify it under the 
% terms of the GNU General Public License as published by the Free Software 
% Foundation, either version 3 of the License, or (at your option) any later 
% version.

% Dcorr is distributed in the hope that it will be useful, but WITHOUT ANY 
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
% FOR A PARTICULAR PURPOSE.  See the GNU Public License for more details.

% You should have received a copy of the GNU Public License along with Dcorr.
% If not, see <http://www.gnu.org/licenses/>.


% the equivalent of the analytical solution that is nested in 
% Zaveri_Eq26Eq31_pde_D_dependent, so that a fair comparison can be made
% to the analytical solution there, with its effective D derived from the
% pde solution, rather than use analyt_eq31_D_dependent, which is a 
% slightly different code for the analytical solution because it considers
% volume change during the Newton step method for iteration

% different to other versions as it has a time loop that continues until
% the e-folding state is reached.  This version does contain an option to
% use a composition-dependent diffusion coefficient

function [eq_time, Rpm0store, times, Cratstore, Cg0store, Cratstore3, Dscalestore] = ...
    analyt_eq31_simple5(D, no_ts, Cstar, Rpm, pdetime, pdeRpm0, ...
    Dcorr_it, Nm, Aieq, pdetimer)
    
    % ---------------------------------------------------------------------
    % inputs:
    % D - self-diffusion coefficients of components (semi-volatile first)
    % eq_time - expected diffusion equilibrium time (usually from the pde
    % ) (s)
    % no_ts - number of time steps to solve diffusion over
    % Cstar - effective saturation concentration of semi-volatile (ug/m^3 
    % (air))
    % ml0 - initial mass loading of condensed phase (ug/m^3 (air))
    % Rpm - initial particle radius (m)
    % Cg - initial gas-phase concentration of semi-volatile (ug/m^3 (air))
    % pdetime - times used by the pde approach (s)
    % pdeRpm0 - radius estimated of the pde (m), only used if iteratively
    % finding the D correction factor
    % Dcorr_it - flag for iteratively finding the D correction value (at a
    % certain time step) needed to match pde estimates of r (1=on)
    % Nm - number concentration of particles (# /m^3 (air))
    % Aieq - equilibrium particle-phase mole fraction of semi-volatile on
    % top row, time this mole fraction attained on the bottom row (s)
    % pdetf - flag for working through prescribed pde times (1=on, 0=off)
    % effD_mark - marker for using corrected D (i.e. scaled to pde) (1=on,
    % 0=off)
    % var_es_flag - flag for working through changes to saturation ratio
    % with time (1 for on, 0=off)
    % ---------------------------------------------------------------------
    
    % outputs: 
    % eq_time - e-folding time (s) as predicted by the analytical approach
    
    % ---------------------------------------------------------------------
    % starting values and constants
    
    ti = 1; % count on number of time steps
    
    times = zeros(1,1); % result times (s)
    
    % equilibrium concentration (mole 
    % fraction) of semi-volatile (ug m^{-3} 
    % (air))
    Aieq_now = (abs(pdetime(1)-Aieq(2,1))/(Aieq(2,2)-Aieq(2,1)))* ...
            Aieq(1,2)+ ...
            (abs(pdetime(1)-Aieq(2,2))/(Aieq(2,2)-Aieq(2,1)))*Aieq(1,1);

    if Aieq(1,2) > Aieq(1,1)
        esinc = 1;
    else
        esinc = 0;
    end
    % starting values and constants
    [kci, Db, Rpm, p, Cstar, Ai, Dg, Ceq, M] = ...
        starting_stuff(times, D, Cstar, Rpm, Aieq(1, 1), Aieq_now);
    
    % x mesh (m) and shell volumes (m^3)
    [Vpt] = spat_arrays(Rpm, Nm);
    
    % records of important values used by the analytical solution through
    % time steps
    Vp0store = zeros(1, 1);
    Rpm0store = zeros(1, 1);
    Dscalestore = zeros(1, 1);
    Cg0store = zeros(1, 1);
    Cratstore = zeros(1, 1);
    Cratstore3 = zeros(2, 1);
    Ca0store = zeros(2, 1); % particle-phase concentration (ug/m^3 (air))
    % remember initial values
    Vp0store(1) = Vpt(1);
    Rpm0store(1) = Rpm;
    
    k2store = zeros(1, 1);
    
    % particle-phase concentration of each component at t = 0 s (ug/m^3 
    % (air))
    Ca0store(:,1) = Ai(:, 1).*(Vpt(1));
    % sum of the concentration of all components in the particle-phase 
    % initially (ug/m^3 (particle))
    Aj = sum(Ai(:, 1));

    % mass of non-volatile and semi-volatile component in a single particle 
    % (ug/(single particle))
    mugnv = Ca0store(2, 1)/Nm;
    
    % time loop conditions
    % follow pdetime
    while_var = ti;
    while_con = no_ts-1;
    
    % equilibrium mole fraction of semi-volatile during run
    Aieq_now = Aieq(1,2);
    
    while while_var <= while_con
        
        % time step (s)
        % when moving through pdetimes
        h = pdetime(ti+1)-pdetime(ti);
        
        % recall the total particle-phase volume and the radius from the 
        % end of the previous time step (m^3)
        Vp0 = Vp0store(ti); 
        Rpm0 = Rpm0store(ti);
         
        % initial bulk particle-phase concentration and first guess of new
        % concentration (ug/m^3 (air))
        Ca0 = Ca0store(:, ti); 
        Ca2 = Ca0store(:, ti); 
        
        % total number of moles in bulk particle-phase (mol)
        xT = (Ca2(1)/M(1)+Ca2(2)/M(2));
        % semi-volatile mole fraction in bulk condensed-phase (fraction)
        xsv = (Ca2(1)/M(1))/xT;
        
        xnv = 1-xsv; % nv mole fraction in condensed-phase (fraction)
        % logarithmic dependence of diffusion coefficient on mole fraction 
        % in particle bulk
        Dx = (Db(1)^xsv*Db(2)^xnv);
        % sigmoidal dependence of diffusion coefficient on mole fraction,
        % change the first value in the alpha equation to change dependence
        % steepness
%         alpha = exp(-6*(1-xsv)^3);
%         Dx = (Db(1)^(alpha*xsv)*Db(2)^(1-alpha*xsv));

        Cg2 = Cstar*Aieq_now;
        % difference between initial and equilibrium mole fraction
        if ti == 1;
            est = abs(Aieq_now-Aieq(1, 1));
        end
        
        % current particle-phase surface concentrations (mol/m^3 (particle)),
		% non-volatile first.  See starting_stuff for equation description
		Cptnv = ((1/Aieq_now)-1)/((M(1)/p(1))+(M(2)/p(2))*((1/Aieq_now)-1));
		Cptsv = (1-Cptnv*(M(2)/p(2)))/(M(1)/p(1));
		if Aieq_now == 0
			Cptnv = p(2)/M(2);
			Cptsv = 0.0;
		end 
		% convert to ug/m^3 (particle) from mol/m^3 (particle)
		Cptnv = Cptnv*M(2);
		Cptsv = Cptsv*M(1);
		Ceq = [Cptsv; Cptnv]; % equilibrium (surface) concentration (ug/m^3 (particle))
			
        % ratio of mole fraction at surface to that in the bulk 
		% assuming instantaneous gas-phase diffusion
        if esinc == 1
%             Cratstore(ti) = (Aieq_now)/ ...
%             (((Ai(1, ti)/M(1)))/((Ai(1, ti)/M(1))+(Ai(2, ti)/M(2))));
            Cratstore(ti) = (Aieq_now)- ...
            (((Ai(1, ti)/M(1)))/((Ai(1, ti)/M(1))+(Ai(2, ti)/M(2))));
        elseif esinc == 0 
            % mole fraction of semi-volatile in bulk
            Cratstore(ti) = (((Ai(1,ti)/M(1)))/ ...
                ((Ai(1,ti)/M(1))+(Ai(2,ti)/M(2)))); 
        end
        
        % get parameters for estimating scaling factor and estimate
        % scaling factor
            

        % **********************************
        % correction on
        Crat = Cratstore(ti);
%         if Ai(1, ti) == 0 
%             Crat = 10;
%         end
        
        Dscale = Dscale_fun(Crat, est, Db, esinc); 
        
%         Dscale = 1.0; % correction off when uncommented
        % **********************************
        Dx = Dx*Dscale; % corrected diffusion coefficient

        
        
        % overall gas-side mass transfer coefficient (m/s) (eq. 20)
        Kg = Kgcalc(Dg, Rpm0, Dx, Cstar, Aj, kci);
        
        % concentration difference factor in eqs. 29/31 (/s)
        k2 = 4*pi*Rpm0^2.*Nm*Kg;
        k2store(ti) = k2;
        % eqs. 31/32 for mass-transfer
        % average concentration in the particle-phase (ug/m^3 (air))
        % slow reaction case (neglect Q):
        errora = k2*(Cg2-xsv*Cstar)-(kci*Ca2(1))-...
                ((Ca2(1)-Ca0(1))/h);
        errorad = k2*(-(Cstar/(Aj*Vp0)))-kci-1/h; % differentiated form
        
        % loop until we get the error acceptably low
        while abs(errora)>(Ca2(1)/1e15) %|| abs(errorg)>(Cg2/1e14)

            % improve concentration at end of time step estimates (ug/m^3 
            % (air))
            Ca2(1) = Ca2(1)-(errora/errorad);

            % total number of moles in condensed-phase (mol)
            xT = (Ca2(1)/M(1)+Ca2(2)/M(2));
            % semi-volatile mole fraction in condensed-phase (fraction)
            xsv = (Ca2(1)/M(1))/xT;

            % revalue error
            errora = k2*(Cg2-xsv*Cstar)-(kci*Ca2(1))-...
                ((Ca2(1)-Ca0(1))/h);

        end

        % new mass of semi-volatile (ug/single particle)
        mugsv = Ca2(1)/Nm;
        
        % new volume of single particle (m^3)
        Vsingle = (mugsv/p(1))+(mugnv/p(2));
        % new particle bulk-average concentrations (ug/m^3 (particle))
        C1 = mugsv/Vsingle;
        Ai(1, ti+1) = C1;
        Cnv = mugnv/Vsingle; 
        Ai(2, ti+1) = Cnv;
        
        % new volume of all particles (m^3/m^3 (air))
        Vp0store(ti+1) = Vsingle*Nm;
        
        % new particle radius (m)
        Rpm0store(ti+1) = ((3*Vp0store(ti+1))/(4*pi*Nm))^(1/3);
        disp(ti)
        
        if Dcorr_it == 1
            % ****** use if iteratively finding D correction whilst stepping
            % through pde time steps ******************************************
%             while abs((pdeRpm0(ti+1)-Rpm0store(ti+1))/pdeRpm0(ti+1))>...
%                     abs(((pdeRpm0(end)-pdeRpm0(1))/pdeRpm0(1))/1e4)
            % first the most relevant pdeRpm
            pdergt0 = pdeRpm0(pdetimer <= pdetime(ti));
            
            if isempty(pdergt0) == 1
                pdeRpmnow = pdeRpm0(1);
            else
                pdeRpmnow = pdergt0(end);
            end
            count = 0;
            
            while abs(pdeRpmnow-Rpm0store(ti+1))/pdeRpm0(1)>1.0e-6
                
%                 if ti == 20 && count >-1
%                     disp('whoop')
%                     disp(pdeRpmnow)
%                     disp(Rpm0store(ti+1))
%                     disp((pdeRpmnow-Rpm0store(ti+1)))
%                     disp(Dx)
%                     disp(est)
%                     if count==10
%                         return
%                     end
%                 end
                % treatment depends on difference in mole fraction at 
                % particle surface (representing the gas-phase saturation 
                % ratio at t>0) and in the particle bulk initially 
                % (representing the gas-phase saturation ratio at t==0) 
                % (i.e deltaes value)
                if est<0.80 && est>0 && esinc == 0
                    Dx = Dx*(1-100.0*((pdeRpmnow-Rpm0store(ti+1))/pdeRpmnow));
                elseif est>=0.80 && esinc == 0
                    Dx = Dx*(1-2.0*((pdeRpmnow-Rpm0store(ti+1))/pdeRpmnow));
                elseif est>=0 && est<=0.05
                    Dx = Dx*(1+1000.0*((pdeRpmnow-Rpm0store(ti+1))/pdeRpmnow));
                elseif est>0.05
                    Dx = Dx*(1+10.0*((pdeRpmnow-Rpm0store(ti+1))/pdeRpmnow));
                end
                
                count = count+1;
                Kg = Kgcalc(Dg, Rpm0, Dx, Cstar, Aj, kci);
    %             k2 = Vp0*(3.0/Rpm0)*Kg;
                k2 = 4*pi*Rpm0^2.*Nm*Kg;
                errora = k2*(Cg2-(Ca2(1)/(Aj*Vp0))*Cstar)-(kci*Ca2(1))-...
                ((Ca2(1)-Ca0(1))/h);
                errorad = k2*(-(Cstar/(Aj*Vp0)))-kci-1/h;
                
                while abs(errora)>(Ca2(1)/1e4)%abs(errora)>(Ca2(1)/(eq_time*1e11))
                    
                    Ca2(1) = Ca2(1)-errora/errorad;
                    % total number of moles in condensed-phase (mol)
                    xT = (Ca2(1)/M(1)+Ca2(2)/M(2));
                    % semi-volatile mole fraction in condensed-phase (fraction)
                    xsv = (Ca2(1)/M(1))/xT;

                    % revalue error
                    errora = k2*(Cg2-xsv*Cstar)-(kci*Ca2(1))-...
                            ((Ca2(1)-Ca0(1))/h);
                end
                mugsv = Ca2(1)/Nm; Vsingle = (mugsv/p(1))+(mugnv/p(2));
                C1 = mugsv/Vsingle; Ai(1,ti+1) = C1; Cnv = mugnv/Vsingle;
                Ai(2,ti+1) = Cnv; Vp0store(ti+1) = Vsingle*Nm;
                Rpm0store(ti+1) = ((3*Vp0store(ti+1))/(4*pi*Nm))^(1/3);
%                 if ti<4;
%                     break
%                 end
            end
            
            
            % *************************************************************
        end
        Dscalestore(ti+1) = Dx/(Db(1)^xsv*Db(2)^xnv);
%         Dx_correct(2,ti) = Dx;
        % new concentration of condensed sv and nv (ug/m^3 (air))
        Ca2(1) = C1*Vp0store(ti+1);
        Ca2(2) = Cnv*Vp0store(ti+1);
        
        % store particle-phase concentration (ug/m^3 (air))
        Ca0store(:, ti+1) = Ca2;

        % density of the particle-phase (ug/m^3 (particle))
        Aj = sum(Ai(:, ti+1));

        times(ti+1) = times(ti)+h; % record time through simulation (s)
        ti = ti+1; % count on number of time steps
        while_var = ti;

        
    end % time loop

    eq_time = pdetime(end); % going off pde results (s)
    
end % analytical function

% -------------------------------------------------------------------------
% nested function for scaling factor of diffusion coefficient
function [Dscale] = Dscale_fun(Crat, est, Db, esinc)

    % increase in saturation ratio of the semi-volatile  
    if esinc == 1

         
        esx = [0.049999999 0.20 0.35 0.65 0.80 0.88];
        Draty = [1.000000001 1e-2 1e-4  1e-6 1e-8 1e-10 0.99999999e-12];


        p1 =  [1.50     1.55	1.60    1.65	1.70    1.75	1.80;
               1.75     1.80	1.85    1.90	1.95    2.00	2.05;
               2.00     2.00	2.00  	2.00	1.90	1.80    1.67;
               2.00     2.00	1.70    1.50	1.40    1.30	1.25;
               2.00     1.70	1.30	1.23	1.19	1.14	1.13;
               2.6      1.35	1.22     1.1 	1.08    1.07	1.13];

        p2 = [  150.0	185.0   228.0   285.0	352.0	450.0   580.0;
                30.0	40.0    57.0    77.0	105.0	135.0   180.0;
                15.0	24.0    36.0    51.0	56.0	61.0    61.0;
                6.0	    12.0    16.0    19.2	23.5    26.5	29.3;
                5.3 	10.2 	12.4	16.2	20.2	23.3	25.9;
                4.0     7.4      11.4 	16.0 	19.9 	22.6	25.3];
            
        p3 = [  0.7     0.7     0.7     0.7 	0.7    	0.7     0.7;
                0.4 	0.4     0.4     0.4 	0.4 	0.4     0.4;
                0.1	    0.1     0.1     0.1     0.1 	0.1     0.1;
                -0.3	-0.4    -0.4    -0.4	0.2    0.2 	    0.2;
                -2.3	-2.5	-1.5	-1.2	-0.8	-0.3	0.1;
                -2.5    -2.8    -1.5 	-1.5 	-1.5 	-1.5	0.0];
        
%         p1n(1) = interp1(esx, p1(:,1), est, 'linear');
%         p1n(2) = interp1(esx, p1(:,2), est, 'linear');
%         p1n(3) = interp1(esx, p1(:,3), est, 'linear');
%         p1n(4) = interp1(esx, p1(:,4), est, 'linear');
%         p1n(5) = interp1(esx, p1(:,5), est, 'linear');
%         p1n(6) = interp1(esx, p1(:,6), est, 'linear');
%         p1n(7) = interp1(esx, p1(:,7), est, 'linear');
% 
%         p2n(1) = interp1(esx, p2(:,1), est, 'linear');
%         p2n(2) = interp1(esx, p2(:,2), est, 'linear');
%         p2n(3) = interp1(esx, p2(:,3), est, 'linear');
%         p2n(4) = interp1(esx, p2(:,4), est, 'linear');
%         p2n(5) = interp1(esx, p2(:,5), est, 'linear');
%         p2n(6) = interp1(esx, p2(:,6), est, 'linear');
%         p2n(7) = interp1(esx, p2(:,7), est, 'linear');
%         
%         p3n(1) = interp1(esx, p3(:,1), est, 'linear');
%         p3n(2) = interp1(esx, p3(:,2), est, 'linear');
%         p3n(3) = interp1(esx, p3(:,3), est, 'linear');
%         p3n(4) = interp1(esx, p3(:,4), est, 'linear');
%         p3n(5) = interp1(esx, p3(:,5), est, 'linear');
%         p3n(6) = interp1(esx, p3(:,6), est, 'linear');
%         p3n(7) = interp1(esx, p3(:,7), est, 'linear');
%         
%         p1n = interp1(Draty, p1n, Db(2)/Db(1),'linear');
%         p2n = interp1(Draty, p2n, Db(2)/Db(1),'linear');
%         p3n = interp1(Draty, p3n, Db(2)/Db(1),'linear');

% interpolate with respect to the change in saturation ratio first
        p1n = zeros(1, 7);
        p2n = zeros(1, 7);
        p3n = zeros(1, 7);

        p1n(1) = (interp1((esx), (p1(:, 1)), (est), 'spline'));
        p1n(2) = (interp1((esx), (p1(:, 2)), (est), 'spline'));
        p1n(3) = (interp1((esx), (p1(:, 3)), (est), 'spline'));
        p1n(4) = (interp1((esx), (p1(:, 4)), (est), 'spline'));
        p1n(5) = (interp1((esx), (p1(:, 5)), (est), 'spline'));
        p1n(6) = (interp1((esx), (p1(:, 6)), (est), 'spline'));
        p1n(7) = (interp1((esx), (p1(:, 7)), (est), 'spline'));
        
        p2n(1) = 10.^(interp1(log10(esx), log10(p2(:, 1)), log10(est),'linear'));
        p2n(2) = 10.^(interp1(log10(esx), log10(p2(:, 2)), log10(est),'linear'));
        p2n(3) = 10.^(interp1(log10(esx), log10(p2(:, 3)), log10(est),'linear'));
        p2n(4) = 10.^(interp1(log10(esx), log10(p2(:, 4)), log10(est), 'linear'));
        p2n(5) = 10.^(interp1(log10(esx), log10(p2(:, 5)), log10(est), 'linear'));
        p2n(6) = 10.^(interp1(log10(esx), log10(p2(:, 6)), log10(est), 'linear'));
        p2n(7) = 10.^(interp1(log10(esx), log10(p2(:, 7)), log10(est), 'linear'));
        
        p3n(1) = (interp1(esx, (p3(:, 1)), est, 'spline'));
        p3n(2) = (interp1(esx, (p3(:, 2)), est, 'spline'));
        p3n(3) = (interp1(esx, (p3(:, 3)), est, 'spline'));
        p3n(4) = (interp1(esx, (p3(:, 4)), est, 'linear'));
        p3n(5) = (interp1(esx, (p3(:, 5)), est, 'spline'));
        p3n(6) = (interp1(esx, (p3(:, 6)), est, 'spline'));
        p3n(7) = (interp1(esx, (p3(:, 7)), est, 'spline'));
        
        
        
        % interpolate with respect to the ratio of self-diffusion 
        % coefficients
        if Db(2)/Db(1)>1e-3
            p1n = interp1(log10(Draty), p1n, log10(Db(2)/Db(1)), 'linear');
        else
            p1n = interp1(log10(Draty), p1n, log10(Db(2)/Db(1)), 'spline');
        end
        if est < 0.27
            if Db(2)/Db(1)<1e-8
                p2n = interp1(log10(Draty), p2n, log10(Db(2)/Db(1)), 'spline');
            else
                p2n = interp1(log10(Draty), (p2n), log10(Db(2.0)/Db(1)), 'linear');
            end
        elseif est >= 0.27 && est < 0.65
            if Db(2)/Db(1)<1e-4
                p2n = 10.^interp1(log10(Draty), log10(p2n+2.0), log10(Db(2)/Db(1)), 'linear')-2;
            else
                p2n = interp1(log10(Draty), (p2n), log10(Db(2)/Db(1)), 'linear');
            end
        elseif est >= 0.65
            if Db(2)/Db(1)>=1e-6 && Db(2)/Db(1)<=1e-4
                
                p2n = 10.^interp1(log10(Draty), log10(p2n), log10(Db(2)/Db(1)), 'linear');
            else
                
                p2n = 10.^(interp1(log10(Draty), log10(p2n), log10(Db(2)/Db(1)), 'linear'));
            end
        end
        
        p3n = (interp1(log10(Draty), (p3n), log10(Db(2)/Db(1)), 'spline'));
        

        Dscale = exp((Crat^p1n)*p2n)-p3n;
        
        % decrease in saturation ratio of the semi-volatile    
    elseif esinc == 0
        
        p1 = [2.81	2.86	2.92	3;
              3.23  3.53    3.46    2;
              3.65	4.4	    4.0	    2;
              5	    8	    5	    2;
              6	    11	    7	    1.9];

        p2 =  [8000	8000	8000	8000;
               350  300     100     -1.6;
               100	50	    -1.0	-1.6;
               23	12	    -1.0	-0.40;
               7	3.0	    0.55	-0.20];

        p3 = [0.40	0.42	0.40	0.42;
              0.32  0.41    0.50    0.52;
              0.25	0.40    0.58	0.62;
              0.00  0.50    0.67	0.76;
              -0.1	0.58	0.78	0.85];

        est = abs(est); % absolute value of difference in saturation ratio
        
        esx = [0.049999999 0.20 0.35 0.65 0.8800000000001 ];
        Draty = [1.000000001 1.0e-4  1.0e-8 0.99999999e-12];
        
        if est<esx(1) % for low values of est
            est = 0.05;
        end
        
        % interpolate with respect to the change in saturation ratio first
        p1n = zeros(1,4);
        p2n = zeros(1,4);
        p3n = zeros(1,4);

        p1n(1) = (interp1((esx), (p1(:,1)), (est), 'linear'));
        p1n(2) = (interp1((esx), (p1(:,2)), (est), 'linear'));
        p1n(3) = (interp1((esx), (p1(:,3)), (est), 'linear'));
        p1n(4) = (interp1((esx), (p1(:,4)), (est), 'linear'));
        
        p2n(1) = 10.^(interp1(log10(esx), log10(p2(:,1)), log10(est),'linear'));
        p2n(2) = 10.^(interp1(log10(esx), log10(p2(:,2)), log10(est),'linear'));
        p2n(3) = 10.^(interp1(log10(esx), log10(p2(:,3)+2.0), log10(est),'linear'))-2;
        p2n(4) = 10.^(interp1((esx), log10(p2(:,4)+2.0), (est), 'linear'))-2;
        
        p3n(1) = (interp1(esx, (p3(:,1)), est, 'linear'));
        p3n(2) = (interp1(esx, (p3(:,2)), est, 'linear'));
        p3n(3) = (interp1(esx, (p3(:,3)), est, 'linear'));
        p3n(4) = (interp1(esx, (p3(:,4)), est, 'linear'));
        
        
        
        % interpolate with respect to the ratio of self-diffusion 
        % coefficients

        p1n = interp1(log10(Draty),p1n,log10(Db(2)/Db(1)),'linear');
        if est < 0.27
            if Db(2)/Db(1)<1e-8
                p2n = 10.^interp1(log10(Draty),log10(p2n+2.0),log10(Db(2)/Db(1)),'spline')-2;
            else
                p2n = interp1(log10(Draty),(p2n),log10(Db(2.0)/Db(1)),'linear');
            end
        elseif est >= 0.27 && est < 0.65
            if Db(2)/Db(1)<1e-4
                p2n = 10.^interp1(log10(Draty),log10(p2n+2.0),log10(Db(2)/Db(1)),'spline')-2;
            else
                p2n = interp1(log10(Draty),(p2n),log10(Db(2.0)/Db(1)),'linear');
            end
        elseif est >= 0.65
            if Db(2)/Db(1)>=1e-6 && Db(2)/Db(1)<=1e-4
                
                p2n = 10.^interp1(log10(Draty), log10(p2n+1.1), log10(Db(2)/Db(1)), 'spline')-1.1;
            else
                
                p2n = 10.^(interp1(log10(Draty), log10(p2n+2.0), log10(Db(2)/Db(1)), 'linear'))-2;
            end
        end
        
        p3n = 10.^(interp1(log10(Draty), log10(p3n+2.0), log10(Db(2)/Db(1)), 'spline'))-2.0;
        
        
        Dscale = exp((Crat^p1n)*p2n)-p3n;
        
    end
    
end
% -------------------------------------------------------------------------
% nested function containing constants for input to the partitioning
% equation
function [kci, D0, Rpm, p, Cstar, Ai, Dg, Ceq, M] = ...
    starting_stuff(times, D, Cstar, Rpm, Asv0, Asvt)

    % ---------------------------------------------------------------------
    % inputs:
    % times - time points to solve diffusion at (s)
    % D - self-diffusion coefficients of components,semi-volatile first
    % (ug/m^3 (air))
    % Cstar - effective saturation concentration of semi-volatile (ug/m^3 
    % (air))
    % Rpm - radius (m)
    % Cg - gas-phase bulk concentration of semi-volatile (ug/m^3 (air))
    % Asv0 - initial mole fraction of semi-volatile in the particle-phase
    % Asvt - equilibrium mole fraction of semi-volatile in the particle-phase
    % ---------------------------------------------------------------------
    
    % particle-phase bulk concentration (ug/m^3 (particle)) (semi-volatile 
    % first) 
    Ai = zeros(2, length(times));
    
    % component densities semi-volatile (first) (ug/m^3 (particle))
    p = [1.0e12; 1.0e12]; 
    M = [1.0e8; 1.0e8]; % component molar masses (ug/mol), semi-volatile 1st
    kci = 0.0e-2; % pseudo-first-order reaction constant (/s)
    % self-diffusion coefficients of components, semi-volatile first, then
    % non-volatile (m^2/s)
    D0 = D; 
    
    % initial particle-phase concentrations (mol/m^3 (particle)),
    % non-volatile first.  These equations are derived from the mole
    % fraction equation (xsv = nsv/(nsv+nnv)) and the volume continuity
    % equation (nsv*(Msv/psv)+nnv*(Mnv/pnv)= total V (with total V then set
    % to 1m^3))
    Cp0nv = ((1/Asv0)-1)/((M(1)/p(1))+(M(2)/p(2))*((1/Asv0)-1));
    Cp0sv = (1-Cp0nv*(M(2)/p(2)))/(M(1)/p(1));
    
    if Asv0 == 0
        Cp0nv = p(2)/M(2);
        Cp0sv = 0.0;
    end 
    % convert to ug/m^3 (particle) from mol/m^3 (particle)
    Cp0nv = Cp0nv*M(2);
    Cp0sv = Cp0sv*M(1);
    Ai(:,1) = [Cp0sv; Cp0nv]; 
    
    % equilibrium particle-phase concentrations (mol/m^3 (particle)),
    % non-volatile first.  These equations are derived from the mole
    % fraction equation (xsv = nsv/(nsv+nnv)) and the volume continuity
    % equation (nsv*(Msv/psv)+nnv*(Mnv/pnv)= total V (with total V then set
    % to 1m^3))
    Cptnv = ((1/Asvt)-1)/((M(1)/p(1))+(M(2)/p(2))*((1/Asvt)-1));
    Cptsv = (1-Cptnv*(M(2)/p(2)))/(M(1)/p(1));
    if Asvt == 0
        Cptnv = p(2)/M(2);
        Cptsv = 0.0;
    end 
    % convert to ug/m^3 (particle) from mol/m^3 (particle)
    Cptnv = Cptnv*M(2);
    Cptsv = Cptsv*M(1);
    Ceq = [Cptsv; Cptnv]; % equilibrium concentrations (ug/m^3 (particle))

    % gas-phase diffusion coefficient (typical value taken from pp. 5159 of
    % Zaveri et al. (2014)) (m^2/s)
    Dg = 1.0e10;%Dg = 5.0e-6;
    
end
% -------------------------------------------------------------------------
% nested function for total particle volume
function [Vpt] = spat_arrays(Rpm,Nm)
    
    % total particle-phase volume (m^3/m^3(air))
    Vpt(1) = ((4.0/3.0)*pi*Rpm^3.0*Nm); 
end
% -------------------------------------------------------------------------
% nested function to calculate gas-side mass transfer coefficient (eq. 20 
% of Zaveri et al. 2014) 
function Kg = Kgcalc(~,Rp,Db,Cstar,Aj,kci)

    if kci == 0     
        % particle-side mass transfer coefficient (m/s) (eq. 24)
        kp = 5.0*(Db/Rp);     
    else
        qi = Rp*((kci/Db)^0.5);
        [Qi] = Qcalc(kci,Db,Rp);
        % particle-side mass transfer coefficient (m/s)  (eq. 23)
        kp = (Db/Rp)*((qi*coth(qi)-1.0)/(1.0-Qi)); 
    end
    
    % *********************************************************************
    % gas-side mass transfer coefficient (m/s)
    
%     % accommodation coefficient (pp. 5 of Zaveri et al. (2008)) (-)
%     alpha = 0.1;
%     % speed of semi-volatile vapour molecules (what MOSAIC gets when iv==12 in line 9146 
%     % of mosaic_box.25.f90) (m/s) (for a M of the semi-volatile of 100 g/mol
%     % and a temperature of 298.0 K (has units of cm/s in fortran MOSAIC))
%     speed = 25117.1943098747e-2;
%     % gas-phase mean free path (m) (has units of cm in fortran MOSAIC)
%     lambda = 3*Dg/speed;
%     Kn = lambda/Rp; % Knudsen number (-)
%     % transition regime correction factor
%     f = (0.75*alpha*(1+Kn))/(Kn*(1+Kn)+0.283*alpha*Kn+0.75*alpha);
% 
%     kg = (Dg/Rp)*f; % gas-side mass transfer coefficient (m/s)

    % ********************************************************************* 
    SR = (Cstar/(Aj)); % coefficient for kp term (eq. 20)
    
    % when gas-phase diffusion significant
%     if kci <1.e-2 % slow reaction
%         % overall gas-side mass transfer coefficient (m/s) (eq. 20)
%         Kg = 1/(1/(kg)+(1/(kp)*SR)); 
%     else % quick reaction
%         Kg = kg;
%     end
    % when gas-phase diffusion instantaneous (eq. 20)
    Kg = 1/(0+((1/kp)*SR));
    
end
% *************************************************************************
% nested function to calculate Q (eq. 8 of Zaveri et al. 2014) 
function [Qi] = Qcalc(kci,Dbi,Rpm)
    % ratio of the particle radius to the reacto-diffusive length (-)
    qi = Rpm*((kci/Dbi)^0.5);
    if kci>0
        Qi = 3.0*((qi*coth(qi)-1.0)/(qi^2.0)); % Q (-)
    else
        Qi = 1.0;
    end
end