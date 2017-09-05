% Copyright 2017 Simon O'Meara.  Program distributed under the terms of the 
% GNU General Public License

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



% Function to model simultaneous diffusion through a sphere using Fick's
% second law.  All components are allowed to diffuse, one is assumed to
% be non-volatile, while the other is semi-volatile.  Diffusion
% coefficients are allowed to vary simultaneously within the pde solver.

% The spatial mesh updates after every time step.  Made to be compatible
% with Dw_valuation_PDE_nonideal.  Double precision for variables input to
% the diffusion equation.  Like pde_ideal2 but with time step change turned
% on (e.g., in reponse to changes in radius and ratio of bulk to surface 
% concentration). 

%******************DOUBLE PRECISION, time adaptive VERSION*****************

% Simon O'Meara, University of Manchester, August 2015.
function [Zrec, timepde, efoldtime, rcpde] = pde_ideal4(D, ts, Dmethod2,...
    SN_part, Dp, RHt, kcsv)

    tic
    Temper = 298.15; % temperature (K)
    % ts, the integration time step (s) - initial guess at input should be 
    % low (e.g. 1.0e-14s)
    t = linspace(0, ts, 3);
    timepde = zeros(1, 1); % time steps results output at (s)
    rcpde = zeros(1, 1); % store for total particle radius (cm)
    % diffusion coefficients of components through pure solutions of the 
    % other component if using Dmethod 1 or 2 (see below) - water first
    % (cm^2/s)
    Db = D; 
    
    % molar mass and density of water and sucrose (in that order, from CRC 
    % online handbook)
    M = [1.0 1.0];%[18.015 342.296];
    p = [1.0 1.0];%1.5805];
    nchem = length(M); % number of components
    % the mutual diffusion coefficient estimation method to be used: 1) for
    % Vignes (1966) equation 10, 2) for Vignes (1966) equation 15, 3) for
    % Vrentas and Vrentas equation 35 (polymer-solvent systems)
    Dmethod = Dmethod2; 
   
    % ---------------------------------------------------------------------
    % spatial considerations
    
    % number of monolayers
    mlno = 1;
    % spatial arrays
    [clc,V,deltaz] = spat_arrays(mlno,SN_part,Dp);
    clcold = ((3*sum(V(1:SN_part)))/(4*pi))^(1/3);
    % ---------------------------------------------------------------------
    % empty concentration matrix, with initial concentrations inputted
    % (mol/cm^3).  Time in the first dimension, layers in the second one
    % and chemicals in the third (water first)
    [C,Zrec,Vcomp0] = concs(RHt,M,SN_part,timepde,p,V);
    
    % save first clc (cm)
    C(1,1:SN_part+1,3) = clc;
    % average bulk concentration of water at start (mol/cm^3)
    avZ2 = (sum(Zrec(1,SN_part+1:(SN_part*2-1),1)))/(sum(V(1:SN_part-1)));

    % initial difference in concentration between average and surface
    % (mol/cm^3)
    deltaZ0 = abs(C(1,SN_part+1, 1)-avZ2);
    tn = 0;
    Cratiores = ones(1, 1);
    % error tolerances for the pde
    options = odeset('RelTol', 1.0e-6, 'AbsTol', 1.0e-6);
    % remember the initial number of shells
    SN_part_orig = SN_part;
    
    % ---------------------------------------------------------------------
    while (deltaZ0/abs(C(1, SN_part_orig+1, 1)-avZ2))<exp(1) % time loop
        
        tn = tn+1; % time step number
        
        % values at start of time step
        clc0 = clc;
        C0 = C(tn,1:SN_part+1,1:length(M));
        SN_part0 = SN_part;
        V0 = V;
        Vcomp00 = Vcomp0;
        
        timepde(tn+1) = timepde(tn)+t(end);
        
        
        
        % diffusion
        sol = pdepe(2,@eqn,@ic,@bc,clc0,t,options,clc0,C0,Db,nchem,M, ...
            p,V0,Dmethod,Temper,timepde(tn+1),kcsv);
        
        
        % new concentrations (mol/cm^3)
        C(tn+1, 1:size(sol,2), 1:length(M)) = sol(end,:,:); 
        
        % -----------------------------------------------------------------
        % revalue spatial variables
        t0 = t;
        [C, V, clc, Zrec, Vcomp0, SN_part, t, sshellmark] = var_rev(C, V0, M, ...
            tn, p, Zrec, SN_part0, Vcomp00, deltaz, t0, clc0);
        
        % if we've had to change the time step, then repeat
        % diffusion with the new time step.  This first repeat can only
        % result from excess organic leaving the surface shell
        while t0(end) ~= t(end)
            
            timepde(tn+1) = timepde(tn)+t(end);
            sol = pdepe(2,@eqn,@ic,@bc,clc0,t,options,clc0,C0,Db, ...
                nchem,M,p,V0,Dmethod,Temper,timepde(tn+1),kcsv);
            % new concentrations
            C(tn+1,1:size(sol,2),1:length(M)) = sol(end,:,:); 
            t0 = t;
            [C,V,clc,Zrec,Vcomp0,SN_part,t,sshellmark] = var_rev(C,V0, ...
                M,tn,p,Zrec,SN_part0,Vcomp00,deltaz,t0,clc0);

        end
        
        % decrease if the change in radius from the pevious step is
        % more than 1% or increase if it's less than 0.1%
        clcnew = ((3*sum(V(1:SN_part)))/(4*pi))^(1/3);
        % flag for time step having been decreased due to insufficient
        % organic in the surface shell (cm^3)
        if sshellmark == 1
            sshellmark_perm = 1;
        else
            sshellmark_perm = 0;
        end
        % flag for decreasing time step
        t_decrease = 0;
        while (abs(clcnew-clcold)/clcold)*100.0>1.0
            disp('t decrease due to r change')
            t_decrease = 1;
            t = linspace(0,t(end)*0.1,3);

            timepde(tn+1) = timepde(tn)+t(end);
            
            sol = pdepe(2,@eqn,@ic,@bc,clc0,t,options,clc0,C0,Db,nchem,M,p,V0,Dmethod,Temper,timepde(tn+1),kcsv);
            C(tn+1,1:size(sol,2),1:length(M)) = sol(end,:,:); % new concentrations
            t0 = t;
            [C,V,clc,Zrec,Vcomp0,SN_part,t,sshellmark] = var_rev(C,V0,M,tn,p,Zrec,SN_part0,Vcomp00,deltaz,t0,clc0);
            clcnew = ((3*sum(V(1:SN_part)))/(4*pi))^(1/3);
            % if one surface shell marker is noted, then keep sshellmark at
            % a value of 1
            if sshellmark == 1
                sshellmark_perm = 1;
            end
        end
        
        % if the ratio of bulk to surface concentration changes too
        % greatly then decrease time step
        avZ2 = (sum(Zrec(tn+1,SN_part+1:(SN_part*2-1),1)))/(sum(V(1:SN_part-1)));
        while abs((deltaZ0/abs(C(1,SN_part_orig+1,1)-avZ2))-Cratiores(tn))>1.0e-1
            disp('t decrease')
            t_decrease = 1;
            t = linspace(0,t(end)*0.9,3);
            timepde(tn+1) = timepde(tn)+t(end);
            
            sol = pdepe(2,@eqn,@ic,@bc,clc0,t,options,clc0,C0,Db, ...
                nchem,M,p,V0,Dmethod,Temper,timepde(tn+1),kcsv);
            % new concentrations (mol/cm^3)
            C(tn+1,1:size(sol,2),1:length(M)) = sol(end,:,:); 
            t0 = t;
            [C,V,clc,Zrec,Vcomp0,SN_part,t,sshellmark] = var_rev(C,V0,M,tn,p,Zrec,SN_part0,Vcomp00,deltaz,t0,clc0);
            clcnew = ((3*sum(V(1:SN_part)))/(4*pi))^(1/3);
            % if one surface shell marker is noted, then keep sshellmark at
            % a value of 1
            if sshellmark == 1
                sshellmark_perm = 1;
            end
            % ratio of the surface to bulk concentration
            avZ2 = (sum(Zrec(tn+1,SN_part+1:(SN_part*2-1),1)))/(sum(V(1:SN_part-1)));
            
        end
        
        if t_decrease==0 && sshellmark_perm==0
            
            while (abs(clcnew-clcold)/clcold)*100.0<0.005
                disp('t increase')
                t = linspace(0, t(end)*1.1, 3);
                timepde(tn+1) = timepde(tn)+t(end);
            
                sol = pdepe(2, @eqn, @ic, @bc, clc0, t, options, clc0, C0, Db, nchem, M, p, V0, Dmethod, Temper, timepde(tn+1), kcsv);
                C(tn+1, 1:size(sol, 2), 1:length(M)) = sol(end, :, :); % new concentrations
                t0 = t;
                [C, V, clc, Zrec, Vcomp0, SN_part, t, sshellmark] = ...
                    var_rev(C, V0, M, tn, p, Zrec, SN_part0, Vcomp00, deltaz, ...
                    t0, clc0);

                % average concentration of water (mol/cm^3)
                clcnew = ((3*sum(V(1:SN_part)))/(4*pi))^(1/3);
                % if the surface shell marker is noted, then break.
                if sshellmark == 1
                    break
                end
                % if the ratio of bulk to surface concentration changes too
                % greatly then stop increasing the time step
                avZ2 = (sum(Zrec(tn+1,SN_part+1:(SN_part*2-1),1)))/(sum(V(1:SN_part-1)));
                if abs((deltaZ0/abs(C(1,SN_part_orig+1,1)-avZ2))-Cratiores(tn))>1.0e-2
                    break
                end
                
            end
        end
        
        % check if surface shell needs splitting due to excess growth
        % *****************************************************************
        if ((3*sum(V(1:SN_part)))/(4*pi))^(1/3)- ...
                ((3*sum(V(1:SN_part-1)))/(4*pi))^(1/3)>deltaz*2.0
            [SN_part,Zrec,C,V,clc,Vcomp0] = sl_split(SN_part,C,V,M,p, ...
                tn,Zrec,deltaz,Vcomp0);
        end
        
        % *****************************************************************
        % average concentration of water (mol/cm^3)
        avZ2 = (sum(Zrec(tn+1,SN_part+1:(SN_part*2-1),1)))/...
            (sum(V(1:SN_part-1)));
        clcnew = ((3*sum(V(1:SN_part)))/(4*pi))^(1/3);
        % revalue the last time step's ratio of difference in 
        % concentration
        clcold = clcnew;
        Cratiores(tn+1) = (deltaZ0/abs(C(1,SN_part_orig+1,1)-avZ2));
        rcpde(tn) = clcnew;
        disp(deltaZ0/abs(C(1,SN_part_orig+1,1)-avZ2))
        disp(tn)

    end % time loop
    efoldtime = interp1(Cratiores,timepde,exp(1)); % e-folding time (s)
    disp('efoldtime = ')
    disp(efoldtime)
    toc % display computer run time (s)
    % ---------------------------------------------------------------------
    % nested function for splitting the surface shell in two
    function [SN_part,Zrec,C,V,clc,Vcomp0] = ...
            sl_split(SN_part0,C,V0,M,p,tn,Zrec,deltaz,Vcomp0)
        
        % shell end points (cm)
        rc = zeros(1,length(V0));
        % first shell (cm)
        rc(1) = (((3.0*sum(V0(1)))/(4.0*pi))^(1/3));
        for ir2 = 2:length(V0)
            rc(ir2) = ((3.0*sum(V0(1:ir2)))/(4.0*pi))^(1/3);
        end
        clear ir2
        % volume of surface shell needed to return to the original width
        % (cm^3)
        Vreq = ((4.0/3.0)*pi)*(rc(end)^3.0-(rc(end)-deltaz)^3.0);
        
        % if the near-surface shell is already quite large then create a
        % new shell
        if rc(SN_part0-1)-rc(SN_part0-2)>(rc(SN_part0-1)/SN_part0)
            
            % volume of the new near-surface shell (cm^3)
            Vnsl = V0(SN_part0)-Vreq;
            disp(V0(SN_part0).*1e-6)
            % new shell volume array (cm^3)
            V = horzcat(V0(1:SN_part0-1),Vnsl,Vreq);
            % new concentration array (mol/cm^3)
            C(tn+1,1:SN_part0+2,:) = horzcat(C(tn+1,1:SN_part0+1,:),C(tn+1,SN_part0+1,:));
            % increase number of shells
            SN_part = SN_part0+1;
            % new number of moles - concentrations are constant (mol)
            N = zeros(length(M),length(V));
            for ic2 = 1:length(M)
                N(ic2,:) = C(tn+1,2:SN_part+1,ic2).*V;
            end
            clear ic2
            % volumes of each component per shell (cm^3), water first
            Vcomp0 = zeros(2,SN_part);
            Vcomp0(1,:) = N(1,:)*(M(1)/p(1));
            Vcomp0(2,:) = N(2,:)*(M(2)/p(2));
        % otherwise if the near-surface shell is thin, just transfer
        % material from the surface shell to this
        else
            SN_part = SN_part0;
            % ratio of organic to total volume
            Vproporg2 = Vcomp0(2,SN_part0)/(sum(Vcomp0(:,SN_part0)));
            % new volumes of each component (cm^3)
            Vcomp0(1,SN_part0-1:SN_part0) = horzcat(Vcomp0(1,SN_part0-1)+(V0(SN_part0)-Vreq)*(1-Vproporg2),Vreq*(1-Vproporg2));
            Vcomp0(2,SN_part0-1:SN_part0) = horzcat(Vcomp0(2,SN_part0-1)+(V0(SN_part0)-Vreq)*(Vproporg2),(Vreq*(Vproporg2)));
            % new total volumes of shells (cm^3)
            V = V0;
            V(SN_part0-1) = V0(SN_part0-1)+(V0(SN_part0)-Vreq);
            V(SN_part0) = Vreq;

            % new number of moles - water first (mol)
            N = zeros(length(M),length(V));
            for ic2 = 1:length(M)
                N(ic2,:) = (Vcomp0(ic2,1:SN_part0).*(p(ic2)/M(ic2)));
            end
            clear ic2
            % new concentrations, water first (mol/cm^3)
            for ic2 = 1:2
                C(tn+1,2:SN_part0+1,ic2)=N(ic2,1:SN_part0)./(V);
            end
            clear ic2 
        end
        
        % new shell centre coordinates (cm).
        clc = zeros(1,length(V)+1);
        % new shell end points (cm)
        rc = zeros(1,length(V));
        % first shell (cm)
        clc(2) = (((3*sum(V(1)))/(4*pi))^(1/3))/2.0;
        rc(1) = (((3*sum(V(1)))/(4*pi))^(1/3));
        for ir2 = 2:length(V)
            clc(ir2+1) = ((3*sum(V(1:ir2-1)))/(4*pi))^(1/3)+...
                (((3*sum(V(1:ir2)))/(4*pi))^(1/3)-((3*sum(V(1:ir2-1)))/(4*pi))^(1/3))/2.0;
            rc(ir2) = ((3*sum(V(1:ir2)))/(4*pi))^(1/3); 
        end
        clear ir2
        % fix concentrations at the centre to equal those in the core shell
        % (mol/cm^3)
        C(tn+1,1,:) = C(tn+1,2,:);
        % save clc (cm)
        C(tn+1,1:SN_part+1,3) = clc;
        % save new number of moles and shell volumes
        Zrec(tn+1,1:size(N,2)*2,1) = horzcat(N(2,:),N(1,:));
        Zrec(tn+1,1:length(V),2) = V;
        
    end
    % ---------------------------------------------------------------------
    % nested function for revaluing variables after diffusion
    function [C,V,clc,Zrec,Vcomp0,SN_part,t,sshellmark] = var_rev(C,V,M, ...
            tn,p,Zrec,SN_part,Vcomp0,deltaz,t,clc)
        
        sshellmark = 0;
        % new number of moles
        N = zeros(length(M),length(V));
        for ic2 = 1:length(M)
            N(ic2,:) = C(tn+1, 2:SN_part+1, ic2).*V;
        end
        clear ic2
        % volumes of each component per shell (cm^3), water first
        Vcompn = zeros(2,SN_part);
        Vcompn(1, :) = N(1, :)*(M(1)/p(1));
        Vcompn(2, :) = N(2, :)*(M(2)/p(2));
        
        % volume of non-volatile in surface shell (cm^3)
        slvolnv = Vcomp0(2, end)+sum(Vcomp0(2, 1:SN_part-1))-sum(Vcompn(2, 1:SN_part-1));
        % concentration difference between the near-surface and surface
        % shell.
        Cdiffnss = (C(tn+1,SN_part,1)-(C(tn+1,SN_part+1,1)));
        Cdiffcs = (C(tn+1,1,1)-(C(tn+1,SN_part+1,1)));
        % if the near-surface shell's concentration has gone past the
        % equilibrium concentration then decrease time step and repeat
        % diffusion.
        if Cdiffcs>0 && Cdiffnss<0 || Cdiffcs<0 && Cdiffnss>0
            t = t.*0.5;
            sshellmark = 1;
            return
        end
        if slvolnv<0
            t = t.*0.5;
            sshellmark = 1;
            return
        end
        
        % corresponding volume of semi-volatile required to maintain
        % equilibrium (cm^3).
        slvolsv = slvolnv*(Vcomp0(1, end)/Vcomp0(2, end));
        % revalue new volumes (cm^3)
        Vcomp0(2, 1:SN_part-1) = Vcompn(2, 1:SN_part-1);
        Vcomp0(1, 1:SN_part-1) = Vcompn(1, 1:SN_part-1);
        Vcomp0(2, end) = slvolnv;
        Vcomp0(1, end) = slvolsv;
        
        % new total volumes per shell (cm^3)
        V = sum(Vcomp0, 1);
        % new shell centre coordinates (cm)
        clc = zeros(1, length(V)+1);
        % new shell end points (cm)
        rc = zeros(1, length(V));
        % first shell (cm)
        clc(2) = (((3*sum(V(1)))/(4*pi))^(1/3))/2.0;
        rc(1) = (((3*sum(V(1)))/(4*pi))^(1/3));
        for ir2 = 2:length(V)
            clc(ir2+1) = ((3*sum(V(1:ir2-1)))/(4*pi))^(1/3)+...
                (((3*sum(V(1:ir2)))/(4*pi))^(1/3)-((3*sum(V(1:ir2-1)))/(4*pi))^(1/3))/2.0;
            rc(ir2) = ((3*sum(V(1:ir2)))/(4*pi))^(1/3); 
        end
        clear ir2
        
        % new concentration of shells (mol/cm^3)
        for ic2 = 1:length(M)
            C(tn+1,2:SN_part+1,ic2) = (Vcomp0(ic2,:)*(p(ic2)/M(ic2)))./V;
            N(ic2,:) = Vcomp0(ic2,:)*(p(ic2)/M(ic2));
        end
        clear ic2
        
        
        
        % check that the surface shell hasn't shrunk excessively
        if (rc(end)-rc(SN_part-1))<(deltaz/10.0)

            % if it has then first find the volume required to return it to
            % a width of deltaz (cm^3).
            Vreq = (((4/3)*pi*(rc(SN_part-1)+deltaz)^3))- ...
                (((4/3)*pi*(rc(SN_part-1)^3)));
            Vfrac = 1.0;
            % Additional volume needed to achieve this volume (cm^3).
            Vreq1=(Vreq-(sum(V(:))-sum(V(1:SN_part-1))))/Vfrac;
            % If this volume is greater than the volume of the near-surface
            % layer, then reduce the volume.
            while Vreq1>V(SN_part-1)
                Vfrac=Vfrac*2.0;
                % How much needs to be added to achieve this volume (cm^3)?
                Vreq1=(Vreq-(sum(V(:))-sum(V(1:SN_part-1))))/Vfrac;
            end
            % The proportion of this that should be organic (fraction).
            Vproporg=(Vcomp0(2,end)/(Vcomp0(1,end)+Vcomp0(2,end)));
            % The volume of organic and water needed (cm^3).
            Vorgreq=Vreq1*Vproporg;
            Vwreq=Vreq1*(1-Vproporg);
            % Get the current volume ratio of water to organic in the
            % near-surface layer.
            Vratworg=(Vcomp0(1,SN_part-1)/Vcomp0(2,SN_part-1));
            % Take the required organic from the next layer in (cm^3), and take the
            % required amount of water from this layer to maintain the
            % concentration there.
            Vcomp0(2,SN_part-1)=Vcomp0(2,SN_part-1)-Vorgreq;
            Vcomp0(1,SN_part-1)=Vcomp0(1,SN_part-1)-(Vorgreq*Vratworg);
            % Add the required volumes to the surface layer, and any extra
            % water is presumed to evaporate (cm^3).
            Vcomp0(2,end)=Vcomp0(2,end)+Vorgreq;
            Vcomp0(1,end)=Vcomp0(1,end)+Vwreq;
            
            % New total volumes per layer (cm^3).
            V=sum(Vcomp0,1);
            % New layer centre coordinates (cm).
            clc=zeros(1,length(V)+1);
            % New layer end points (cm)
            rc=zeros(1,length(V));
            % First layer (cm).
            clc(2)=(((3*sum(V(1)))/(4*pi))^(1/3))/2.0;
            rc(1)=(((3*sum(V(1)))/(4*pi))^(1/3));
            for ir2=2:length(V)
                clc(ir2+1)=((3*sum(V(1:ir2-1)))/(4*pi))^(1/3)+...
                    (((3*sum(V(1:ir2)))/(4*pi))^(1/3)-((3*sum(V(1:ir2-1)))/(4*pi))^(1/3))/2.0;
                rc(ir2)=((3*sum(V(1:ir2)))/(4*pi))^(1/3); 
            end
            clear ir2
            N=zeros(length(M),length(V));
            % New concentration of layers (mol/cm^3).
            for ic2=1:length(M)
                C(tn+1,2:SN_part+1,ic2)=(Vcomp0(ic2,:)*(p(ic2)/M(ic2)))./V;
                N(ic2,:)=Vcomp0(ic2,:)*(p(ic2)/M(ic2));
            end
            clear ic2
            % -------------------------------------------------------------
            % if the near-surface shell's width falls below the width of
            % the outer shell, then combine the two shells into one new
            % outer shell

            if (rc(SN_part-1)-rc(SN_part-2))<(deltaz/10.0)
                % the current ratio of water to organic volume in the
                % surface shell
                Vratworg = (Vcomp0(1,SN_part)/Vcomp0(2,SN_part));
                % Add any organic in the near surface layer to the outer 
                % layer (cm^3), then add the required volume of water to
                % maintain the equilibrium concentration (cm^3).
                Vcomp0(2,end) = Vcomp0(2,end)+Vcomp0(2,SN_part-1);
                Vcomp0(1,end) = Vcomp0(2,end)*Vratworg;
                % now remove the near-surface shell
                Vcomp0 = horzcat(Vcomp0(:,1:SN_part-2),Vcomp0(:,end));
                SN_part = SN_part-1;
                
                % new total volumes per shell (cm^3)
                V = sum(Vcomp0,1);
                % new shell centre coordinates (cm)
                clc = zeros(1,length(V)+1);
                % new shell end points (cm)
                rc = zeros(1,length(V));
                % first shell (cm)
                clc(2) = (((3*sum(V(1)))/(4*pi))^(1/3))/2.0;
                rc(1) = (((3*sum(V(1)))/(4*pi))^(1/3));
                for ir2 = 2:length(V)
                    clc(ir2+1) = ((3*sum(V(1:ir2-1)))/(4*pi))^(1/3)+ ...
                        (((3*sum(V(1:ir2)))/(4*pi))^(1/3)- ...
                        ((3*sum(V(1:ir2-1)))/(4*pi))^(1/3))/2.0;
                    rc(ir2)=((3*sum(V(1:ir2)))/(4*pi))^(1/3); 
                end
                clear ir2
                % new concentration of shells (mol/cm^3)
                N = zeros(length(M),length(V));
                for ic2=1:length(M)
                    C(tn+1,2:SN_part+1,ic2) = (Vcomp0(ic2,:)* ...
                        (p(ic2)/M(ic2)))./V;
                    N(ic2,:) = Vcomp0(ic2,:)*(p(ic2)/M(ic2));
                end
                clear ic2
            end
            
        end
        
        % fix concentrations at the centre to equal those in the core shell
        % (mol/cm^3).
        C(tn+1,1,1:length(M)) = C(tn+1, 2, 1:length(M));
        % save clc (cm).
        C(tn+1,1:SN_part+1,3) = clc;
        % save new number of moles and shell volumes
        Zrec(tn+1, 1:size(N,2)*2, 1) = horzcat(N(2,:), N(1,:));
        Zrec(tn+1, 1:length(V), 2) = V;
        
    end
    
    % ---------------------------------------------------------------------
    % nested function to value spatial arrays
    function [clc,V0,deltaz] = spat_arrays(mlno,SN_part,Dp)
        
        % diameter of water (cm), as chosen by Zobrist et al. (2011).
%         deltaz=3.0e-8; %(Dp/2.0)/(SN_part);  % ;
        deltaz = (Dp/2.0)*(1.0e-3);
        r0 = (Dp/2.0)-(mean(deltaz)*mlno); % initial radius of bulk (cm)
        rc = zeros(SN_part,1); % cumulative layer widths (cm)
        delta = zeros(SN_part,1); % radial width of layers (cm)

        for ir2 = 1:SN_part-mlno
            rc(ir2) = (r0/(SN_part-mlno))*ir2;
            delta(ir2) = (r0/(SN_part-mlno));
        end
        clear ir2
        % 2 monolayers
        if mlno == 2
            rc(SN_part-1:SN_part) = [r0+deltaz;r0+deltaz*2];
            delta(SN_part-1:SN_part) = [deltaz;deltaz];
        elseif mlno==1
        % one monolayer
            rc(end) = (Dp/2.0);
            delta(end) = deltaz;
        end
        % cumulative layer centre points (cm)
        clc = zeros(SN_part+1,1);
        clc(3:end) = rc(1:length(rc)-1)+delta(2:end)/2.0;
        clc(2) = rc(1)/2.0;
        % volume of individual layers (cm^3)
        V0 = horzcat((4/3)*pi*(rc(1).')^3,((4/3)*pi)*((rc(2:end).').^3 ...
            -(rc(1:length(rc)-1).').^3));
    end
    % ---------------------------------------------------------------------
    % nested function to calculate initial abundances
    function [C,Zrec,Vcomp0] = concs(RH,M,SN_part,timepde,p,V)
        
        % in the ideal case the mole fraction of water equals the
        % fractional relative humidity.
        % Rearrange the volume continuity equation, and the mole fraction
        % equation to find equations for concentrations:
        Zorg0 = 1/((M(2)/p(2))+(RH(1)/(1-RH(1)))*(M(1)/p(1)));
        Zw0 = Zorg0*(RH(1)/(1-RH(1)));
        Zorgeq = 1/((M(2)/p(2))+(RH(2)/(1-RH(2)))*(M(1)/p(1)));
        Zweq = Zorgeq*(RH(2)/(1-RH(2)));
        
        % empty concentration results matrix (mol/cm^3) (water first).
        % Make the number of columns in the C matrix double the inital
        % number of shells + colum for C at the centre, to allow for the 
        % possibility of particle growth and therefore shell addition
        C = zeros(length(timepde),(SN_part+1)*2,length(M)+1);
        C(1,1:SN_part+1,1) = horzcat(ones(1,SN_part)*Zw0,Zweq);
        C(1,1:SN_part+1,2) = horzcat(ones(1,SN_part)*Zorg0,Zorgeq);
        % store number of moles of each component (mol) and shell volume 
        % (cm^3)
        Zrec=zeros(length(timepde),SN_part*2,2);
        Zrec(1,:,1)=horzcat((C(1,2:SN_part+1,2).*V),(C(1,2:SN_part+1,1).*V));
        Zrec(1,1:SN_part,2)=V;
        
        % Volumes of each component per shell (cm^3), water first.
        Vcomp0=zeros(2,SN_part);
        Vcomp0(1,:)=Zrec(1,SN_part+1:end,1).*(M(1)/p(1));
        Vcomp0(2,:)=Zrec(1,1:SN_part,1).*(M(2)/p(2));
        clear Zorg0 Zorgeq Zw0 Zweq
    end
    % ---------------------------------------------------------------------
    function [c,f,s] = eqn(~,~,u,DuDx,~,~,Db,~,M,p,~,Dmethod,Temper,time,kcsv)      
        x = [u(1)/sum(u),u(2)/sum(u)]; % mole fraction
        % diffusion coefficient (cm^2/s)
        Dpde = diff_coeff(Dmethod, Db, x, M, p, Temper, time, u);
        % c and f equations used previously, suitable for constant
        % diffusion coefficients, but not those that vary as a function of
        % radius
%         c = 1./Dpde(:); % diffusion coefficient (cm^2/s)
%         f = DuDx; % concentration vs. radius gradient
        c = [1;1]; % diffusion coefficient (cm^2/s)
        f = Dpde(:).*DuDx; % concentration vs. radius gradient
        % source/sink term: semi-volatile first, then non-volatile
        s = [kcsv*u(1);0]; 
    end
    % initial conditions
    function [C0] = ic(x,x2,Carray,~,~,~,~,~,~,~,~,~,~) 
        rn = (x2(:)==x);
        % withdraw preceding concentrations (mol/cm^3)
        C0 = squeeze(squeeze(Carray(1,rn,:))); 
    end
    % boundary conditions
    function [pl,ql,pr,qr]=bc(~,ul,~,ur,~,~,Carray,~,~,~,~,~,~,~,~,~,~) 
        pl=ul;
        ql=[0;0];
        pr=ur-squeeze(Carray(1,end,:)); % right b.c. (mol/cm^3).
        qr=[0;0];
    end
end