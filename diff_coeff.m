% function to calculate the mutual diffusion coefficient
% Simon O'Meara, University of Manchester, December 2014

function [D] = diff_coeff(method,Db,x,M,p,Temper,time,u)
    
    % call upon the correct estimation, given the method
    if method==4
        D = Vignes_eq10(Db,x);
    elseif method == 2
        D = Vignes_eq15(Db,x);
    elseif method == 5
        D = Vrentas_eq35(Db,x,M,p);
    elseif method == 1
        D2 = constant(Db);
    elseif method == 3
        D = sig(Db,x);
    elseif method == 6
        D2 = diff_coeffETH(Temper,u,M,time);
    elseif method == 3.5
        D = sigETHres(Db,x);
    end
        
    % duplicate the mutual diffusion coefficient, so the pde solver can see
    % one value per component
    if method>1 && method<6
        D = horzcat(D,D);
    end

    if method==1 || method==6
        D = D2;
    end
    % ---------------------------------------------------------------------
    % nested function for equation 10 of Vignes (1966) - linear dependence
    % on mole fraction
    function [D] = Vignes_eq10(Db,x)
        % multiply the mole fraction of each component by its
        % self-diffusion coefficient
        D = sum(Db.*x);
    end
    % ---------------------------------------------------------------------
    % nested function for equation 15 of Vignes (1966) - logarithmic
    % dependence on mole fraction
    function [D] = Vignes_eq15(Db,x)
        D = prod(Db.^x);
    end
    % ---------------------------------------------------------------------
    % nested function for equation 35 of Vrentas and Vrentas (2000) -
    % solvent-polymer system
    function [D] = Vrentas_eq35(Db,x,M,p)
        W = Db(1)./Db(2); % Eq. 31
        % volume fraction, semi-volatile
        vf = (x(1)*V*(M(1)/p(1)))/sum((x(:)*V).'.*(M./p)); 
        D = Db(1).*((1+W+vf(1)*(W-1))/... % Eq. 35.
            (1+W-vf(1)*(W-1)));
    end
    % ---------------------------------------------------------------------
    % nested function for maintaining the diffusion coefficients of
    % individual components
    function [D2] = constant(Db)
        D2 = Db;
    end
    % ---------------------------------------------------------------------
    % nested function for sigmoidal dependence on mole fraction
    function [Dbi] = sig(Db,x)
        % note, semi-volatile should be the first x value in x
        
        % correction parameter coefficients
        C = -3.105;
        Dpar = 3.3;
    	cp = exp((1.0-x(1))^2.0*(C+3.0*Dpar-4.0*Dpar*(1.0-x(1))));
    	% exponents
    	exponen = [x(1)*cp,(1-x(1)*cp)];
        Dbi = prod(Db.^exponen);
        
    end
    % ---------------------------------------------------------------------
    function [D] = sigETHres(Db,x)
        
        Dw = Db(1);
        Dorg = Db(2);
        xw = x(1);
        C = -3.0;
        D = -5.0;
        D = (Dw^(xw*(exp(((1-xw)^2)*(D+3*C-4*C*(1-xw))))))*...
            (Dorg^(1-xw*((exp(((1-xw)^2)*(D+3*C-4*C*(1-xw)))))));
        
    end
    % ---------------------------------------------------------------------
    % nested function for Zobrist (2011) diffusion coefficient
    function [D] = diff_coeffETH(Temper,Z1,M,time)
    	% diffusion coefficient of semi-volatile per boundary (cm^2/s)
        % weight fraction of solute in this shell
        ws=(Z1(2).*M(2))./(Z1(2).*M(2)+Z1(1).*M(1));
        % temperature at this time (K)
        if time<=Temper(2,1)
            Temper3=Temper(1,1);
        else
            Temper3=Temper(1,2);
        end
        
		% water activity in this shell
        [aw] = aw_calc(ws,Temper3);
        
        % diffusion coefficients (m^2/s)
        a=7+0.175*(1-46.46.*(1-aw));
        b=262.867*(1+10.53.*(1-aw)-0.3.*(1-aw).^2.0);
        T0=127.9*(1+0.4514.*(1-aw)-0.5.*(1-aw).^1.7);
        Dw=10.^(-(a+b./(290.0-T0))); % semi-volatile
        Dorg=10.^(-(a+b./(292.0-T0))); % non-volatile
        % convert diffusion coefficient to cm^2/s (first element) and
        % include D of non-volatile (2nd element)
        D=[Dw.*1.0e4 Dorg.*1.0e4];
        
    end
    % ---------------------------------------------------------------------
    % nested function to replicate the water activity equation (eq. 10) of 
    % Zobrist et al. 2011
    function [aw]=aw_calc(ws,Temper3)
        
        % water activity parameters values
        a = -1;
        b = -0.99721;
        c = 0.13599;
        d = 0.001688;
        e = -0.005151;
        f = 0.009607;
        g = -0.006142;
        T0 = 298.15;
       
        aw = ((1+a.*ws)./(1+b.*ws+c.*ws.^2))+(Temper3-T0)*(d.*ws+e.*ws.^2+f.*ws.^3+g.*ws.^4);
        
        if length(aw) == 1
            if ws == 1.0 ||aw<0.0
                aw = 0.0;
            elseif ws == 0.0 || aw>1.0;
                aw = 1.0;
            end
        else
            searchwshi = ws(:)>=1.0;
            searchawlo = aw(:)<0.0;
            searchwslo = ws(:)<=0.0;
            searchawhi = aw(:)>=1.0;
    %         ensure limits are kept to
            aw(searchwshi) = 0.0;
            aw(searchawlo) = 0.0;

            aw(searchwslo) = 1.0;
            aw(searchawhi) = 1.0;
        end
    end
    % ---------------------------------------------------------------------
    
                                    
end