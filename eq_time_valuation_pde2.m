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

% function to find e-folding times for PDE, remember to change file names
% to reflect what's saved inside them
function [timepde, rcpde] = eq_time_valuation_pde2()

    % ---------------------------------------------------------------------
    % inputs
    
    % Dorg (low to high) - self-diffusion coefficients of organics (cm^2/s)
    % SN_part - number of shells (one for each organic self-diffusion 
    % coefficient value), see .xls document
    % 'Model_shell_res_sensitivity' in 'Diffusion' folder
    % Dp (high to low) - initial particle diameters (cm)
    % Treat_no_min - D dependency method number minimum
    % Treat_no_max - D dependency method number maximum
    % RH - instant change in saturation ratio (fraction) (e.g. [0.01 0.90]
    % for 1-90%)
    % kcsv - condensed-phase reaction rate (/s)
    % note, the output saved files should have their names below manually
    % altered to reflect the saturation ratio change
    % ---------------------------------------------------------------------
    
    Dorg = 1.0e-10;
    SN_part = 100;
    Dp = 2.0e-2; % (cm)
    Treat_no_min = 2;
    Treat_no_max =2;
    RH = [0.00 0.88];
    kcsv = 0.0;
    
    % number of organic self-diffusion coefficients to loop through
    Dorg_no = length(Dorg); 
    
    % results array
    timepde = zeros(length(Dp), Dorg_no, (Treat_no_max-Treat_no_min+1));
    
    % efold time first guess, particle diameters represented along rows and
    % D along columns (s)
    ts0 = 9e-3;
    
    % organic diffusion coefficient loop
    for Dorgi = 1:Dorg_no
        
        % non-volatile self-diffusion coefficient (cm^2/s)
        Dorgcon = Dorg(Dorgi);
        % shell number corresponds to organic self-diffusion coefficient
        SN_partn = SN_part(Dorgi);  

        % particle diameter loop (cm)
        for Dpi = 1:length(Dp)
            
            % diameter (cm)
            Dpn = Dp(Dpi);
            % count on D dependency number
            Treat_noi = 0;
            % D dependency loop
            for Treat_no = Treat_no_min:Treat_no_max
                
                % count on D dependency
                Treat_noi = Treat_noi+1;
                % display inputs to diffusion simulation
                
                disp('Dp = ')
                disp(Dpn)
                disp('Dmethod=')
                disp(Treat_no)

                ts = ts0(Dpi, Dorgi)/1000;
                disp('ts = ')
                disp(ts)
                % in case we are using the constant D treatment
                if Treat_no == 1
                    Dw = Dorg;   
                else
                    % diffusion coefficient of water (cm^2/s) 
                    % (Scala et al. 2000)
                    Dw = 1.0e-4; 
%                     
                end
                D = [Dw, Dorgcon];
                disp('D (Dsv first) = ')
                disp(D)
                
                % simulate diffusion
                [Zrec, timepde, efoldtime, rcpde] = pde_ideal4(D, ts, Treat_no, ...
                    SN_partn, Dpn, RH, kcsv);
                % pde times and radii
                
                save('efold_times_logD_e00_et88_r4_Dnv14_Dsv8_pde_1e_100shell.mat', 'timepde', 'rcpde')
                
                % ---------------------------------------------------------
            
            end % D dependency loop
        end % diameter loop    
    end % Dorg loop
    
end % function