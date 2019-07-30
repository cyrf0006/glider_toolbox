function cleaning_options = configDataCleaning()
%CONFIGDATACLEANING  Configure glider data cleaning.
%   WORKING COPY!!!!
%  Syntax:
%    CLEANING_OPTIONS = CONFIGDATACLEANING()
%
%  Description:
%    CLEANING_OPTIONS = CONFIGDATACLEANING() should return a struct setting the 
%    options for glider data cleaning as needed by the function GRIDGLIDERDATA.
%
%  Examples:
%    cleaning_options = configDataCleaning()
%
%  See also:
%    CLEANGLIDERDATA
%
%  Authors:
%    Frederic Cyr  <Frederic.Cyr@mio.osupytheas.fr>


%  Copyright (C) 2013-2015
%  ICTS SOCIB - Servei d'observacio i prediccio costaner de les Illes Balears
%  <http://www.socib.es>
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  error(nargchk(0, 0, nargin, 'struct'));

  cleaning_options = struct();
  
  cleaning_options.time_list = {'time'};
  
  cleaning_options.variable_list = {
    'conductivity'
    'temperature'
    'pressure'
    'chlorophyll'
    'turbidity'
    'cdom'
    'oxygen_concentration'
    'oxygen_saturation'
    'conductivity_corrected_thermal'
    'temperature_corrected_thermal'
    'salinity'
    'density'
    'salinity_corrected_thermal'
    'density_corrected_thermal'
    'mfl_v1'
    'mfl_v2'
    'mfl_v3'
    'mfl_v4'
    'methane_volt'
    'mfl_methane_scaled'
    'mets_temp_volt'
    'mets_temp_scaled'
    'backscatter_700'
    'oxygen_frequency'                  
  };

  
  cleaning_options.lowpass_default = {0.05}; % in Hz
  cleaning_options.std_threshold_default = {1}; % 
  
  cleaning_options.non_default_parameters = {}; % 

end
