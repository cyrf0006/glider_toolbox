function data = getDeploymentInfoLocal(dplFile, varargin)
%GETDEPLOYMENTINFOLOCAL  Get deployment information from deployment local file.
%
%  Syntax:
%    DATA = GETDEPLOYMENTINFOLocal(DPLFILE)
%    DATA = GETDEPLOYMENTINFOLocal(DPFILE, OPTIONS)
%    DATA = GETDEPLOYMENTINFOLocal(DPLFILE, OPT1, VAL1, ...)
%
%  Description:
%    DATA = GETDEPLOYMENTINFOLocal(DPLFILE) reads a
%    deployment file where meta_data on the deployment are manually
%    edited (DPLFILE). Here is a minimalistic example of such file:
%     glider_name = {'gliderName'} 
%     deployment_id = {78}
%     deployment_name = {'MarsPlanet'}
%     deployment_start = {'2015-10-28 09:36:55'}
%     deployment_end = {'2015-11-10 09:12:40'}
%     glider_serial: '001'
%     glider_model: 'SeaExplorer'
%     project = {'MARS2015'}
%     abstract = {'Some description'}
%     project_url = {'http://myproject.org/'}
%     principal_investigator = {'Pi Name'}
%     principal_investigator_email = {'pi@glider.edu'}
%     author = {'Glider Respo'}
%     author_email = {'respo@glider.edu'}
%   
%   This function also reads the SeaExplorer configuration file
%   'seapayload.cfg' (PAYLOADFILE). The function returns a struct
%   DATA with fields given by corresponding columns in the input
%   files. 
%
%    DATA = GETDEPLOYMENTINFOLocal(DPLFILE, PAYLOADFILE, OPTIONS) and 
%    DATA = GETDEPLOYMENTINFOLocal(DPLFILE, PAYLOADFILE, OPT1, VAL1, ...) accept the 
%    following options given in key-value pairs OPT1, VAL1... or in struct
%    OPTIONS with field names as option keys and field values as option values:
%      FIELDS: database column renaming.
%        String cell array with alternative field names for output structure.
%        It should have the same number of elements than selected columns.
%        If empty, no renaming is done and column names are used as field names.
%        Default value: [] (do not rename columns)
%      TIME_FIELDS: timestamp fields.
%        String cell array with the name of the output fields to be converted to
%        from timestamp string to serial date number.
%        Default value: {'deployment_start' 'deployment_end'}
%      TIME_FORMAT: timestamp field format.
%        String with the format of the timestamp columns returned by the query.
%        Default value: 'yyyy-mm-dd HH:MM:SS' (ISO 8601 format)
%      CFGFILE*: Extra information from the configuration file.
%        String giving the name of deployment/config file of type
%        seapayload.cfg. This file should exist since it is used to
%        configure the gliders prior deployment.    
%        Default value: 'none' (in this case the info will  not be
%        present in final struct). If present, the info will be
%        stored in meta.headers.config
%         *Specific to SeaExplorer gliders
%    
%    The returned struct DATA should have the following fields to be considered 
%    a deployment structure:
%      DEPLOYMENT_ID: deployment identifier (invariant over time).
%      DEPLOYMENT_NAME: deployment name (may eventually change).
%      DEPLOYMENT_START: deployment start date (see note on time format).
%      DEPLOYMENT_END: deployment end date (see note on time format).
%      GLIDER_NAME: glider platform name (present in Slocum file names).
%      GLIDER_SERIAL: glider serial code (present in Seaglider file names).
%      GLIDER_MODEL: glider model name (like Slocum G1, Slocum G2, Seaglider).
%    The returned structure may include other fields, which are considered to be
%    global deployment attributes by functions generating final products like
%    GENERATEOUTPUTNETCDF.
%
%  Notes:
%    Time columns selected in the query should be returned as UTC timestamp
%    strings in ISO 8601 format ('yyyy-mm-dd HH:MM:SS') or other format accepted
%    by DATENUM, and are converted to serial date number format. 
%    Null entries are set to invalid (NaN).
%
%  Example (here 'Info.dpl' and 'seapayload.cfg' MUST be the local folder):
%    deployment_list = getDeploymentInfoLocal('Info.dpl', 'cfgfile', 'seapayload.cfg');
%
%  See also:
%    GENERATEOUTPUTNETCDF
%    FETCH
%    DATENUM
%
%  Authors:
%    Frederic Cyr  <Frederic.Cyr@mio.osupytheas.fr>
%    Joan Pau Beltran  <joanpau.beltran@socib.cat>

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

    error(nargchk(1, 10, nargin, 'struct'));

    
    %% Set options and default values.
    options.fields = {};
    options.time_fields = {'deployment_start' 'deployment_end'};
    options.time_format = 'yyyy-mm-dd HH:MM:SS';
    options.cfgfile = 'none';

    
    %% Parse optional arguments.
    % Get option key-value pairs in any accepted call signature.
    argopts = varargin;
    if isscalar(argopts) && isstruct(argopts{1})
        % Options passed as a single option struct argument:
        % field names are option keys and field values are option values.
        opt_key_list = fieldnames(argopts{1});
        opt_val_list = struct2cell(argopts{1});
    elseif mod(numel(argopts), 2) == 0
        % Options passed as key-value argument pairs.
        opt_key_list = argopts(1:2:end);
        opt_val_list = argopts(2:2:end);
    else
        error('glider_toolbox:getDeploymentInfoLocal:InvalidOptions', ...
              'Invalid optional arguments (neither key-value pairs nor struct).');
    end
    % Overwrite default options with values given in extra arguments.
    for opt_idx = 1:numel(opt_key_list)
        opt = lower(opt_key_list{opt_idx});
        val = opt_val_list{opt_idx};
        if isfield(options, opt)
            options.(opt) = val;
        else
            error('glider_toolbox:getDeploymentInfoLocal:InvalidOption', ...
                  'Invalid option: %s.', opt);
        end
    end
    
    %% Retrieve deployment info from local file as a structure.  
    fid=fopen(dplFile,'rt');
    configFile = {};
    while ~feof(fid) % Dump whole file in cell array
        str=fgetl(fid); 
        configFile = [configFile; {str}];
    end
    data = struct();
    for i = 1:numel(configFile)
        str=configFile{i}; % read line-by-line

        if ~isempty(str)     
            % if a comment line
            if str(1) == '#' ||  str(1) == '%'
                continue
            elseif ~isempty(strfind(str, '=')); 
                command = ['data.' str ' ;'];
                eval(command);
            end
        else
            continue
        end
    end
 
    %% Retrieve data from payload file as a structure (if asked).
    % (In a later version, it should read and store content)
    if strcmp(options.cfgfile, 'none') ~= 1
        fid=fopen(options.cfgfile,'rt');
        if fid ~=-1
            config = {};
            while ~feof(fid) % Dump whole file in cell array
                str=fgetl(fid); 
                config = [config; str];
            end    
            data.configuration_file = sprintf('%s\n', config{:});
        else
            disp(['[WARNING!] No configuration file found. Field left empty']);
        end 
    end
         
    %% Convert to cell array for postprocessing.
    % MATLAB is not consistent when the DataReturnFormat is structure.
    % If no rows match the selected query, an empty array is returned instead.
    if isstruct(data)
        if isempty(options.fields)
            fields = fieldnames(data);
        else
            fields = cellstr(options.fields(:));
        end
        data = struct2cell(data);
        for i = 1:size(data, 1)
            if isnumeric(data{i})
                data{i} = num2cell(data{i});
            end
        end
        data = horzcat(data{:});
    else
        fields = cellstr(options.fields);
        data = cell(0, numel(fields));
    end


    %% Convert time fields from timestamp string to serial date number.
    time_format = options.time_format;
    time_fields = cellstr(options.time_fields);
    time_field_columns = ismember(fields, time_fields);
    time_data = data(:,time_field_columns);
   
    if iscellstr(time_data)
        % DATENUM does not handle empty date string cell arrays properly.
        time_data_null = strcmp('null', time_data);
        if any(~time_data_null(:))
            time_data(~time_data_null) = ...
                num2cell(datenum(time_data(~time_data_null), time_format));
        end
        time_data(time_data_null) = {NaN};
    else
        error('glider_toolbox:db_tools:TimeFieldError', ...
              'Wrong time data type (not a timestamp string).');
    end
    data(:, time_field_columns) = time_data;
        

    %% Convert back to structure array with new field names.
    data = cell2struct(data, fields, 2);
    
end
