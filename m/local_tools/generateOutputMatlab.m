function mat = generateOutputMatlab(filename, data, meta, data_grid, meta_grid, deployment)
%GENERATEOUTPUTMATLAB  Generate Matlab output for glider deployment data.
%
%  Syntax:
%    MAT = GENERATEOUTPUTNETCDF(FILENAME, DATA, META, DATA_GRID, META_GRID, DEPLOYMENT)
%
%  Notes:
%    Usually input data is the output of LOADSLOCUMDATA or LOADSEAGLIDERDATA, 
%    PROCESSGLIDERDATA or GRIDGLIDERDATA.
%
%  Examples:
%    mat = generateOutputNetCDF(filename, data, meta, data_grid, meta_grid, deployment)
%
%  See also:
%    SAVENC
%    POSIXTIME2UTC
%    POSIXTIME
%    NMEA2DEG
%    LOADSLOCUMDATA
%    LOADSEAGLIDERDATA
%    PREPROCESSGLIDERDATA
%    PROCESSGLIDERDATA
%    GRIDGLIDERDATA
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

  
    error(nargchk(6, 6, nargin, 'struct'));
  
  
    %% Fill structure 
    s = struct();
    s.deployment = deployment;
    s.data_raw = data;
    s.meta_raw = meta;
    s.data_grid = data_grid;
    s.meta_grid = meta_grid;

    %% Rename struct
    struct_name = deployment.deployment_name;
    command = sprintf('%s = s', struct_name);
    eval(command);

    %% Create base directory of target file if needed.
    % This seems to be the best way to check if a relative path points to
    % an existing directory (EXIST checks for existance in the whole load path).
    [file_dir, ~, ~] = fileparts(filename);
    [status, attrout] = fileattrib(file_dir);
    if ~status
        [success, message] = mkdir(file_dir);
        if ~success
            error('glider_toolbox:generateOutputMatlab:MatlabDirectoryError', ...
                  'Could not create directory %s: %s.', file_dir, message);
        end
    elseif ~attrout.directory
        error('glider_toolbox:generateOutputMatlab:MatlabDirectoryError', ...
              'Not a directory: %s.', attrout.Name);
    end
    
    
    %% Generate the file.
    command = sprintf('save(filename, ''%s'')', struct_name);
    eval(command)
  
    
    %% Return the absolute name of the generated file.
    [status, attrout, ~] = fileattrib(filename);
    if status==0
        % We should never get here (if Matlab creation succeed, file must exist).
        error('glider_toolbox:generateOutputMatlab:MatlabFileError', ...
              'Matlab file generation succeed but problems with output file %s: %s.', ...
              filename, attrout);
    end
    mat = attrout.Name;

end
