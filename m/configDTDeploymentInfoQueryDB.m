function [sql_query, deployment_fields] = configDTDeploymentInfoQueryDB()
%CONFIGDTDEPLOYMENTINFOQUERYDB  Configure the query to retrieve delayed time glider deployment information.
%
%  Syntax:
%    [SQL_QUERY, DEPLOYMENT_FIELDS] = CONFIGDTDEPLOYMENTINFOQUERYDB()
%
%  Description:
%    [SQL_QUERY, DEPLOYMENT_FIELDS] = CONFIGDTDEPLOYMENTINFOQUERYDB() should 
%    return the SQL query to retrieve the information about glider deployments
%    to be processed in delayed time. String SQL_QUERY is the query to execute. 
%    The mapping between deployment fields and data base table columns is given 
%    by the string cell array DEPLOYMENT_FIELDS. Deployment fields are described 
%    in GETDEPLOYMENTINFODB.
%
%  Notes:
%    Edit this file filling in the field mapping of your data base and the
%    query that returns that fields for each deployment.
%
%  Examples:
%    [sql_query, deployment_fields] = configDTDeploymentInfoQueryDB()
%
%  See also:
%    GETDEPLOYMENTINFODB
%
%  Authors:
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

  error(nargchk(0, 0, nargin, 'struct'));

  % Select the identifiers of deployments to process.
  deployment_ids = {'1' '2'};
  
  % Select the deployment fields.
  % First column is deployment field
  % Second column is column in data base table.
  fields_map = {
    'deployment_id'          'deployment_id'
    'deployment_name'        'deployment_name'
    'deployment_start'       'deployment_initial_date'
    'deployment_end'         'deployment_end_date'
    'glider_name'            'platform_name'
    'glider_serial'          'instrument_serial'
    'glider_model'           'instrument_model'
  };

  deployment_fields = fields_map(:,1)';
  db_fields = fields_map(:,2)';

  % Build the query.
  db_fields_str = [sprintf('%s, ', db_fields{1:end-1}) db_fields{end}];

  sql_query = ['select ' db_fields_str ...
               '  from instrumentation.deployment' ...
               '  inner join instrumentation.instrument' ...
               '    on (deployment_instrument_id=instrument_id)' ...
               '  inner join instrumentation.instrument_type' ...
               '    on (instrument_instrument_type_id=instrument_type_id)' ...
               '  inner join instrumentation.instrument_platform' ...
               '    on (instrument_platform_instrument_id=instrument_id and (instrument_platform_installation_date < deployment_end_date or deployment_end_date is null) and (deployment_initial_date < instrument_platform_uninstallation_date or instrument_platform_uninstallation_date is null))' ...
               '  inner join instrumentation.platform' ...
               '    on (instrument_platform_platform_id = platform_id)' ...
               '  inner join instrumentation.institution' ...
               '    on (deployment_institution_id=institution_id)' ...
               '  where (deployment_id in (' deployment_ids '));'];

end
