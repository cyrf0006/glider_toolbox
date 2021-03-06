function dbar = bar2dbar(bar)
%BAR2DBAR  Convert pressure from bars to decibars.
%
%  Syntax:
%    DBAR = BAR2DBAR(BAR)
%
%  Description:
%    DBAR = BAR2DBAR(BAR) converts pressure readings in array BAR from bars to
%    decibars by multiplying by 10.
%
%  Notes:
%    This is simply a convenience function to call the conversion with an
%    explicit name.
%
%  Examples:
%    dbar = bar2dbar(bar)
%
%  Authors:
%    Joan Pau Beltran  <joanpau.beltran@socib.cat>

%  Copyright (C) 2014-2015
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

  error(nargchk(nargin, 1, 1, 'struct'));

  dbar = bar * 10;

end
