function theResult = nc_load_bfm(theNetCDFFile, varargin)

% nc_load_bfm -- Load NetCDF variables from BFM output files
%  nc_load_bfm('theNetCDFFile', 'var1', 'var2', ...) loads the
%   given variables of 'theNetCDFFile' into the Matlab
%   workspace of the "caller" of this routine.  If no names
%   are given, all variables are loaded.  The names of the
%   loaded variables are returned or assigned to "ans".
%   Attributes corresponding to 'units' and 'long_name' are 
%   loaded in a structure ending with the suffix '_att'.
%
%  nc_load_bfm('theNetCDFFile', 'prefix', 'str' 'var1', 'var2', ...)
%   prepend the prefix str to all the variable names.
%   'prefix', 'pref' and 'pre' are synonyms.
%
%  nc_load_bfm('theNetCDFFile', 'suffix', 'str' 'var1', 'var2', ...)
%   append the suffix str to all the variable names.
%   'suffix', 'suff' and 'suf' are synonyms.
%
% Examples:
%  nc_load_bfm('file.nc','suffix','_new','var1')
%   creates a variable named 'var1_new' in the workspace.
%
%  nc_load_bfm('file.nc','prefix','single_','var1')
%   creates a variable named 'single_var1' in the workspace.
%
% nc_load_bfm requires snctools installed on your system
%
% Modified by P.G. Fogli and M. Vichi, CMCC Bologna, Italy.

error(nargchk(1,1000,nargin));
%
result = [];
if nargout > 0, theResult = result; end
%
if (~ischar(theNetCDFFile))
    error(['Not a valid file name: ''' theNetCDFFile ''' is not a string.']);
end
%
if (exist(theNetCDFFile,'file')~=2)
    error(['Unable to read file ''' theNetCDFFile ''': No such file.']);
end
%
pref=[];
suf=[];
point='_';
vars=cell(0);
%
i=1;
j=1;
while i<=nargin-1
   switch varargin{i}
       case {'prefix','pref','pre'}
           if ~ischar(varargin{i+1})
               error(['Bad value for property: ''' varargin{i} '''.']);
           elseif (any(strcmp({'suffix';'suff';'suf';'point';'.'},varargin{i+1})))
               error(['Bad value for property ''' varargin{i} '''.']);
           else
               pref=varargin{i+1};
               i=i+2;
           end
       case {'suffix','suff','suf'}
           if ~ischar(varargin{i+1})
               error(['Bad value for property: ''' varargin{i} '''.']);
           elseif (any(strcmp({'prefix';'pref';'pre';'point';'.'},varargin{i+1})))
               error(['Bad value for property ''' varargin{i} '''.']);
           else
               suf=varargin{i+1};
               i=i+2;
           end
       case {'point','.'}
           if ~ischar(varargin{i+1})
               error(['Bad value for property: ''' varargin{i} '''.']);
           elseif (any(strcmp({'suffix';'suff';'suf';'prefix';'pref';'pre'},varargin{i+1})))
               error(['Bad value for property ''' varargin{i} '''.']);
           else
               point=varargin{i+1};
               i=i+2;
           end
       otherwise
           vars{j}=varargin{i};
           j=j+1;
           i=i+1;
   end
end
% Get all the names if none is given
if (isempty(vars))
    info = nc_info(theNetCDFFile);
    NV = numel(info.Dataset);
    vars = cell(1,NV);
    for n=1:NV
        vars{n} = info.Dataset(n).Name;
    end
end
%
for i = 1:numel(vars)
    thedata = nc_varget(theNetCDFFile,vars{i});
    if (isa(thedata,'single'))
        thedata = double(thedata);
    end
    varname = regexprep(vars{i},{'\.','\(','\)','-'},{point,'_','','_'});
    assignin('caller',[pref varname suf], ...
        thedata);
    vinfo = nc_getvarinfo(theNetCDFFile,vars{i});
    assignin('caller', [pref varname suf '_att'], vinfo.Attribute)
end
%
result = vars;
% change time to the matlab format
time=nc_varget(theNetCDFFile,'time');
vinfo=nc_getvarinfo(theNetCDFFile,'time');
time_att=vinfo.Attribute;
start=datenum(time_att.Value(15:end));
time_num=start+time/86400.;
assignin('caller','time_num',time_num)
%
if nargout > 0
   theResult=result;
end
%
return
