function ConvertNCtoText
% TOM
clc;
clear;

% INPUTS
fold='/Users/lovato/GIT/bfm/run/standalone/co2test';
fname='CMIP5_Historical_GHG_1765_2005.nc';
target='CO2';

% Read Netcdf
ncid = netcdf.open(strcat(fold,'/',fname), 'NC_NOWRITE');
[~,nvars,~,unlimdimid] = netcdf.inq(ncid);

% Retrieve timeline
time = netcdf.getVar(ncid,unlimdimid);

% Retrieve target variable
for varid = 1 , nvars;
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncid,varid);
    if strcmp(varname,target) == 1
        data = netcdf.getVar(ncid,varid);
        display ('ok');
    end;
end
netcdf.close(ncid)

% Manual set of new timeline
newtime = datenum((1764:2004),1,1);
timest = datestr(newtime,'yyyy-mm-dd HH:MM');

% save to ascii file
fid = fopen('CMIP5_Historical_GHG_1765_2005.dat','w+');

fprintf(fid,'"%16s" , "%15s" \n','Time',target);

for i = 1 : length(data);
fprintf(fid,'"%16s" , %15.8f \n', timest(i,:),data(i));
end
fclose(fid);

display ('FATTOOOOOO');

end