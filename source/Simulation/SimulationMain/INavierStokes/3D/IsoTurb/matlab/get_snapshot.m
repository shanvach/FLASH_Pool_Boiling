function [DATA]=get_snapshot(file)

MDIM = 3;

%---Define precision format.
iformat = 'int32';
rformat = 'float64';

%---Read file
fid = fopen(file, 'r');

% Read Time:
temp = fread(fid, 1, iformat);
vars = fread(fid, [1], rformat);
temp = fread(fid, 1, iformat);

DATA.time = vars;

% Read istep,blockCount,NXB,NYB,NZB:
temp = fread(fid, 1, iformat);
vars = fread(fid, [5], iformat);
temp = fread(fid, 1, iformat);

DATA.istep =      vars(1);
DATA.blockCount = vars(2);
NXB = vars(3); NYB = vars(4); NZB = vars(5); 

DATA.NXB =            NXB;
DATA.NYB =            NYB;
DATA.NZB =            NZB;

% Allocate:
DATA.del   = zeros(MDIM,DATA.blockCount);
DATA.coord = zeros(MDIM,DATA.blockCount);
DATA.bsize = zeros(MDIM,DATA.blockCount);

DATA.zc    = zeros(NZB,DATA.blockCount);
DATA.uavg    = zeros(NZB,DATA.blockCount);
DATA.vavg    = zeros(NZB,DATA.blockCount);
DATA.wavg    = zeros(NZB,DATA.blockCount);
DATA.pavg    = zeros(NZB,DATA.blockCount);

DATA.omgx2   = zeros(NZB,DATA.blockCount);
DATA.omgy2   = zeros(NZB,DATA.blockCount);
DATA.omgz2   = zeros(NZB,DATA.blockCount);

DATA.q2      = zeros(NZB,DATA.blockCount);
DATA.urms2   = zeros(NZB,DATA.blockCount);
DATA.vrms2   = zeros(NZB,DATA.blockCount);
DATA.wrms2   = zeros(NZB,DATA.blockCount);
DATA.prms2   = zeros(NZB,DATA.blockCount);

DATA.eps_nu  = zeros(NZB,DATA.blockCount);
DATA.lambdaf = zeros(NZB,DATA.blockCount);
DATA.lambdag = zeros(NZB,DATA.blockCount);

DATA.uu  = zeros(NZB,DATA.blockCount);
DATA.vv  = zeros(NZB,DATA.blockCount);
DATA.ww  = zeros(NZB,DATA.blockCount);
DATA.uv  = zeros(NZB,DATA.blockCount);
DATA.uw  = zeros(NZB,DATA.blockCount);
DATA.vw  = zeros(NZB,DATA.blockCount);

% Now Start Block loop
for lb=1:DATA.blockCount
    
    % Read deltas, coord, bsize:
    temp = fread(fid, 1, iformat);
    vars = fread(fid, [MDIM], rformat);
    temp = fread(fid, 1, iformat);
    DATA.del(1:MDIM,lb)  = vars(1:MDIM);
    
    temp = fread(fid, 1, iformat);
    vars = fread(fid, [MDIM], rformat);
    temp = fread(fid, 1, iformat);
    DATA.coord(1:MDIM,lb)= vars(1:MDIM);
    
    temp = fread(fid, 1, iformat);
    vars = fread(fid, [MDIM], rformat);
    temp = fread(fid, 1, iformat);
    DATA.bsize(1:MDIM,lb)= vars(1:MDIM);

    % Read zc location of block lb:
    temp = fread(fid, 1, iformat);
    vars = fread(fid, [NZB], rformat);
    temp = fread(fid, 1, iformat);
    DATA.zc(1:NZB,lb) = vars(1:NZB)';
    
    % Read mean velocities and pressure:
    temp = fread(fid, 1, iformat);
    vars = fread(fid, [NZB], rformat);
    temp = fread(fid, 1, iformat);
    DATA.uavg(1:NZB,lb) = vars(1:NZB)';
    
    temp = fread(fid, 1, iformat);
    vars = fread(fid, [NZB], rformat);
    temp = fread(fid, 1, iformat);
    DATA.vavg(1:NZB,lb) = vars(1:NZB)';
    
    temp = fread(fid, 1, iformat);
    vars = fread(fid, [NZB], rformat);
    temp = fread(fid, 1, iformat);
    DATA.wavg(1:NZB,lb) = vars(1:NZB)';
    
    temp = fread(fid, 1, iformat);
    vars = fread(fid, [NZB], rformat);
    temp = fread(fid, 1, iformat);
    DATA.pavg(1:NZB,lb) = vars(1:NZB)';
    
    % Read vorticities squared:
    temp = fread(fid, 1, iformat);
    vars = fread(fid, [NZB], rformat);
    temp = fread(fid, 1, iformat);
    DATA.omgx2(1:NZB,lb) = vars(1:NZB)';
    
    temp = fread(fid, 1, iformat);
    vars = fread(fid, [NZB], rformat);
    temp = fread(fid, 1, iformat);
    DATA.omgy2(1:NZB,lb) = vars(1:NZB)';

    temp = fread(fid, 1, iformat);
    vars = fread(fid, [NZB], rformat);
    temp = fread(fid, 1, iformat);
    DATA.omgz2(1:NZB,lb) = vars(1:NZB)';

    % Read q^2 and rms fluctuations:
    temp = fread(fid, 1, iformat);
    vars = fread(fid, [NZB], rformat);
    temp = fread(fid, 1, iformat);
    DATA.q2(1:NZB,lb) = vars(1:NZB)';
    
    temp = fread(fid, 1, iformat);
    vars = fread(fid, [NZB], rformat);
    temp = fread(fid, 1, iformat);
    DATA.urms2(1:NZB,lb) = vars(1:NZB)';
    
    temp = fread(fid, 1, iformat);
    vars = fread(fid, [NZB], rformat);
    temp = fread(fid, 1, iformat);
    DATA.vrms2(1:NZB,lb) = vars(1:NZB)';

    temp = fread(fid, 1, iformat);
    vars = fread(fid, [NZB], rformat);
    temp = fread(fid, 1, iformat);
    DATA.wrms2(1:NZB,lb) = vars(1:NZB)';

    temp = fread(fid, 1, iformat);
    vars = fread(fid, [NZB], rformat);
    temp = fread(fid, 1, iformat);
    DATA.prms2(1:NZB,lb) = vars(1:NZB)';

    
    % Read Dissipation and microscales:
    temp = fread(fid, 1, iformat);
    vars = fread(fid, [NZB], rformat);
    temp = fread(fid, 1, iformat);
    DATA.eps_nu(1:NZB,lb) = vars(1:NZB)';
    
    temp = fread(fid, 1, iformat);
    vars = fread(fid, [NZB], rformat);
    temp = fread(fid, 1, iformat);
    DATA.lambdaf(1:NZB,lb) = vars(1:NZB)';

    temp = fread(fid, 1, iformat);
    vars = fread(fid, [NZB], rformat);
    temp = fread(fid, 1, iformat);
    DATA.lambdag(1:NZB,lb) = vars(1:NZB)';

    % Read Reynolds stresses:
    temp = fread(fid, 1, iformat);
    vars = fread(fid, [NZB], rformat);
    temp = fread(fid, 1, iformat);
    DATA.uu(1:NZB,lb) = vars(1:NZB)';
    
    temp = fread(fid, 1, iformat);
    vars = fread(fid, [NZB], rformat);
    temp = fread(fid, 1, iformat);
    DATA.vv(1:NZB,lb) = vars(1:NZB)';

    temp = fread(fid, 1, iformat);
    vars = fread(fid, [NZB], rformat);
    temp = fread(fid, 1, iformat);
    DATA.ww(1:NZB,lb) = vars(1:NZB)';

    temp = fread(fid, 1, iformat);
    vars = fread(fid, [NZB], rformat);
    temp = fread(fid, 1, iformat);
    DATA.uv(1:NZB,lb) = vars(1:NZB)';

    temp = fread(fid, 1, iformat);
    vars = fread(fid, [NZB], rformat);
    temp = fread(fid, 1, iformat);
    DATA.uw(1:NZB,lb) = vars(1:NZB)';

    temp = fread(fid, 1, iformat);
    vars = fread(fid, [NZB], rformat);
    temp = fread(fid, 1, iformat);
    DATA.vw(1:NZB,lb) = vars(1:NZB)';

    
end

fclose(fid);

return