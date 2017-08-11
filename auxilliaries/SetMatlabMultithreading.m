warning('off','MATLAB:maxNumCompThreads:Deprecated'); %Suppress deprecation warning for 'maxNumCompThreads()'
%-----PREPARE TO DO PARALLEL CALCULATION-----
%Resources available on cluster node
numProc_cluster = str2num(getenv('PBS_NUM_PPN'));
numProc_private = 4;

%Find current poolsize
if any(strcmpi(version('-release'),{'2010a','2010b','2011a','2011b','2012a','2012b','2013a'}))
    CurMatPool = matlabpool('size');
elseif any(strcmpi(version('-release'),{'2013b','2014a','2014b','2015a','2015b','2016a','2016b','2017a','2017b','2018a','2018b','2019a','2019b'}))
    poolobj = gcp('nocreate'); %Will be empty if no pool
    if isempty(poolobj); CurMatPool = 0; else CurMatPool = poolobj.NumWorkers; end;
else
    error(['Edit script "SetMatlabMultithreading.m" to choose work with version ( ' version('-release') ')'])
end
fprintf('Current Matlab pool is %g (0 is fine - no parallel loops)\n',CurMatPool);


%Display current number of active computational threads
fprintf('Current number of (active) computational threads is %g\n',maxNumCompThreads())


%Opens additional computation threads (=numProc_cluster) or opens
%matlabpool if on private computer
if ~isempty(numProc_cluster)
    fprintf('Number of available ppn on cluster is %g\n',numProc_cluster)
    
    %matlabpool('open',numProc_cluster);
    if numProc_cluster ~= 1
        LASTppn = maxNumCompThreads(numProc_cluster);
        fprintf(['     %g maxNumCompThreads opened [for implicit'...
            'parallelism (previously set to %g)].\n\n'],numProc_cluster,LASTppn)
        
    else
        fprintf(['     Cluster automatically initiates with -singleCompThread ' ...
            'option - this fits requested ppn; no action taken.\n\n'])
    end
    
    
else
    fprintf('No cluster registered - running on personal computer\n\n')
    if CurMatPool == 0 && ~( exist('nopool','var') && nopool == 1 )
        if any(strcmpi(version('-release'),{'2010a','2010b','2011a','2011b','2012a','2012b','2013a'}))
            matlabpool('open',numProc_private);
        elseif any(strcmpi(version('-release'),{'2013b','2014a','2014b','2015a','2015b','2016a','2016b','2017a','2017b','2018a','2018b','2019a','2019b'}))
            CurMatPool = parpool(numProc_private);
        else
            error(['Edit script "SetMatlabMultithreading.m" to choose work with version ( ' version('-release') ')'])
        end
    end
end

if exist('nopool','var'); clear nopool; end %clear this variable (which may be supplied to avoid opening a worker pool on a personal computer)