function parallelisation = initParallelisation(parallelisation)

if ~strcmpi(parallelisation, 'min') && ~strcmpi(parallelisation, 'max')
  parallelisation = str2double(parallelisation);
  if isnan(parallelisation)
    errmsg = 'Error: Number of parallel cores is NaN';
    msgbox(errmsg,'Error','Error');
    return
  end
  if parallelisation > 1
    scheduler = parcluster(parallel.defaultClusterProfile);
    if parallelisation > scheduler.NumWorkers
      errmsg = 'Error: not enough processor cores to support the required level of parallelisation. Please check your system settings.';
      msgbox(errmsg,'Error','Error');
      return
    end
  elseif parallelisation < 1
    errmsg = 'Error: Number of parallel cores is smaller than 1';
    msgbox(errmsg,'Error','Error');
    return
  end
end

if ischar(parallelisation) || parallelisation > 1
  parallelPool = gcp('nocreate');
  if isempty(parallelPool)
    openWorkers = 1;
  else
    openWorkers = parallelPool.NumWorkers;
  end
  scheduler = parcluster(parallel.defaultClusterProfile);
  switch parallelisation
    case {'min', 1}
      parallelisation = 1;
      if parallelisation ~= openWorkers
        delete(parallelPool)
      end
    case 'max'
      parallelisation = scheduler.NumWorkers;
      if parallelisation ~= openWorkers
        delete(parallelPool)
        parpool(parallelisation);
      end
    otherwise
      parallelisation = floor(parallelisation);
      if parallelisation ~= openWorkers
        if parallelisation < 1
          parallelisation = 1;
          delete(parallelPool)
        elseif parallelisation > scheduler.NumWorkers
          parallelisation = scheduler.NumWorkers;
          delete(parallelPool)
          parpool(parallelisation);
        else
          delete(parallelPool)
          parpool(parallelisation);
        end
      end
  end
  disp(['Using ' num2str(parallelisation) ' parallel cores.']);
end