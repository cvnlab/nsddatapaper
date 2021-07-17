function analyzebehavior_nsd_filecheck

nsdsetup;

% loop over subject
for p=1:length(allnsdids), allnsdids{p}

  % find all scan directories
  dirs0 = matchfiles(sprintf('%s/rawdata/*-%s-nsd??',nsddir,allnsdids{p}));
  
  % for each scan session
  for q=1:length(dirs0), dirs0{q}
    
    % get the session number from the directory name
    temp = regexp(dirs0{q},'.+-nsd(\d+)?','tokens');

    % check that the session number is the expected one
    assert(q==str2double(temp{1}{1}));
    fprintf('.');
    
    % find all .mat files
    files0 = matchfiles([dirs0{q} '/mat_files_from_scan/*nsd*subj*run*.mat']);
    assert(length(files0)==12);
    fprintf('$');
    
    % loop over .mat files
    for r=1:length(files0)
    
      % extract values from the filename
      temp = regexp(files0{r},'.+_nsd_subj(\d+)_run(\d+)_exp(\d+).mat','tokens');
      
      % check that the subject number is expected
      assert(str2double(temp{1}{1})==str2double(allnsdids{p}(4:end)));
      fprintf('!');
      
      % check that experiment number is the scan session number
      assert(str2double(temp{1}{3})==q);
      fprintf('!');
      
      % check that run number is expected
      assert(str2double(temp{1}{2})==r);
      fprintf('!');

    end
    
    fprintf('\n');

  end

  fprintf('\n\n');

end
