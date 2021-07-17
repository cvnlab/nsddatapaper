% inherited from "05b fsaverage transfer mapping.m"

% This script creates and saves lookups to transfer to/from subject surfaces and fsaverage (both orig and DENSE)

%% Part 1: create subject fsaverage transfers

% define
subjects={'fsaverage'};  % DID NOT RE-RUN, BUT INCLUDE HERE FOR COMPLETENESS!
subjects=arrayfun(@(x) sprintf('subj%02d',x),1:8,'UniformOutput',0);
hemis={'lh','rh'};
surfsuffixes={'orig','DENSE'};

% do it
for s = 1:numel(subjects)
    subject=subjects{s};
    
    for ss = 1:numel(surfsuffixes)
        surfsuffix=surfsuffixes{ss};
        
        
        if(isequal(surfsuffix,'orig'))
            fsavgsuffix='';
        else
            fsavgsuffix=surfsuffix;
        end
        
        % e.g., if lh.whiteDENSE does not exist (e.g. because we did not do dense processing), just skip this case
        if ~exist(sprintf('%s/%s/surf/lh.white%s',cvnpath('freesurfer'),subject,fsavgsuffix),'file')
          continue;
        end
        
        
        for h = 1:numel(hemis)
            hemi=hemis{h};
            
            numvals=cvnreadsurface(subject,hemi,'sphere',surfsuffix,'justcount',true);
            
            %forward transform
            validix=cvntransfertosubject(subject,'fsaverage',(1:numvals)',hemi,'nearest',surfsuffix,surfsuffix);
            xferfile=sprintf('%s/%s/surf/%s.%s_to_fsaverage%s.mat',cvnpath('freesurfer'),subject,hemi,surfsuffix,fsavgsuffix);
            save(xferfile,'validix');
            
            %backward transform
            numvals_fsvals=numel(validix);
            validix=cvntransfertosubject('fsaverage',subject,(1:numvals_fsvals)',hemi,'nearest',surfsuffix,surfsuffix);
            xferfile=sprintf('%s/%s/surf/%s.fsaverage%s_to_%s.mat',cvnpath('freesurfer'),subject,hemi,fsavgsuffix,surfsuffix);
            save(xferfile,'validix');
            
        end
    end
end




%% NOTE: Part 2 is omitted here, but was indeed run on June 4 2019; see "05b fsaverage transfer mapping.m"
