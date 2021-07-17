% inherit from "05d flatgrids.m"
% Note: we assume we already ran "05d flatgrids fsaverage only.m"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This document has steps that have to be routinely performed when doing
% anatomical pre-processing.

%%%%% For each subject, construct patches (based on the fsaverage-defined ROIs)
%%%%% and then flatten these patches.
%%%%% [THIS STEP NEEDS TO BE DONE IN ROUTINE ANATOMICAL PRE-PROCESSING.]

% The end result of this process will be things like:
%   surf/lh.gVTC.patch.3d
%   surf/lh.gVTC.flat.patch.3d
%   surf/lh.DENSE.gVTC.patch.3d
%   surf/lh.DENSE.gVTC.flat.patch.3d
%   surf/lh.DENSETRUNCpt.gVTC.patch.3d
%   surf/lh.DENSETRUNCpt.gVTC.flat.patch.3d
%
% Note that fsaverage is processed as if it were a subject.
% Note also that some subjects get only 'orig' processing (not 'DENSE', not 'DENSETRUNCpt').
%
% Note that the flattening process sometimes drops vertices (and these seem to be
% vertices on the outskirts/edges of the cut patch), so be careful!!
%
% Note that we now explicitly use "-dilate 2" so the above comment is even more true!!

% define

% dense processing:
subjectids = {'fsaverage'};   % WE DO NOT RE-RUN, BUT INCLUDE HERE FOR COMPLETENESS
subjectids = arrayfun(@(x) sprintf('subj%02d',x),1:8,'UniformOutput',0);

% define more
hemis = {'lh' 'rh'};
surfsuffixes = {'orig'};   % DROP: 'DENSETRUNCpt' 'DENSE'};

% % non-dense processing:
% subjectids = {'CVNS013' 'CVNS014'};
% subjectids = {'CVCS001','CVCS002','CVCS003','CVCS004','CVCS005','CVCS006','CVCS007','CVCS008','CVCS009','CVCS010'};
% hemis = {'lh' 'rh'};
% surfsuffixes = {'orig'};

% what to do
patches = {'gVTC' 'hVTC' 'gEVC'};

%%%%% Copy+paste this section directly into matlab (and not evalme - issue with parfor)
% do it
for zz=1:length(patches)
  parfor p=1:length(subjectids)
    for r=1:length(surfsuffixes)
      for q=1:length(hemis)
        subjectid = subjectids{p}
  
        % read in files
        surf0 = cvnreadsurface(subjectid,hemis{q},'inflated',surfsuffixes{r});
        roimask0 = cvnroimask(subjectid,hemis{q},patches{zz},[],surfsuffixes{r},'collapsebinary');

        % construct patch information
        patch = struct;
        patch.x = surf0.vertices(roimask0==1,1);
        patch.y = surf0.vertices(roimask0==1,2);
        patch.z = surf0.vertices(roimask0==1,3);
        patch.npts = sum(roimask0==1);
        patch.ind = find(roimask0==1);
  
        % determine faces that involve exactly 2 vertices that are all in our ROI (logical, FACES x 1)
        faceix = sum(ismember(surf0.faces,patch.ind),2)==2;
  
        % determine the vertices involved in the above-mentioned faces that are also part of the ROI
        patch.edge = intersect(union(flatten(surf0.faces(faceix,:)),[]),flatten(patch.ind));  % row vector, vertex indices
  
        % calc
        if isequal(surfsuffixes{r},'orig')
          suffix0 = '';
        else
          suffix0 = ['.' surfsuffixes{r}];
        end
      
        % some filenames
        file0 = sprintf('%s/%s/surf/%s%s.%s.patch.3d',         cvnpath('freesurfer'),subjectid,hemis{q},suffix0,patches{zz});
        file1 = sprintf('%s/%s/surf/%s%s.%s.flat.patch.3d',    cvnpath('freesurfer'),subjectid,hemis{q},suffix0,patches{zz});
          % some notes:
          % ??file2 = sprintf('%s/%s/surf/%s.gVTC.flat.patch.3d.asc',cvnpath('freesurfer'),subjectid,hemis{q});
          % file0 = sprintf('surf/%s.%s.flat.patch.3d.asc',prefix0,patchtypes{q});
          % surf = read_patch_asc(file0);
          %     %        faces: [298081x3 double]
          %     %     vertices: [149552x1 double]
          %     %            x: [149552x1 double]
          %     %            y: [149552x1 double]
          %     %            z: [149552x1 double]
          %     % faces are 0-indexes relative to the original faces
          %     % vertices are 0-indexes indicating which vertices in the original

        % save the patch to a file
        fast_write_patch_kj(file0,patch);
      
        % flatten the patch
        if isequal(surfsuffixes{r},'orig')

          % in this case, we can just do the command directly
          unix_wrapper(sprintf('mris_flatten -dilate 2 %s %s',file0,file1));   % -w 10
      
        else
      
          % in this case, we have to simulate the DENSE or DENSETRUNCpt surface
          file2 = sprintf('%s/%s/surf/%s.smoothwm%s',           cvnpath('freesurfer'),subjectid,hemis{q},surfsuffixes{r});
          dir0 = maketempdir;     % make a temp directory
          copyfile(file0,dir0);   % put the patch into that directory
          copyfile(file2,sprintf('%s/%s.smoothwm',dir0,hemis{q}));  % put the smoothwmDENSE[TRUNCpt] as smoothwm in that directory
          unix_wrapper(sprintf('mris_flatten -dilate 2 %s/%s %s',dir0,stripfile(file0,1),file1));   % -w 10
          rmdir(dir0,'s');

        end

      end
    end
  end
end

%%%%% For all subjects (including fsaverage!), automatically determine the best rotation
%%%%% to match the template values. Then actually perform this rotation and
%%%%% overwrite the patches.
%%%%% [THIS STEP NEEDS TO BE DONE IN ROUTINE ANATOMICAL PRE-PROCESSING.]

% define

% dense processing:
subjectids = {'fsaverage'};   % WE DO NOT RE-RUN, BUT INCLUDE HERE FOR COMPLETENESS
subjectids = arrayfun(@(x) sprintf('subj%02d',x),1:8,'UniformOutput',0);
hemis = {'lh' 'rh'};
surfsuffixes = {'orig'};   % DROPPED: 'DENSETRUNCpt' 'DENSE'};

% % non-dense processing:
% subjectids = {'CVNS013' 'CVNS014'};
% subjectids = {'CVCS001','CVCS002','CVCS003','CVCS004','CVCS005','CVCS006','CVCS007','CVCS008','CVCS009','CVCS010'};
% hemis = {'lh' 'rh'};
% surfsuffixes = {'orig'};

% which
patches = {'gVTC' 'hVTC' 'gEVC'};

% do it
for zz=1:length(patches)

  % load appropriate file
  load(sprintf('/home/stone/generic/Dropbox/cvnlabsmall/code/analysis/05d flatgrids/fsaverage rotations/%s/record.mat',patches{zz}),'avals','tangles');

  for p=1:length(subjectids)
    for r=1:length(surfsuffixes)
      for q=1:length(hemis)
        subjectid = subjectids{p}

        % calc
        if isequal(surfsuffixes{r},'orig')
          suffix0 = '';
          suffix1 = '';
        else
          suffix0 = ['.' surfsuffixes{r}];
          suffix1 = surfsuffixes{r};
        end
      
        % load stuff
        patchfile = sprintf('%s/%s/surf/%s%s.%s.flat.patch.3d',cvnpath('freesurfer'),subjectid,hemis{q},suffix0,patches{zz});
        a1 = fast_read_patch_kj(patchfile);
        vals = cvnroimask(subjectid,hemis{q},'aparc.a2009s',[],surfsuffixes{r},'collapsevals');
        %%curv = cvnreadsurfacemetric(subjectid,hemis{q},'curv','',surfsuffixes{r});
        roi = logical(copymatrix(zeros(size(vals)),double(a1.vno),1));

        % calc stuff
        valsA = vals(roi==1);  % aparc
        %%curvA = 2*((curv(roi==1) < 0) - .5);  % curvature

        % place COM at origin
        a1.x = zeromean(a1.x);
        a1.y = zeromean(a1.y);

        % calculate errors for all integer rotations
        err = [];
        for rot=1:360
      
          % copy a1 to a2 and rotate it    [MAYBE WE SHOULD USE procrustes.m!]
          a2 = a1;
          ang = rot/180*pi;
          rotmat = [cos(ang) -sin(ang) 0;
                    sin(ang)  cos(ang) 0;
                         0       0     1];
          temp = rotmat * [a2.x; a2.y; ones(1,length(a2.x))];
          a2.x = temp(1,:);
          a2.y = temp(2,:);

          % measure angles
          tangles0 = NaN*zeros(length(avals{q,r}));
          nums = [];
          for aaa=1:length(avals{q,r})
            for bbb=1:length(avals{q,r})
              if aaa>=bbb
                continue;
              end
              ix1 = find(valsA==avals{q,r}(aaa));
              ix2 = find(valsA==avals{q,r}(bbb));
              nums(aaa) = length(ix1);  % just the first for-loop
              if isempty(ix1) || isempty(ix2)
                % do nothing
              else
                tangles0(aaa,bbb) = atan2(mean(a2.y(ix2))-mean(a2.y(ix1)),mean(a2.x(ix2))-mean(a2.x(ix1)));
              end
            end
          end
          badix = nums < (median(nums) * (1/2));
          tangles0(badix,:) = NaN;  % small clusters are DELETED!
          tangles0(:,badix) = NaN;  % small clusters are DELETED!
        
          % summarize the error (L1)
          err(rot) = nansum(flatten(abs(circulardiff(tangles0,tangles{q,r},2*pi))));

        end
      
        % which rotation minimizes error?
        [dd,mnix] = min(err);
          %figure;plot(err);drawnow;
      
        % apply that rotation to a1
        ang = mnix/180*pi;
        rotmat = [cos(ang) -sin(ang) 0;
                  sin(ang)  cos(ang) 0;
                       0       0     1];
        temp = rotmat * [a1.x; a1.y; ones(1,length(a1.x))];
        a1.x = temp(1,:);
        a1.y = temp(2,:);
      
        % save the patch (NOTE: the only thing we change is .x and .y)
        fast_write_patch_kj(patchfile,a1);
        
        % calculate badness and write out these files
        flatsurf = cvnreadsurface(subjectid,hemis{q},sprintf('%s.flat.patch.3d',patches{zz}),surfsuffixes{r});  % load the flat patch
          % calculate face normals and determine which ones seem bad
        fN=facenormals(flatsurf.vertices,flatsurf.faces);
        fN=fN(:,3)>0;
        badfaces=fN~=round(mean(fN));
          % determine a badness matrix S (vertices x faces, filled with 1s were badness occurs)
        fidx=repmat((1:size(flatsurf.faces,1))',[1 3]);
        S=sparse(flatsurf.faces(:),fidx(:),badfaces(fidx(:)),...
            size(flatsurf.vertices,1),size(flatsurf.faces,1));
          % which vertices are involved in any bad face? (column vector of 0s/1s)
        bad = full(max(S,[],2))>0;
          % save this information (like, "label/rhDENSETRUNCpt.gVTC.patch.flat.3d.badness.mgz")
        cvnwritemgz(subjectid,sprintf('%s.flat.patch.3d.badness',patches{zz}),bad,[hemis{q} suffix1], ...
          [cvnpath('freesurfer') '/' subjectid '/label/']);

        % inspect the results
        figureprep([100 100 600 600]);
        scattersparse(a1.x,a1.y,5000,1,16,valsA,'o'); axis equal; colormap(jet); %caxis([-2 2]);
        xiqr = iqr(a1.x);
        yiqr = iqr(a1.y);
        ix0 = ismember(double(a1.ind),double(a1.edge));
        scatter(a1.x(ix0),a1.y(ix0),'k.');
        ix1 = bad(double(a1.vno));
        scatter(a1.x(ix1),a1.y(ix1),49,'kx');
        axis([-2*xiqr 2*xiqr -2*yiqr 2*yiqr]);
        figurewrite(sprintf('%s_%s_%s',subjectid,hemis{q},surfsuffixes{r}),[],[], ...
          sprintf('/home/stone/generic/Dropbox/cvnlabsmall/code/analysis/05d flatgrids/rotationresults/%s',patches{zz}));
 
      end
    end
  end
end

%%%%% Notes:

% - FS often creates outlier vertices (vertices that are placed far away from the bulk
%   of the vertices) in the flattening process, so beware of range issues.
%   These outlier vertices are probably just weird vertices located on the
%   boundary of the patch, and so aren't really very important. cvnlookupimages.m and
%   cvnreadsurface.m internally just identifies and ignores these outlier vertices.
%
% - FS's 'edge' field that it places in the flat patch files seems to leave the 
%   truncated edge unmarked. We don't really make use of the 'edge' field, 
%   so this weirdness doesn't matter.
%
% - FS sometimes creates loops and other topological defects in the flat patches.
%   These defects are relatively small in size but occur in about 25% of subjects/patches.
%   Our 'badness' calculation is a reasonable check as to where these defects are,
%   and the strategy will be to anatviz the badness values (which are saved in files like
%   "label/rhDENSETRUNCpt.gVTC.patch.flat.3d.badness.mgz"). Loops won't actually show
%   up after cvnreadsurface.m's clean-up operations (so any values associated with those
%   loops can't been seen), whereas other topological defects (squashed warts) may 
%   still show up (these show up in funny ways - it is just a nearest-neigbhor operation,
%   so you could get a weird mixture of wart and below-wart surface values).
%
% - The surface area of flat patch triangles appear to be systematically smaller than what
%   they are on *.white. This seems to be because mris_flatten works in reference to 
%   *.smoothwm, which is smaller in surface area (~10% less) compared to *.white.
%   The way we handle all of this is that cvnreadsurface.m, when reading in flat patches,
%   automatically scales the x- y-coordinates in the flat patch by a factor equal to
%   surface-area(white)/surface-area(flat) (calculated over the faces present in
%   the flat patch). After this scaling, we can assume that (on average across space)
%   there is a direct relationship between pixel distances and real physical units (mm)
%   measured on the *.white surfaces.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JUNK BELOW:

%       % what are the indices of these, relative to the patch.ind vertices?
%       [~,patch.edge] = ismember(edgeverts,patch.ind);

% a1 = fast_read_patch_kj('lh.gVTC.flat.patch.3d');
% figure; plot(a1.x,a1.y,'ro','MarkerSize',1); axis equal;
% a1 = fast_read_patch_kj('lh.DENSE.gVTC.flat.patch.3d');
% figure; plot(a1.x,a1.y,'ro','MarkerSize',1); axis equal;
% a1 = fast_read_patch_kj('lh.DENSETRUNCpt.gVTC.flat.patch.3d');
% figure; plot(a1.x,a1.y,'ro','MarkerSize',1); axis equal;
