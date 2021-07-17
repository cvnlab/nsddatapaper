% modeled against "analysis_freesurfer_flatgrids.m"

% NOTE:
% First, run for fsaverage.
% Then, run for the subjects.

% dense processing:
subjectids = {'fsaverage'};
subjectids = arrayfun(@(x) sprintf('subj%02d',x),1:8,'UniformOutput',0);
hemis = {'lh' 'rh'};
surfsuffixes = {'orig'};  % 'DENSETRUNCpt' 'DENSE'};

% which
patches = {'full'};

% do it
for zz=1:length(patches)

  for p=1:length(subjectids)
    subjectid = subjectids{p};

    % load appropriate file
    if isequal(subjectid,'fsaverage')
      tangles = {};
    else
      load('~/nsd/ppdata/flattenedsurfaces.mat','tangles');
    end

    for r=1:length(surfsuffixes)
      for q=1:length(hemis)

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
        vals = cvnroimask(subjectid,hemis{q},'Kastner2015*',[],surfsuffixes{r},'collapsevals');  % aparc.a2009s
        %%curv = cvnreadsurfacemetric(subjectid,hemis{q},'curv','',surfsuffixes{r});
        roi = logical(copymatrix(zeros(size(vals)),double(a1.vno),1));

        % calc stuff
        valsA = vals(roi==1);  % aparc
        nn = 25;
        assert(max(valsA)==nn);
        %%curvA = 2*((curv(roi==1) < 0) - .5);  % curvature

        % place COM at origin
        a1.x = zeromean(a1.x);
        a1.y = zeromean(a1.y);

        % calculate errors for all integer rotations
        err = [];
        if isequal(subjectid,'fsaverage')
          rotvals = 360;
        else
          rotvals = 1:360;
        end
        for rot=rotvals, rot
      
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
%          whtodo = 1:25;  %picksubset(1:nn,10);
          tangles0 = NaN*zeros(nn,nn);  %length(avals{q,r}));
          nums = [];
          for aaa=1:nn%10%74%length(avals{q,r})
            for bbb=1:nn%10%74%length(avals{q,r})
              if aaa>=bbb
                continue;
              end
              ix1 = find(valsA==aaa);%avals{q,r}(aaa));
              ix2 = find(valsA==bbb);%avals{q,r}(bbb));
%              nums(aaa) = length(ix1);  % just the first for-loop
              if isempty(ix1) || isempty(ix2)
                % do nothing
              else
                tangles0(aaa,bbb) = atan2(mean(a2.y(ix2))-mean(a2.y(ix1)),mean(a2.x(ix2))-mean(a2.x(ix1)));
              end
            end
          end
%           badix = nums < (median(nums) * (1/2));
%           tangles0(badix,:) = NaN;  % small clusters are DELETED!
%           tangles0(:,badix) = NaN;  % small clusters are DELETED!
          if isequal(subjectid,'fsaverage')
            tangles{q,r} = tangles0;
          end
        
          % summarize the error (L1)
          err(rot) = nansum(flatten(abs(circulardiff(tangles0,tangles{q,r},2*pi))));

        end
        
        if isequal(subjectid,'fsaverage')
          continue;
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
        copyfile(patchfile,[patchfile '.backup']);  % delete when done
        fast_write_patch_kj(patchfile,a1);
        
 
      end
    end
  end
end

% save
if isequal(subjectid,'fsaverage')
  save('~/nsd/ppdata/flattenedsurfaces.mat','tangles');
end
