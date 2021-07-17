for yy=1:length(viewswanted)
  Lookup = [];  % for each new view, need to reset the cache

  for zz=1:length(atlases)

    % define
    if ischar(atlases{zz}{1})
      dataname = atlases{zz}{1};
      if isequal(atlases{zz}{1},'curvature')
        data = [];  % in this case, we just allow default mechanisms to take hold
      else
        inputfile = sprintf([nsd_datalocation '/freesurfer/%s/label/??.%s.mgz'],dataid,atlases{zz}{1});
        data = cvnloadmgz(inputfile);
      end
    else
      dataname = atlases{zz}{1}{1};
      data = atlases{zz}{1}{2};
    end
    
    if isempty(atlases{zz}{5})
      outname = dataname;
    else
      outname = atlases{zz}{5};
    end
  
    outputfile = sprintf('%s/%s%s_%s_%s.png',outputdir,subjectid,viewnames{yy},outname,dataid);
    outputfile2 = sprintf('%s/%s%s_%s_names_%s.png',outputdir,subjectid,viewnames{yy},outname,dataid);
    outputfile3 = sprintf('%s/%s%s_%s_nothresh_%s.png',outputdir,subjectid,viewnames{yy},outname,dataid);

    % perform lookup and save to .png file
    [rawimg,Lookup,rgbimg,himg] = cvnlookup(subjectid,viewswanted(yy),data,atlases{zz}{2},atlases{zz}{3},atlases{zz}{4},Lookup,0,extraopts);
    imwrite(rgbimg,outputfile);
    
    if ischar(atlases{zz}{1}) && ~isequal(dataname,'curvature')

      % add roinames and save to another .png file
      [~,roinames,~] = cvnroimask(dataid,hemis,[dataname '*'],[],'orig','cell');
      roinames = regexprep(roinames{1},'@.+','');
      rgbimg = drawroinames(rawimg,rgbimg,Lookup,1:numel(roinames),cleantext(roinames));
      imwrite(rgbimg,outputfile2);

      % add no thresholding and no shading version
      [rawimg,Lookup,rgbimg,himg] = cvnlookup(subjectid,viewswanted(yy),data,atlases{zz}{2},atlases{zz}{3},[],Lookup,0,[extraopts {'surfshading' false}]);
      imwrite(rgbimg,outputfile3);

    end
  end

end
