% data directory names
dirs = { ...
'DATECENSORED-CVNS001-readingF001'
'DATECENSORED-CVNS001-readingF002'
'DATECENSORED-CVNS001-readingF003'
'DATECENSORED-CVNS001-readingF004'
'DATECENSORED-CVNS001-readingF005'
'DATECENSORED-CVNS001-readingF006'
'DATECENSORED-NSD139-prffloc'  % session 1
'DATECENSORED-NSD149-prffloc'
'DATECENSORED-NSD254-prffloc'
'DATECENSORED-NSD814-prffloc'
'DATECENSORED-NSD244-prffloc'
'DATECENSORED-NSD929-prffloc'
'DATECENSORED-NSD433-prffloc'
'DATECENSORED-NSD841-prffloc'
'DATECENSORED-NSD832-prffloc'
'DATECENSORED-NSD400-prffloc'
'DATECENSORED-NSD258-prffloc'
'DATECENSORED-NSD134-prffloc'
'DATECENSORED-NSD432-prffloc'
'DATECENSORED-NSD168-prffloc'
'DATECENSORED-NSD139-nsd01'
'DATECENSORED-NSD814-nsd01'
'DATECENSORED-NSD134-nsd01'
'DATECENSORED-NSD400-nsd01'
'DATECENSORED-NSD929-nsd01'
'DATECENSORED-NSD814-nsd02'
'DATECENSORED-NSD258-nsd01'
'DATECENSORED-NSD149-nsd01'
'DATECENSORED-NSD139-nsd02'
'DATECENSORED-NSD400-nsd02'
'DATECENSORED-NSD168-nsd01'
'DATECENSORED-NSD814-nsd03'
'DATECENSORED-NSD258-nsd02'
'DATECENSORED-NSD139-nsd03'
'DATECENSORED-NSD134-nsd02'
'DATECENSORED-NSD149-nsd02'
'DATECENSORED-NSD139-nsd04'
'DATECENSORED-NSD134-nsd03'
'DATECENSORED-NSD149-nsd03'
'DATECENSORED-NSD929-nsd02'
'DATECENSORED-NSD400-nsd03'
'DATECENSORED-NSD168-nsd02'
'DATECENSORED-NSD814-nsd04'
'DATECENSORED-NSD258-nsd03'
'DATECENSORED-NSD134-nsd04'
'DATECENSORED-NSD149-nsd04'
'DATECENSORED-NSD929-nsd03'
'DATECENSORED-NSD400-nsd04'
'DATECENSORED-NSD168-nsd03'
'DATECENSORED-NSD929-nsd04'
'DATECENSORED-NSD814-nsd05'
'DATECENSORED-NSD139-nsd05'
'DATECENSORED-NSD134-nsd05'
'DATECENSORED-NSD149-nsd05'
'DATECENSORED-NSD139-nsd06'
'DATECENSORED-NSD929-nsd05'
'DATECENSORED-NSD400-nsd05'
'DATECENSORED-NSD168-nsd04'
'DATECENSORED-NSD814-nsd06'
'DATECENSORED-NSD258-nsd04'
'DATECENSORED-NSD134-nsd06'
'DATECENSORED-NSD149-nsd06'
'DATECENSORED-NSD814-nsd07'
'DATECENSORED-NSD400-nsd06'
'DATECENSORED-NSD168-nsd05'
'DATECENSORED-NSD814-nsd08'
'DATECENSORED-NSD258-nsd05'
'DATECENSORED-NSD139-nsd07'
'DATECENSORED-NSD134-nsd07'
'DATECENSORED-NSD139-nsd08'
'DATECENSORED-NSD134-nsd08'
'DATECENSORED-NSD149-nsd07'
'DATECENSORED-NSD814-nsd09'
'DATECENSORED-NSD168-nsd06'
'DATECENSORED-NSD168-nsd07'
'DATECENSORED-NSD814-nsd10'
'DATECENSORED-NSD258-nsd06'
'DATECENSORED-NSD139-nsd09'
'DATECENSORED-NSD134-nsd09'
'DATECENSORED-NSD149-nsd08'
'DATECENSORED-NSD929-nsd06'
'DATECENSORED-NSD258-nsd07'
'DATECENSORED-NSD929-nsd07'
'DATECENSORED-NSD814-nsd11'
'DATECENSORED-NSD258-nsd08'
'DATECENSORED-NSD139-nsd10'
'DATECENSORED-NSD134-nsd10'
'DATECENSORED-NSD149-nsd09'
'DATECENSORED-NSD168-nsd08'
'DATECENSORED-NSD400-nsd07'
'DATECENSORED-NSD168-nsd09'
'DATECENSORED-NSD814-nsd12'
'DATECENSORED-NSD258-nsd09'
'DATECENSORED-NSD139-nsd11'
'DATECENSORED-NSD134-nsd11'
'DATECENSORED-NSD149-nsd10'
'DATECENSORED-NSD929-nsd08'
'DATECENSORED-NSD168-nsd10'
'DATECENSORED-NSD929-nsd09'
'DATECENSORED-NSD400-nsd08'
'DATECENSORED-NSD168-nsd11'
'DATECENSORED-NSD929-nsd10'
'DATECENSORED-NSD258-nsd10'
'DATECENSORED-NSD139-nsd12'
'DATECENSORED-NSD134-nsd12'
'DATECENSORED-NSD149-nsd11'
'DATECENSORED-NSD168-nsd12'
'DATECENSORED-NSD400-nsd09'
'DATECENSORED-NSD168-nsd13'
'DATECENSORED-NSD929-nsd11'
'DATECENSORED-NSD400-nsd10'
'DATECENSORED-NSD258-nsd11'
'DATECENSORED-NSD139-nsd13'
'DATECENSORED-NSD134-nsd13'
'DATECENSORED-NSD149-nsd12'
'DATECENSORED-NSD168-nsd14'
'DATECENSORED-NSD400-nsd11'
'DATECENSORED-NSD168-nsd15'
'DATECENSORED-NSD814-nsd13'
'DATECENSORED-NSD258-nsd12'
'DATECENSORED-NSD139-nsd14'
'DATECENSORED-NSD134-nsd14'
'DATECENSORED-NSD149-nsd13'
'DATECENSORED-NSD168-nsd16'
'DATECENSORED-NSD814-nsd14'
'DATECENSORED-NSD400-nsd12'
'DATECENSORED-NSD168-nsd17'
'DATECENSORED-NSD814-nsd15'
'DATECENSORED-NSD258-nsd13'
'DATECENSORED-NSD139-nsd15'
'DATECENSORED-NSD134-nsd15'
'DATECENSORED-NSD149-nsd14'
'DATECENSORED-NSD929-nsd12'
'DATECENSORED-NSD929-nsd13'
'DATECENSORED-NSD400-nsd13'
'DATECENSORED-NSD258-nsd14'
'DATECENSORED-NSD929-nsd14'
'DATECENSORED-NSD814-nsd16'
'DATECENSORED-NSD134-nsd16'
'DATECENSORED-NSD139-nsd16'
'DATECENSORED-NSD149-nsd15'
'DATECENSORED-NSD929-nsd15'
'DATECENSORED-NSD400-nsd14'
'DATECENSORED-NSD168-nsd18'
'DATECENSORED-NSD814-nsd17'
'DATECENSORED-NSD139-nsd17'
'DATECENSORED-NSD258-nsd15'
'DATECENSORED-NSD149-nsd16'
'DATECENSORED-NSD168-nsd19'
'DATECENSORED-NSD929-nsd16'
'DATECENSORED-NSD400-nsd15'
'DATECENSORED-NSD814-nsd18'
'DATECENSORED-NSD139-nsd18'
'DATECENSORED-NSD134-nsd17'
'DATECENSORED-NSD149-nsd17'
'DATECENSORED-NSD258-nsd16'
'DATECENSORED-NSD814-nsd19'
'DATECENSORED-NSD258-nsd17'
'DATECENSORED-NSD134-nsd18'
'DATECENSORED-NSD149-nsd18'
'DATECENSORED-NSD258-nsd18'
'DATECENSORED-NSD139-nsd19'
'DATECENSORED-NSD134-nsd19'
'DATECENSORED-NSD149-nsd19'
'DATECENSORED-NSD258-nsd19'
'DATECENSORED-NSD149-nsd20'
'DATECENSORED-NSD814-nsd20'
'DATECENSORED-NSD258-nsd20'
'DATECENSORED-NSD134-nsd20'
'DATECENSORED-NSD134-nsd21'
'DATECENSORED-NSD400-nsd16'
'DATECENSORED-NSD258-nsd21'
'DATECENSORED-NSD149-nsd21'
'DATECENSORED-NSD814-nsd21'
'DATECENSORED-NSD400-nsd17'
'DATECENSORED-NSD134-nsd22'
'DATECENSORED-NSD134-nsd23'
'DATECENSORED-NSD400-nsd18'
'DATECENSORED-NSD258-nsd22'
'DATECENSORED-NSD814-nsd22'
'DATECENSORED-NSD139-nsd20'
'DATECENSORED-NSD134-nsd24'
'DATECENSORED-NSD149-nsd22'
'DATECENSORED-NSD400-nsd19'
'DATECENSORED-NSD139-nsd21'
'DATECENSORED-NSD400-nsd20'
'DATECENSORED-NSD258-nsd23'
'DATECENSORED-NSD149-nsd23'
'DATECENSORED-NSD139-nsd22'
'DATECENSORED-NSD139-nsd23'
'DATECENSORED-NSD134-nsd25'
'DATECENSORED-NSD258-nsd24'
'DATECENSORED-NSD814-nsd23'
'DATECENSORED-NSD400-nsd21'
'DATECENSORED-NSD258-nsd25'
'DATECENSORED-NSD139-nsd24'
'DATECENSORED-NSD814-nsd24'   % END OF RELEASE
'DATECENSORED-NSD400-nsd22'
'DATECENSORED-NSD134-nsd26'
'DATECENSORED-NSD400-nsd23'
'DATECENSORED-NSD258-nsd26'
'DATECENSORED-NSD149-nsd24'
'DATECENSORED-NSD814-nsd25'
'DATECENSORED-NSD139-nsd25'
'DATECENSORED-NSD134-nsd27'
'DATECENSORED-NSD149-nsd25'
'DATECENSORED-NSD258-nsd27'
'DATECENSORED-NSD149-nsd26'
'DATECENSORED-NSD814-nsd26'
'DATECENSORED-NSD139-nsd26'
'DATECENSORED-NSD134-nsd28'
'DATECENSORED-NSD134-nsd29'
'DATECENSORED-NSD400-nsd24'
'DATECENSORED-NSD258-nsd28'
'DATECENSORED-NSD814-nsd27'
'DATECENSORED-NSD400-nsd25'
'DATECENSORED-NSD139-nsd27'
'DATECENSORED-NSD134-nsd30'
'DATECENSORED-NSD400-nsd26'
'DATECENSORED-NSD258-nsd29'
'DATECENSORED-NSD149-nsd27'
'DATECENSORED-NSD814-nsd28'
'DATECENSORED-NSD139-nsd28'
'DATECENSORED-NSD134-nsd31'   % END OF RELEASE
'DATECENSORED-NSD814-nsd29'
'DATECENSORED-NSD149-nsd28'
'DATECENSORED-NSD258-nsd30'
'DATECENSORED-NSD814-nsd30'
'DATECENSORED-NSD139-nsd29'
'DATECENSORED-NSD134-nsd32'
'DATECENSORED-NSD149-nsd29'
'DATECENSORED-NSD258-nsd31'
'DATECENSORED-NSD814-nsd31'
'DATECENSORED-NSD929-nsd17'
'DATECENSORED-NSD139-nsd30'
'DATECENSORED-NSD134-nsd33'
'DATECENSORED-NSD134-nsd34'
'DATECENSORED-NSD929-nsd18'
'DATECENSORED-NSD258-nsd32'   % END OF RELEASE
'DATECENSORED-NSD149-nsd30'   %.catchup
'DATECENSORED-NSD929-nsd19'   %
'DATECENSORED-NSD139-nsd31'   %
'DATECENSORED-NSD258-nsd33'   %.ENDcatchup
'DATECENSORED-NSD134-nsd35'
'DATECENSORED-NSD134-nsd36'
'DATECENSORED-NSD814-nsd32'
'DATECENSORED-NSD258-nsd34'
'DATECENSORED-NSD929-nsd20'
'DATECENSORED-NSD929-nsd21'
'DATECENSORED-NSD139-nsd32'
'DATECENSORED-NSD134-nsd37'
'DATECENSORED-NSD400-nsd27'
'DATECENSORED-NSD814-nsd33'
'DATECENSORED-NSD258-nsd35'
'DATECENSORED-NSD929-nsd22'
'DATECENSORED-NSD258-nsd36'
'DATECENSORED-NSD134-nsd38'
'DATECENSORED-NSD814-nsd34'
'DATECENSORED-NSD139-nsd33'
'DATECENSORED-NSD929-nsd23'
'DATECENSORED-NSD168-nsd20'
'DATECENSORED-NSD929-nsd24'
'DATECENSORED-NSD149-nsd31'
'DATECENSORED-NSD400-nsd28'
'DATECENSORED-NSD134-nsd39'
'DATECENSORED-NSD168-nsd21'
'DATECENSORED-NSD134-nsd40'
'DATECENSORED-NSD814-nsd35'
'DATECENSORED-NSD149-nsd32'
'DATECENSORED-NSD929-nsd25'
'DATECENSORED-NSD139-nsd34'
'DATECENSORED-NSD168-nsd22'
'DATECENSORED-NSD139-nsd35'
'DATECENSORED-NSD258-nsd37'
'DATECENSORED-NSD168-nsd23'
'DATECENSORED-NSD400-nsd29'
'DATECENSORED-NSD814-nsd36'
'DATECENSORED-NSD929-nsd26'
'DATECENSORED-NSD139-nsd36'
'DATECENSORED-NSD168-nsd24'
'DATECENSORED-NSD929-nsd27'
'DATECENSORED-NSD139-nsd37'
'DATECENSORED-NSD400-nsd30'
'DATECENSORED-NSD168-nsd25'
'DATECENSORED-NSD168-nsd26'
'DATECENSORED-NSD929-nsd28'
'DATECENSORED-NSD168-nsd27'
'DATECENSORED-NSD929-nsd29'
'DATECENSORED-NSD139-nsd38'
'DATECENSORED-NSD258-nsd38'
'DATECENSORED-NSD258-nsd39'
'DATECENSORED-NSD814-nsd37'
'DATECENSORED-NSD139-nsd39'
'DATECENSORED-NSD258-nsd40'
'DATECENSORED-NSD400-nsd31'
'DATECENSORED-NSD400-nsd32'
'DATECENSORED-NSD814-nsd38'
'DATECENSORED-NSD929-nsd30'
'DATECENSORED-NSD139-nsd40'
'DATECENSORED-NSD168-nsd28'
'DATECENSORED-NSD814-nsd39'
'DATECENSORED-NSD168-nsd29'
'DATECENSORED-NSD134-nsdsynthetic'
'DATECENSORED-NSD258-nsdsynthetic'
'DATECENSORED-NSD149-nsdsynthetic'
'DATECENSORED-NSD929-nsdsynthetic'
'DATECENSORED-NSD168-nsd30'
'DATECENSORED-NSD814-nsd40'
'DATECENSORED-NSD258-nsdimagery'
'DATECENSORED-NSD134-nsdimagery'
'DATECENSORED-NSD400-nsdsynthetic'
'DATECENSORED-NSD139-nsdsynthetic'
'DATECENSORED-NSD400-nsdimagery'
'DATECENSORED-NSD139-nsdimagery'
'DATECENSORED-NSD149-nsdimagery'
'DATECENSORED-NSD814-nsdsynthetic'
'DATECENSORED-NSD168-nsdsynthetic'
'DATECENSORED-NSD168-nsdimagery'
'DATECENSORED-NSD814-nsdimagery'
'DATECENSORED-NSD929-nsdimagery'
};
% 'DATECENSORED-NSD134-nsdsyntheticpilot'
% 'DATECENSORED-NSD134-nsdimageryFAILED'

% preprocessing scripts
projectnames = {'readingF' 'nsd'};
projectnums = 2*ones(1,length(dirs));
projectnums(1:6) = 1;

% calc
datadirs = {};
ppdirs = {};
glmdirs = {};
for p=1:length(dirs)
  datadirs{p} = sprintf('/home/surly-raid1/kendrick-data/%s/rawdata/%s/',projectnames{projectnums(p)},dirs{p});
  ppdirs{p} =   sprintf('/home/surly-raid1/kendrick-data/%s/ppdata/%s/', projectnames{projectnums(p)},dirs{p});
  glmdirs{p} =  sprintf('/home/surly-raid1/kendrick-data/%s/glmdata/%s/', projectnames{projectnums(p)},dirs{p});
end

% special
nsddir =          '/home/surly-raid1/kendrick-data/nsd/';
nsddatadir =      '/home/surly-raid4/kendrick-data/nsd/nsddata/';
nsddatabetasdir = '/home/surly-raid4/kendrick-data/nsd/nsddata_betas/';
nsddatatimeseriesdir = '/home/surly-raid3/kendrick-data/nsd/nsddata_timeseries/';
nsdfsdir =        '/home/surly-raid4/kendrick-data/nsd/nsddata/freesurfer/';
  %nsdfsdir =        '/home/stone-ext1/freesurfer/subjects/';

% Freesurfer IDs and dates
fsids = {};
scandates = [];
for p=1:length(dirs)
  temp = regexp(dirs{p},'(?<date>.+?)-(?<id>.+?)-.+','names');
  fsids{p} = temp.id;
  scandates(p) = str2double(temp.date);
end
fsids(1:6) = repmat({'C1051'},[1 6]);

% all NSD IDs
allnsdidsNUM = [134 139 149 168 258 400 814 929];
allnsdids = arrayfun(@(x) sprintf('NSD%d',x),allnsdidsNUM,'UniformOutput',0);

% structural sessions directory names
structuraldirs = {};
for p=1:length(allnsdids)
  structuraldirs{p} = matchfiles(sprintf('/home/surly-raid1/kendrick-data/nsd/rawdata/*-%s-structural*',allnsdids{p}));
end

% crops for NSD subjects
nsdcropranges = {};
nsdcropranges{1} = {[9 112] [20 100] [2 84]};
nsdcropranges{2} = {[9 114] [21 102] [1 84]};
nsdcropranges{3} = {[8 113] [22 102] [1 82]};
nsdcropranges{4} = {[12 110] [20 104] [4 83]};
nsdcropranges{5} = {[15 111] [21 99] [4 81]};
nsdcropranges{6} = {[6 118] [19 103] [2 84]};
nsdcropranges{7} = {[12 106] [23 100] [2 82]};
nsdcropranges{8} = {[11 113] [20 99] [4 81]};  %OLD nsdcropranges{8} = {[3 114] [15 104] [1 82]};

% the matrix size for the prepared EPI data.
% the two preps share the same "first" corner voxel (Anterior, Right, Inferior)
nsdmatrixsize = {};
nsdmatrixsizeLOW = {};
for p=1:length(nsdcropranges)
  crop0 = nsdcropranges{p};
  [xx,yy,zz] = ndgrid(crop0{1}(1):1/1.8:crop0{1}(2), ...
                      crop0{2}(1):1/1.8:crop0{2}(2), ...
                      crop0{3}(1):1/1.8:crop0{3}(2));
  nsdmatrixsize{p} = size(xx);
  [xx,yy,zz] = ndgrid(crop0{1}(1):1:crop0{1}(2), ...
                      crop0{2}(1):1:crop0{2}(2), ...
                      crop0{3}(1):1:crop0{3}(2));
  nsdmatrixsizeLOW{p} = size(xx);
end

% which session index is the reference/synthetic/imagery for each subject?
nsdrefsession = [23 21 28 31 27 24 22 25];
nsdsyntheticsession = [303 312 305 317 304 311 316 306];
nsdimagerysession = [310 314 315 318 309 313 319 320];

% final nsd sessions
totalnsdsessions = [40 40 32 30 40 32 40 30];

% make session matrix (8 x 40 with index of the session dir)
sessmatrix = [];
for p=1:length(dirs)
  temp = regexp(dirs{p},'.+-NSD(\d+)-nsd(\d+)','tokens');
  if ~isempty(temp)
    sessmatrix(find(allnsdidsNUM==str2double(temp{1}{1})), ...
               str2double(temp{1}{2})) = p;
  end
end

% make screening session matrix
screeningmatrix = [18 7 8 20 17 16 10 12];

% more
nsdsubjixfun = @(sessix) find(ismember(allnsdids,fsids{sessix}));
nsdidfun = @(sessix) allnsdids{nsdsubjixfun(sessix)};
nsdsessnumfun = @(sessix) str2double(getfield(regexp(dirs{sessix},'.+-nsd(?<id>\d+?)$','names'),'id'));
nsdscreeningfun = @(sessix) ~isempty(regexp(dirs{sessix},'.+-prffloc$'));
nsdscreeningfun2 = @(subjix) datadirs{screeningmatrix(subjix)};
nsdallsessfun = @(sessix) datadirs(filterout(sessmatrix(nsdsubjixfun(sessix),:),0));
nsdallsessfun2 = @(subjix) datadirs(filterout(sessmatrix(subjix,:),0));
nsdallsessfun3 = @(subjix) filterout(sessmatrix(subjix,:),0);

% DEPRECATED:
% % calc
% nsdsyntheticsessions = [];
% nsdimagerysessions = [];
% for p=1:length(dirs)
%   temp = regexp(dirs{p},'.+-NSD(\d+)-nsdsynthetic$','tokens');
%   if ~isempty(temp)
%     nsdsyntheticsessions = [nsdsyntheticsessions p];
%   end
%   temp = regexp(dirs{p},'.+-NSD(\d+)-nsdimagery$','tokens');
%   if ~isempty(temp)
%     nsdimagerysessions = [nsdimagerysessions p];
%   end
% end

% experiment parameters
expecteddim_readingF =    {{[0.8  0.8  0.8] 2.2   [162 200 84] repmat([211],[1 10])} ...
                           {[2    2    2.4] 0.391 [ 72  96 28] ones(1,4)}};
expecteddim_readingFalt = {{[0.8  0.8  0.8] 2.2   [162 200 84] repmat([160],[1 12])} ...
                           {[2    2    2.4] 0.391 [ 72  96 28] ones(1,4)}};
expecteddim_prffloc =     {{[1.8  1.8  1.8] 1.6   [120 120 84] repmat([188 188 195 195],[1 3])} ...
                           {[2.16 2.16 3.6] 0.51  [100 100 42] ones(1,4)}};
expecteddim_nsd =         {{[1.8  1.8  1.8] 1.6   [120 120 84] repmat(188,[1 12])} ...
                           {[2.16 2.16 3.6] 0.51  [100 100 42] ones(1,4)}};
expecteddim_nsd2 =        {{[1.8  1.8  1.8] 1.6   [120 120 84] repmat(188,[1 12])} ...
                           {[2.16 2.16 3.6] 0.51  [100 100 42] ones(1,5)}};
expecteddim_nsd3 =        {{[1.8  1.8  1.8] 1.6   [120 120 84] [188*ones(1,3) 111 70 188*ones(1,8)]} ...
                           {[2.16 2.16 3.6] 0.51  [100 100 42] ones(1,5)}};
expecteddim_nsd4 =        {{[1.8  1.8  1.8] 1.6   [120 120 84] repmat(188,[1 12])} ...
                           {[2.16 2.16 3.6] 0.51  [100 100 42] ones(1,3)}};
expecteddim_nsd5 =        {{[1.8  1.8  1.8] 1.6   [120 120 84] [repmat(181,[1 2]) repmat(188,[1 10])]} ...
                           {[2.16 2.16 3.6] 0.51  [100 100 42] ones(1,5)}};
expecteddim_nsdsynthetic ={{[1.8  1.8  1.8] 1.6   [120 120 84] [repmat(268,[1 8])]} ...
                           {[2.16 2.16 3.6] 0.51  [100 100 42] ones(1,3)}};
expecteddim_nsdimagery =  {{[1.8  1.8  1.8] 1.6   [120 120 84] [150 300 150 150 300 150 150 300 150 150 150 150]} ...
                           {[2.16 2.16 3.6] 0.51  [100 100 42] ones(1,5)}};
expecteddim_nsdrs =         {{[1.8  1.8  1.8] 1.6   [120 120 84] repmat(188,[1 14])} ...
                             {[2.16 2.16 3.6] 0.51  [100 100 42] ones(1,4)}};
expecteddim = {expecteddim_readingF expecteddim_readingFalt expecteddim_prffloc expecteddim_nsd expecteddim_nsd2 expecteddim_nsd3 expecteddim_nsd4 expecteddim_nsdrs expecteddim_nsd5 expecteddim_nsdsynthetic expecteddim_nsdimagery};
sessiontypes = 4*ones(1,length(datadirs));
sessiontypes(1:6) = [1 1 1 2 2 2];
sessiontypes(6+(1:14)) = 3;
sessiontypes(35) = 5;
sessiontypes(40) = 6;
sessiontypes(64) = 5;
sessiontypes(143) = 5;
sessiontypes(152) = 7;
sessiontypes(234) = 5;
sessiontypes(295) = 9;
sessiontypes(nsdsyntheticsession) = 10;
sessiontypes(nsdimagerysession) = 11;
for p=filterout(union(sessmatrix(:),[]),0)
  if nsdsessnumfun(p) >= 21 && nsdsessnumfun(p) <= 30 || ...    % THIS DEFINES WHEN RESTING-STATE HAPPENED
     (ismember(nsdsubjixfun(p),[1 5]) && ismember(nsdsessnumfun(p),31:38))
    sessiontypes(p) = 8;
  end
end

% which runs have behavioral data missing?
sessionbehaviormissing = zeros(8,40,12);
sessionbehaviormissing(6,19,10) = 1;
sessionbehaviormissing(4,23,3) = 1;

% preprocessing scripts
ppscripts = {'preprocess_readingF' 'preprocess_nsd'};
pptypes = 2*ones(1,length(datadirs));
pptypes(1:6) = 1;

% tweaks to preprocessing
fieldmapslicerangeALTs = cell(1,length(datadirs));
fieldmapslicerangeALTs{29} = {[2:100 100] ':' [ones(1,3) 1:84-3]};
fieldmaptweaks = cell(1,length(datadirs));
fieldmaptweaks{35} =  {[1 2 5 8 11] {1 1 1 [2 3] [3 4] [4 5] [5 6] [6 7] [7 8] [8 9] [9 10] [10 11]}};
fieldmaptweaks{40} =  {[1 5 5.5 10 14] {[1 2] [2 3] [3 4] [4 5] [5.5 6] [6 7] [7 8] [8 9] [9 10] [10 11] [11 12] [12 13] [13 14]}};
fieldmaptweaks{62} =  {[1 5 9 13] {2 2 2 2 [5 6] [6 7] [7 8] [8 9] [9 10] [10 11] [11 12] [12 13]}};
fieldmaptweaks{64} =  {[1 5 9 10 12] {[1 2] [2 3] [3 4] [4 5] [5 6] [6 7] [7 8] [8 9] 3 3 [10 11] [11 12]}};
fieldmaptweaks{101} = {[1 5 9 13] {[1 2] [2 3] [3 4] [4 5] [5 6] [6 7] [7 8] [8 9] 3 3 3 3}};
fieldmaptweaks{143} = {[1 5 9 10 14] {[1 2] [2 3] [3 4] [4 5] [5 6] [6 7] [7 8] [8 9]  [10 11] [11 12] [12 13] [13 14]}};
fieldmaptweaks{152} = {[1 5 9] {[1 2] [2 3] [3 4] [4 5] [5 6] [6 7] [7 8] [8 9] 3 3 3 3}};
fieldmaptweaks{234} = {[1 2 4 8 12] {1 1 [2 3] [3 4] [4 5] [5 6] [6 7] [7 8] [8 9] [9 10] [10 11] [11 12]}};
fieldmaptweaks{295} = {[1 2 4 8 12] {1 1 [2 3] [3 4] [4 5] [5 6] [6 7] [7 8] [8 9] [9 10] [10 11] [11 12]}};
for p=filterout(union(sessmatrix(:),[]),0)
  if sessiontypes(p)==8  % if resting-state session
    fieldmaptweaks{p} = {[1 6 10 15] splitmatrix([1:14; 2:15]',1)};
  end
end
for p=1:length(nsdsyntheticsession)
  fieldmaptweaks{nsdsyntheticsession(p)} = {[1 5 9] splitmatrix([1:8; 2:9]',1)};
end
for p=1:length(nsdimagerysession)
  fieldmaptweaks{nsdimagerysession(p)} = {[1 4 7 10 13] splitmatrix([1:12; 2:13]',1)};
end
fixepifuns = cell(1,length(datadirs));
fixepifuns{35} =  @(epis) preprocess_nsd_fixepi(epis,35);
fixepifuns{40} =  @(epis) preprocess_nsd_fixepi(epis,40);
fixepifuns{295} = @(epis) preprocess_nsd_fixepi(epis,295);

% more tweaks.
epitrtweaks = repmat({0.999878},[1 length(datadirs)]);  % for main data and synthetic (Sep 2018 update for new Stim Mac direct to BOLDscreen setup)
epitrtweaks(nsdimagerysession) = repmat({1.000060},[1 length(nsdimagerysession)]);  % for imagery
epitrtweaks{309} = repmat(1.000060,[1 12]);
epitrtweaks{309}(2) = 1.000821;
epitrtweaks{309}(7) = 1.000604;

% experiment information files
nsdexpfile = '/home/surly-raid4/kendrick-data/nsd/nsddata/experiments/nsd/nsd_expdesign.mat';
nsdsyntheticexpfile = '/home/surly-raid4/kendrick-data/nsd/nsddata/experiments/nsdsynthetic/nsdsynthetic_expdesign.mat';
nsdimageryexpfile = cellfun(@(x) sprintf('/home/surly-raid4/kendrick-data/nsd/nsddata/experiments/nsdimagery/%s_dm.mat',x),{'visA' 'attA' 'imgA_1' 'visB' 'attB' 'imgB_1' 'visC' 'attC' 'imgC_1' 'imgA_2' 'imgB_2' 'imgC_2'},'UniformOutput',0);

% any structurals to reject?
rejectt1 = {[] [] [] [] [] [] [] [2 4]};
rejectt2 = {[] [] [] [] [] [] [] []};
