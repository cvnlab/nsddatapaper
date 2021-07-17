function fa = get_sess_fa(sessix,beh)
% fa = get_sess_fa(sessix,beh)
%
% calculates session-specific false alarm rate given session number and beh
% struct taken from e.g. get_adj_hr.m
%
% jbh 2/22/21

sessInds = beh.SESSION == sessix;
newInds = beh.ISOLD == 0;
errorVals = double(beh.ISCORRECT==0);
fa = mean(errorVals(and(sessInds,newInds)));
