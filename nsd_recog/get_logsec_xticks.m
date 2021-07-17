function [xtdvals, xtdlabels] = get_logsec_xticks


trial = 1;
minute = log(60);
hour = log(60*60);
day = log(24*60*60);
week = log(24*60*60*7);
month = log(24*60*60*30);
year = log(24*60*60*365);

xtdvals = [trial minute hour day week month year];
xtdlabels = {'1 Trial','1 Minute','1 Hour','1 Day','1 Week','1 Month','1 Year'};