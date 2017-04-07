function correcttr4(fname)
%fname example './TT3_20160812_133323.tr4'
v=load(fname,'-mat');
v2.whisker=v.whisker;
v2.trpmtrs=v.trpmtrs;
v2.roi=v.roi;
v2.calib=v.calib;
save(fname,'-struct','v2')
clear v v2
end