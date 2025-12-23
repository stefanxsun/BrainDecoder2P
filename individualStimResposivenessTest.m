function [StimVisRestest,StimVisRestestCat, centerResponsive, offcenterResponsive,StimVisResNum] ...
    = individualStimResposivenessTest(state, data, ResDayThreshold,p_threshold, centerResThreshold)

global pairedCellNum
global mouseTag

StimVisRestest = [];
StimVisRestestIdx =[];
centerResponsive ={};
ii=1;
for mousecnt = 1:length(mouseTag)
    for iday=1:size(data,2) 
        % find Index
        stimValue = data(mousecnt,iday).CaCellInfo.stimValue{state, 1};
        if size(stimValue,2) ==1 ; stimValue = [ones(size(stimValue))*0.04 stimValue];end
        sf = unique(stimValue(:,1));
        sze = unique(stimValue(:,2));
        for isf=1:length(sf)
            for isz=1:length(sze)  
                if size(stimValue,2)==2
                    sz = unique(stimValue(:,2));
                    ind = stimValue(:,1)==sf(isf) & stimValue(:,2) ==sz(isz);
                elseif size(stimValue,2)==3
                    sz = unique(stimValue(:,3));
                    ind = stimValue(:,1)==sf(isf) & stimValue(:,3) ==sz(isz)& stimValue(:,2) ==1;
                end
                for icell= 1:size(data(mousecnt,iday).CaCellMean.Ind,2)
                    res=data(mousecnt,iday).CaCellMean.Ind{state, icell}(ind);
                    pre=data(mousecnt,iday).CaCellMean.preInd{state, icell}(ind);
                    [h,~]=ttest(res,pre,'Tail','Right','Alpha',p_threshold);    %do a t test one tailed so that res must bigger than pre    
                    StimVisRestest{mousecnt,iday}{icell}(isz,isf)=h;
                end
            end
        end
        
    end
    
    
    for icell=1:pairedCellNum(mousecnt)
        StimVisRestestAllday  = cellfun(@(x) x{icell},StimVisRestest(mousecnt,:), 'UniformOutput', 0);
        StimVisRestestCount = sum(cat(3,StimVisRestestAllday{:}),3);
        StimVisRestestCountThresholded = StimVisRestestCount>=ResDayThreshold;
        StimVisResNum{mousecnt}(icell)=sum(StimVisRestestCountThresholded(:));
        StimVisRestestCountSizeSum = sum(StimVisRestestCountThresholded,2);
        StimVisRestestCat{mousecnt}(icell) =0;
        if     StimVisRestestCountSizeSum(1)>0
            StimVisRestestCat{mousecnt}(icell) =1;
        elseif StimVisRestestCountSizeSum(2)>0
            StimVisRestestCat{mousecnt}(icell) =2;
        elseif StimVisRestestCountSizeSum(3)>0
            StimVisRestestCat{mousecnt}(icell) =3;
        elseif StimVisRestestCountSizeSum(4)>0
            StimVisRestestCat{mousecnt}(icell) =4;
        else
            StimVisRestestCat{mousecnt}(icell) =5;  %% non responsive
        end
    end
    
    centerResponsive{mousecnt} =find(StimVisRestestCat{mousecnt}<=centerResThreshold);
    offcenterResponsive{mousecnt} =find(StimVisRestestCat{mousecnt}==3|StimVisRestestCat{mousecnt}==4);
end
end