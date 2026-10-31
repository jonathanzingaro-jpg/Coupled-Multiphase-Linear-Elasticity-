function [osc_SGS,smear_SGS] = calcScore(hB_SGS, nodes, H)

numtimes = size(hB_SGS,2);

% % % smear_SGS=min((max(hB_SGS(:,numtimes+1))-min(hB_SGS(:,numtimes+1))-H)/H,0);
smear_SGS=abs(min(hB_SGS(:,numtimes+1))/abs(min(hB_SGS(:,numtimes+1)))*min((max(hB_SGS(:,numtimes+1))-abs(min(hB_SGS(:,numtimes+1)))-H)/H,0));

osc_SGS=0;
for i = 1:length(nodes)-1
    osc_SGS=osc_SGS+min(hB_SGS(i+1,numtimes+1)-hB_SGS(i,numtimes+1),0);
end
osc_SGS=abs(sum(osc_SGS)/H);
