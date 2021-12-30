%% Fucntion that inputs the EV onset time and duration for HCP data
%  and returns which time frames correspond to that

function [index_frames]= FunctionEVs(onset,duration,TR)

%TR=0.72; %in seconds
%duration=27.5; %in seconds
%onset=7.997; %time of event
end_onset=onset+duration;

aux1=TR:TR:301;
aux2=0:TR:301-TR;
ninterv=size(aux2,2);
Interval=zeros(ninterv,2);
Interval(:,1)=aux2';
Interval(:,2)=aux1';
auxInter=zeros(ninterv,2);
for i=1:ninterv
auxInter(i,1)= onset<=Interval(i,2) && onset>Interval(i,1);
auxInter(i,2)= end_onset<=Interval(i,2) && end_onset>Interval(i,1);
end

%Check if the volume was recorded at the onset times, i.e. after mid of the
%interval

index_onset1=find(auxInter(:,1)==1);
index_onset2=find(auxInter(:,2)==1);

%find midpoint of the interval
mid_onset1=Interval(index_onset1,2)-TR/2;
mid_onset2=Interval(index_onset2,2)-TR/2;

% If onset above the midpoint, then frame of the interval counts
Vol_indicator=onset>mid_onset1;
Vol_indicator_end=end_onset>mid_onset2;

% Which frames belong to the onset: onset+duration?
    
Frames=index_onset1+1-Vol_indicator:1:index_onset2-1+Vol_indicator_end;

index_frames=Frames';


