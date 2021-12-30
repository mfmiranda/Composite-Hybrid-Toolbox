%Generate Table 2
Table2=zeros(12,5);
load('ResultsSim100_CHSB_LongTerm')
Table2(1:4,1)=[mean(IntervalWidthArray);mean(MSEArray);mean(RPArray);mean(FPArray)];
Table2(1:4,2)=[std(IntervalWidthArray);std(MSEArray);std(RPArray);std(FPArray)];
Table2(1:4,3)=[min(IntervalWidthArray);min(MSEArray);min(RPArray);min(FPArray)];
Table2(1:4,4)=[median(IntervalWidthArray);median(MSEArray);median(RPArray);median(FPArray)];
Table2(1:4,5)=[max(IntervalWidthArray);max(MSEArray);max(RPArray);max(FPArray)];
clearvars -except Table2
load('ResultsSim100_LSB_LongTerm')
Table2(5:8,1)=[mean(IntervalWidthArray);mean(MSEArray);mean(RPArray);mean(FPArray)];
Table2(5:8,2)=[std(IntervalWidthArray);std(MSEArray);std(RPArray);std(FPArray)];
Table2(5:8,3)=[min(IntervalWidthArray);min(MSEArray);min(RPArray);min(FPArray)];
Table2(5:8,4)=[median(IntervalWidthArray);median(MSEArray);median(RPArray);median(FPArray)];
Table2(5:8,5)=[max(IntervalWidthArray);max(MSEArray);max(RPArray);max(FPArray)];
clearvars -except Table2
load('ResultsSim100_GSB_LongTerm')
Table2(9:12,1)=[mean(IntervalWidthArray);mean(MSEArray);mean(RPArray);mean(FPArray)];
Table2(9:12,2)=[std(IntervalWidthArray);std(MSEArray);std(RPArray);std(FPArray)];
Table2(9:12,3)=[min(IntervalWidthArray);min(MSEArray);min(RPArray);min(FPArray)];
Table2(9:12,4)=[median(IntervalWidthArray);median(MSEArray);median(RPArray);median(FPArray)];
Table2(9:12,5)=[max(IntervalWidthArray);max(MSEArray);max(RPArray);max(FPArray)];
clearvars -except Table2