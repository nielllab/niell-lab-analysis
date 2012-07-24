function [orient freq speed contrast phase TempFreq var1value var2value var3value positionX positionY length eye] = generateVarParams(handles)

Var1 = get(handles.Var1,'Value');
Var2 = get(handles.Var2,'Value');
Var3 = get(handles.Var3,'Value');


Var1_range = str2num(get(handles.Var1Range,'String'));
if Var1>1
    nCond1 = size(Var1_range,2)
else
    nCond1=1;
end

Var2_range = str2num(get(handles.Var2Range,'String'));
if Var2>1
    nCond2 = size(Var2_range,2)
else
    nCond2=1;
end

Var3_range = str2num(get(handles.Var3Range,'String'));
if Var3>1
    nCond3 = size(Var3_range,2)
else
    nCond3=1;
end

Orient0 = str2num(get(handles.Orient0,'String'));
Freq0 = str2num(get(handles.Freq0,'String'));
Speed0 = str2num(get(handles.Speed0,'String'));
Contrast0 = str2num(get(handles.Contrast0,'String'));
Phase0 = str2num(get(handles.Phase0,'String'));
TempFreq0 = str2num(get(handles.TempFreq0,'String'));
PositionX0 =  str2num(get(handles.PositionX0,'String'));
PositionY0 =  str2num(get(handles.PositionY0,'String'));
Length0 =  str2num(get(handles.Length0,'String'));
Eye0 = str2num(get(handles.eyeCond0,'String'));

cond = 0;
for c1 = 1:nCond1
    for c2 = 1:nCond2
        for c3 = 1:nCond3
            cond=cond+1;
            orient(cond)=Orient0;
            freq(cond)=Freq0;
            speed(cond)=Speed0;
            contrast(cond)=Contrast0;
            phase(cond) = Phase0;
            TempFreq(cond)=TempFreq0;
            positionX(cond)=PositionX0;
            positionY(cond)=PositionY0;
            length(cond)=Length0;
            eye(cond) = Eye0;
            if Var1>1
                var1value(cond) = Var1_range(c1);
            else
                var1value(cond) = 0;
            end
            if Var2>1
                var2value(cond)=Var2_range(c2);
            else
                var2value(cond)=0;
            end
            if Var3>1
                var3value(cond)=Var3_range(c3);
            else
                var3value(cond)=0;
            end
            
            if Var1>1
                switch Var1
                    case 2
                        orient(cond)=Var1_range(c1);
                    case 3
                        freq(cond)=Var1_range(c1);
                    case 4
                        speed(cond) =Var1_range(c1);
                    case 5
                        contrast(cond)=Var1_range(c1);
                    case 6
                        TempFreq(cond) = Var1_range(c1);
                    case 7
                        phase(cond)=Var1_range(c1);
                    case 8
                        positionX(cond)=Var1_range(c1);
                    case 9
                        positionY(cond)=Var1_range(c1);
                    case 10
                        eye(cond) = Var1_range(c1);
                end
                if Var2>1
                    switch Var2
                        case 2
                            orient(cond)=Var2_range(c2);
                        case 3
                            freq(cond)=Var2_range(c2);
                        case 4
                            speed(cond) =Var2_range(c2);
                        case 5
                            contrast(cond)=Var2_range(c2);
                        case 6
                            TempFreq(cond) = Var2_range(c2);
                        case 7
                            phase(cond)=Var2_range(c2);
                        case 8
                            positionX(cond)=Var2_range(c2);
                        case 9
                            positionY(cond)=Var2_range(c2);
                        case 10
                            eye(cond)=Var2_range(c2);
                    end
                    if Var3>1
                        switch Var3
                            case 2
                                orient(cond)=Var3_range(c3);
                            case 3
                                freq(cond)=Var3_range(c3);
                            case 4
                                speed(cond) =Var3_range(c3);
                            case 5
                                contrast(cond)=Var3_range(c3);
                            case 6
                                TempFreq(cond) = Var3_range(c3);
                            case 7
                                phase(cond)=Var3_range(c3);
                            case 8
                                positionX(cond)=Var3_range(c3);
                            case 9
                                positionY(cond)=Var3_range(c3);
                            case 10
                                eye(cond)=Var3_range(c3);
                        end
                    end  %if var3
                end  %if var2
            end  %if var1
        end %c3
    end  %c2
end %c1


orient
freq
speed
contrast
phase
TempFreq
positionX
positionY

