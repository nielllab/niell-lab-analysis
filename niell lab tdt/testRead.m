close all
clear all
Block_Name = 'bars8dbw_1'
Tank_Name = '071612_lgn';
max_t = 1000;
max_events = 10^9;
event_code = 'pAll'
ch=1;

TTX = actxcontrol('ttank.x');
if TTX.ConnectServer('Local', 'nlab') ~= 1
  err = 'error connecting to server'
end
if (invoke(TTX, 'OpenTank', Tank_Name, 'R') ~= 1)
  err = 'error opening tank'
end
TTX.SetGlobalV('WavesMemLimit',1e9);

if (invoke(TTX, 'SelectBlock', ['~' Block_Name]) ~= 1)
  err = 'error selecting block'
end

TTX.SetGlobalV('Channel',ch)
TTX.SetGlobalV('T2', max_t)

waves = TTX.ReadWavesV('pAll');

N = TTX.ReadEventsV(   max_events ,event_code,ch, 0, ...
    0,max_t,'All')


invoke(TTX, 'CloseTank');
invoke(TTX, 'ReleaseServer');