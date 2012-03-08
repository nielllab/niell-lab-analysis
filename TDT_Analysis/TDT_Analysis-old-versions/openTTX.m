function TTX  = openTTX(Tank_Name, Block_Name)
% Matlab codes for reading from TTank
% Jianhua Cang 06-27-03; 10-27/03  
% connect to the tank and read the block
% temporary code, use GUI in the future

% changed to function by cmn 10-10/05

Block_Name


TTX = actxcontrol('ttank.x', [1 1 1 1]);
if (invoke(TTX, 'ConnectServer', 'Local', 'jc') ~= 1)
  err = 'error connecting to server'
end
if (invoke(TTX, 'OpenTank', Tank_Name, 'r') ~= 1)
  err = 'error opening tank'
end
if (invoke(TTX, 'SelectBlock', Block_Name) ~= 1)
  err = 'error selecting block'
end
