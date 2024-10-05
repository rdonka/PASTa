function [] = writeNexVarHeader(fid, varHeader)
% [] = writeNexVarHeader(fid, varHeader) -- write nexVarHeader structure
% to the specified .nex file. 
% 
% INPUT:
%   fid - file id
%           
%   varHeader - nexVarHeader structure

fwrite(fid, varHeader.Type, 'int32');
fwrite(fid, varHeader.Version, 'int32');
writeStringPaddedWithZeros(fid, varHeader.Name, 64);
fwrite(fid, varHeader.DataOffset, 'int32');
fwrite(fid, varHeader.Count, 'int32');    
fwrite(fid, varHeader.WireNumber, 'int32');
fwrite(fid, varHeader.UnitNumber, 'int32');
fwrite(fid, varHeader.Gain, 'int32');
fwrite(fid, varHeader.Filter, 'int32');
fwrite(fid, varHeader.XPos, 'double');
fwrite(fid, varHeader.YPos, 'double');
fwrite(fid, varHeader.WFrequency, 'double');
fwrite(fid, varHeader.ADtoMV, 'double');
fwrite(fid, varHeader.NPointsWave, 'int32');
fwrite(fid, varHeader.NMarkers, 'int32');
fwrite(fid, varHeader.MarkerLength, 'int32');
fwrite(fid, varHeader.MVOffset, 'double');
fwrite(fid, varHeader.PrethresholdTimeInSeconds, 'double');
fwrite(fid, char(zeros(1, 52)), 'char');
