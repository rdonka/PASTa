function [] = writeNex5VarHeader(fid, varHeader)
% [] = writeNex5VarHeader(fid, varHeader) -- write nex5VarHeader structure
% to the specified .nex5 file. 
% 
% INPUT:
%   fid - file id
%           
%   varHeader - nex5VarHeader structure

fwrite(fid, varHeader.Type, 'int32');
fwrite(fid, varHeader.Version, 'int32');
writeStringPaddedWithZeros(fid, varHeader.Name, 64);
fwrite(fid, varHeader.DataOffset, 'uint64');
fwrite(fid, varHeader.Count, 'uint64');    
fwrite(fid, varHeader.TimestampDataType, 'int32');
fwrite(fid, varHeader.ContinuousDataType, 'int32');
fwrite(fid, varHeader.SamplingFrequency, 'double');
writeStringPaddedWithZeros(fid, varHeader.Units, 32);
fwrite(fid, varHeader.ADtoUnitsCoefficient, 'double');
fwrite(fid, varHeader.UnitsOffset, 'double');
fwrite(fid, varHeader.NumberOfDataPoints, 'uint64');
fwrite(fid, varHeader.PrethresholdTimeInSeconds, 'double');
fwrite(fid, varHeader.MarkerDataType, 'int32');
fwrite(fid, varHeader.NumberOfMarkerFields, 'int32');
fwrite(fid, varHeader.MarkerLength, 'int32');
fwrite(fid, varHeader.ContinuousIndexOfFirstPointInFramgmentDataType, 'int32');
fwrite(fid, char(zeros(1, 60)), 'char');


