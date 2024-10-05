function [] = writeStringPaddedWithZeros(fid, theString, totalLength)
% [] = writeStringPaddedWithZeros(fid, theString, length) -- write the
% string padded with zeros to the specified  file. 
% 
% INPUT:
%   fid - file id          
%   theString - string to write
%   totalLength - total number of bytes to be written
if isa(theString, 'char') == 0
   error 'invalid data in nexFile: trying to write non-char as char'
end

if length(theString) > totalLength
    theString = theString(1:totalLength);
end

fwrite(fid, theString, 'char');

if totalLength - length(theString) > 0
    padding = char(zeros(1, totalLength - length(theString)));
    fwrite(fid, padding, 'char'); 
end
