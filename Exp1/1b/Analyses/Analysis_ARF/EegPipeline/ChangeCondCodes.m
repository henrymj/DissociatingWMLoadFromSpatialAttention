function erp = ChangeCondCodes(erp)
% write the condition codes to the erp.event structure
%
% Note: this can only handle one kind of marker at a time (i.e., stimulus
% or response) because it discards the 'S' or 'R' at the start of the
% string.
%
% Written by Joshua J. Foster, January 29, 2016

erp.eventCodes(1,1)= NaN; %boundary event (i.e., start of recording)
for ii = 2 : length(erp.eventTimes)
    s = erp.event(1,ii).type(2:4); % get characters 2-4 (i.e. the digits). 1 is actually stored as blank-blank-1.
    s = str2double(s);             % covert string to double
    erp.eventCodes(1,ii) = s;      % save the event code structure
end
