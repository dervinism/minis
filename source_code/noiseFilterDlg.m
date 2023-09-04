function filtering = noiseFilterDlg
% minis helper function for displaying the noise filtering dialogue.

button = questdlg('Band-stop filter the regular (systemic) noise?','Filter Data','Yes','No','No');
if strcmpi(button, 'Yes')
    filtering.state = 'on';
else
    filtering.state = 'off';
end