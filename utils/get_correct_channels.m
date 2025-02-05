function phy = get_correct_channels(data)

%rename structure fields for easier processing
data_cell = struct2cell(data);
channelNames = fieldnames(data);
for f = 1:numel(channelNames)
    channelNames{f} = ['Ch_' num2str(f)];
end
data = cell2struct(data_cell,channelNames);

% find channels for resp, cardiac and trigger
if contains(data.Ch_7.title,'scanner','IgnoreCase',true) || contains(data.Ch_7.title,'trigger','IgnoreCase',true) || contains(data.Ch_7.title,'mri','IgnoreCase',true)  
    phy.cardiac   = data.Ch_1;
    phy.resp      = data.Ch_2;
    phy.trigger   = data.Ch_7;
    fprintf('Channel numbers are assigned as in spike template\n');
else
    %this will be obsolete once everyone is using the spike
    %template
    for i = 1:numel(channelNames)
        if strcmp(data(i).title,'scanner')
            fprintf('scan trigger channel is channel %d\n',i);
            eval(sprintf('phy.trigger = data.Ch%d;',i));
        elseif strcmp(data(i).title,'PULS')
            fprintf('pulse channel is channel %d\n',i);
            eval(sprintf('phy.cardiac = data.Ch%d;',i));
        elseif strcmp(data(i).title,'Resp')
            fprintf('resp channel is channel %d\n',i);
            eval(sprintf('phy.resp = data.Ch%d;',i));
        end
    end
end
