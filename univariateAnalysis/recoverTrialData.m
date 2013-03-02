
% file_name = 'RK-20091230-scan2';

function trial_data = recoverTrialData(file_name)

fh = fopen(file_name, 'r');

trialheader = {'### trialnum:', '%d'};
trialfields = {'trial_struct.task_id',            '%f', 'task_id'
               'trial_struct.task_name',          '%s', 'task_name'
               'trial_struct.wait_FixAcq',        '%f', 'wait_FixAcq'
               'trial_struct.wait_FixAcq2StimOn', '%f', 'wait_FixAcq2StimOn'
               'trial_struct.wait_StimOff2Done',  '%f', 'wait_StimOff2Done'
               'trial_struct.wait_InterTrial',    '%f', 'wait_InterTrial'
               'trial_struct.nofix_timeout',      '%f', 'nofix_timeout'
               'trial_struct.nochoice_timeout',   '%f', 'nochoice_timeout'
               'trial_struct.fixbreak_timeout',   '%f', 'fixbreak_timeout'
               'trial_struct.eyecalib_key',       '%f', 'eyecalib_key'
               'trial_struct.box_response',       '%f', 'box_response'                                                        
               'trial_struct.box_device',         '%f', 'box_device'
               'trial_struct.key_device',         '%f', 'key_device'
               'trial_struct.give_feedback',      '%f', 'give_feedback'
               'trial_struct.stop_feedback_at',   '%f', 'stop_feedback_at'
               'trial_struct.stim_fileName',      '%s', 'stim_fileName'
               'trial_struct.stim_group'          '%f', 'stim_group'
               'trial_struct.stim_seed',          '%f', 'stim_seed'
               'trial_struct.stim_novar_pct',     '%f', 'stim_novar_pct'
               'trial_struct.stim_maxdur',        '%f', 'stim_maxdur'
               'trial_struct.stim_window',        '%f', 'stim_window'
               'trial_struct.fp_pos',             '%f', 'fp_pos' 
               'trial_struct.fp_color',           '%f', 'fp_color'
               'trial_struct.stim_pattern',       '%f', 'stim_pattern'
               'trial_struct.fp_window.type',     '%s', 'fp_window.type'
               'trial_struct.fp_window.spec',     '%f', 'fp_window.spec'
               'trial_struct.recenter_fp_window', '%f', 'recenter_fp_window'
               'trial_struct.boxmap',             '%s', 'boxmap'
               'trial_struct.stim_coh_seq',       '%f', 'stim_coh_seq'
               'trial_struct.stim_coh_time',      '%f', 'stim_coh_time' 
               'trial_struct.stim_pulse_num',     '%f', 'stim_pulse_num' 
               'trial_struct.event_name',         '%s', 'event_name'
               'trial_struct.event_happened',     '%f', 'event_happened' 
               'trial_struct.event_time',         '%f', 'event_time' 
               'trial_struct.response',           '%f', 'response' 
               'trial_struct.response_mode',      '%s', 'response_mode' 
               'trial_struct.result',             '%s', 'result'
               'trial_struct.abort',              '%f', 'abort' 
               'trial_struct.start_t',            '%f', 'start_t' 
               'trial_struct.end_t',              '%f', 'end_t' 
               'trial_struct.rt',                 '%f', 'rt'      
               'trial_struct.cresponse_mode',     '%s', 'cresponse_mode'
               'trial_struct.cresponse',          '%f', 'cresponse'
               'trial_struct.crt',                '%f', 'crt'
               'trial_struct.cresult',            '%s', 'cresult'
               'trial_struct.trial_id',           '%f', 'trial_id' 
               'trial_struct.time0',              '%f', 'time0'  
              };

try
    trial_data = emptyTrial(trialfields);
    while 1,
        line = fgetl(fh);    
        if ~ischar(line),
            break;
        else
            L = strfind(line,trialheader{1});
            if ~isempty(L)
                n = sscanf(line(L+length(trialheader{1}):end), trialheader{2});
                if ~isempty(n) && n>0,
                    [trial_data(n), eof] = getTrial(fh, trialfields, trialheader{1});
                    if eof,
                        break;
                    end;
                end;
            end;    
        end;
    end;
    
    fclose(fh);
    
catch err
    fclose(fh);
    
        %rethrow the error to declare the error that stopped the program
        %note that this should be the last line of the code. Functions
        %following this line weill not be executed. Hey! It is generating
        %an error!
    for i = 1 : length(err.stack),
        fprintf('error in line %d of %s\n', err.stack(i).line, err.stack(i).name);
    end;
    rethrow(err);
end;



%% emptyTrial
    %create an empty struct based on trialfields
function trial_struct = emptyTrial(trialfields)

trial_struct = [];

for i = 1 : size(trialfields,1),
    param_val = [];
        %check for '.' in the field name. if '.' is present we have to deal with a cascade
        %of dynamic fields 
    L = findstr(trialfields{i,3},'.');
    if isempty(L),
        trial_struct.(trialfields{i,3}) = param_val;
    else
        s = 'trial_struct';
        b = 1;
        for j = 1 : length(L),
            s = [s '.(''' trialfields{i,3}(b:L(j)-1) ''')'];
            b = L(j)+1;
        end;
        eval([s '.(''' trialfields{i,3}(L(end)+1:end) ''')=param_val;']);
    end;
end;


        
%% getTrial
    %read the file and populate trials fields
function [trial_struct, eof] = getTrial(fh, trialfields, termination)

eof = 0;
trial_struct = emptyTrial(trialfields);

    %read the file line by line, break when end of file is reached or when a new trial
    %is found in the text file 
while 1,
    line = fgetl(fh);
        %check for end of file
    if ~ischar(line),
        eof = 1;
        break;
    end;
    if ~isempty(line),
            %check for the end of trial_struct
        if ~isempty(strfind(line,termination)),
            fseek(fh, -(length(line)+1), 0);
            break;
        end;
            %assign the line to one of the trial fields, create trials fields dynamically
            %if needed 
        for i = 1 : size(trialfields,1),
                %does the field name match
            L = strfind(line,trialfields{i,1});
            if ~isempty(L),
                L = L(end);     %take the last instance of field name in the line, if the are more than one 
                    %run more test on the field name to ensure we are not using an erroneous match
                    %for a correct match we should have either ':' or '_number' immediately after 
                    %the field name
                if line(L+length(trialfields{i,1})) == ':',
                    cell_ = 0;
                elseif line(L+length(trialfields{i,1}))=='_' && ...
                       ismember(line(L+length(trialfields{i,1})+1),'123456789'), 
                    cell_ = 1;
                else
                    continue;
                end;
                    %get the parameter value. look for the assignment sign ':'
                K = findstr(line(L:end),':');
                if ~isempty(K),
                    K = K(1);   %take the first instance of ':' after the field name  
                    param_val = sscanf(line(L+K:end), trialfields{i,2}); %#ok<NASGU>  suppress unused variable message
                else                                                     %            we use param_val in eval below 
                    param_val = NaN;                                     %#ok<NASGU>
                end;
                    %is a cell check for the counter of the cell field  
                if cell_,
                    cellcounter_ = sscanf(line(L+length(trialfields{i,1})+1:L+K-2), '%d');
                end;
                    %check for '.' in the field name. if '.' is present we have to deal with a cascade
                    %of dynamic fields 
                L = findstr(trialfields{i,3},'.');
                if isempty(L),
                    s = 'trial_struct.(trialfields{i,3})';
                else
                    s = 'trial_struct';
                    b = 1;
                    for j = 1 : length(L),
                        s = [s '.(''' trialfields{i,3}(b:L(j)-1) ''')'];
                        b = L(j)+1;
                    end;
                    s = [s '.(''' trialfields{i,3}(L(end)+1:end) ''')'];
                end;
                    %assign param_val to the field
                if cell_,
                    eval([s '{' num2str(cellcounter_,'%d') '} = param_val;']);
                else
                    eval([s ' = param_val;']);
                end;
            end;
        end;
    end;
end;




