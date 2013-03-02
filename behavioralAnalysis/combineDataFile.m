

function t = combineDataFile(file_names, file_path)

    %list of fields that must exist in all data files. they will be added to the data if
    %they are missing in a data loaded from a file.
field_names = { 'task_id'
                'task_name'
                'max_trial_num'
                'data_file'
                'wait_FixAcq'
                'wait_FixAcq2DotsOn'
                'wait_DotsOff2Done'
                'wait_InterTrial'
                'nofix_timeout'
                'nochoice_timeout'
                'error_timeout'
                'fixbreak_timeout'
                'abort_key'
                'eyecalib_key'
                'eye_response'
                'box_response'
                'key_response'
                'box_device'
                'key_device'
                'give_feedback'
                'stop_feedback_at'
                'stim_fileName'
                'stim_group'
                'stim_seed'
                'stim_novar_pct'
                'stim_maxdur'
                'stim_window'
                'fp_pos'
                'fp_color'
                'targ_pos'
                'targ_color'
                'ctarg_pos'
                'ctarg_color'
                'stim_pattern'
                'fp_window'
                'recenter_fp_window'
                'targ_window'
                'ctarg_window'
                'keymap'
                'boxmap'
                'stim_coh_seq'
                'stim_coh_time'
                'stim_pulse_num'
                'event_name'
                'event_happened'
                'event_time'
                'response'
                'response_mode'
                'result'
                'abort'
                'start_t'
                'end_t'
                'rt'
                'cresponse_mode'
                'cresponse'
                'crt'
                'cresult'
                'trial_id'
                'eye_data'
                'time0'
              };
          
    %combine all data files 
t = [];
for i = 1 : length(file_names),
        %load the data file
    if isequal(file_names{i}(end-3:end),'.mat'),
        load(fullfile(file_path,file_names{i}));
    else
        trial_data = recoverTrialData(fullfile(file_path,file_names{i}));
    end;   
        %check for missing fields and add them if needed
    for f = 1 : length(field_names),
        if ~isfield(trial_data,field_names{f}),
            trial_data(end).(field_names{f}) = [];
        end;            
    end;
        %add a field to define ownership (scan number) in the combined data set
    for j = 1 : length(trial_data),
        trial_data(j).ownership = i;
        if size(trial_data(j).event_time,1)~=1,
            trial_data(j).event_time = trial_data(j).event_time';
        end;
    end;
        %combine the data
    t = [t, trial_data];                                    %#ok<AGROW>
end;







