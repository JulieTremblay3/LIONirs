%NIRS_AvgTrial_Noise_Close
selection = questdlg('Do you really want to close ?',...
    ' ',...
    'Yes Save & Close','Yes but discard change & Close','No Cancel','Yes Save & Close');
switch selection(1:3),
    case 'Yes'
        if strcmp(selection(5:8),'Save')
            guiHOMER = getappdata(0,'gui_SPMnirsHSJ');
            PMI = get(guiHOMER,'UserData');
            currentsub = 1;
            cf = PMI{currentsub}.currentFile;           
            for filenb=1:numel(PMI{currentsub}.filenoise)
                noise = PMI{currentsub}.filenoise{filenb}.noise;
                vmrk_path = PMI{currentsub}.filenoise{filenb}.name
                [label,ind_dur_ch] = read_vmrk_all(vmrk_path);
                
                %remove old bad steps
                list_bad_step =  [];
                for i=1:size(label,1)
                    if strmatch(label{i,1},'bad_step')
                        list_bad_step = [list_bad_step; i];
                    end
                end
                label(list_bad_step,:)=[];
                ind_dur_ch(list_bad_step,:) = [];
                ind_dur_ch_2 = mat2d2ind_dur_ch(noise');
                label_2 = cell(size(ind_dur_ch_2,1),2);
                label_2(:,1)={'bad_step'};
                label_2(:,2)={'manual'};
                ind_dur_ch = [ind_dur_ch;ind_dur_ch_2];
                label = [label;label_2];
                write_vmrk_all(vmrk_path,ind_dur_ch,label)
            end
            'done'
        end
        
        delete(gcf)
    case 'No '
        return
end


