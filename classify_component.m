clear all;

% NOTE: If any patterns need to be added/removed/updated: 
% Update the "prepare_slices.m" and "order_groups.m" scripts. 
% In "prepare_slices.m", the pattern name and corresponding slice numbers 
% is added to the exemplar it belongs to
% In "order_groups.m" the numbers are shifted to include the new pattern 
% See files for further details

% close all;
% Matches the incoming network image to the exemplars using template matching

% The exemplars saved in /Exemplars are the mean images (can be easily updated) 
% from the SAN drive at /Templates/WoodwardNetworks/Exemplars/Templates_Current

% This script uses the slices described in the Jan24 classification guide

%% Begin. Change line 13 to component number, and 14 to the to-be-classified network

% network that we want to classify
fprintf('Please select component to classify.\n');
[file,filepath] = uigetfile('*.img;*.nii','Multiselect','on');

num_images = 1;

if isequal(file,0)
   disp('No file selected');
   return;
elseif isa(file,'cell')
    num_images = size(file,2);
end


for img = 1:num_images
    for flip_val = 0:1

        clearvars -except num_exemplars num_images flip_val img file filepath

        if num_images == 1
            curr_img = file;
        else
            curr_img = file{:,img};
        end

        to_match_name = [filepath,curr_img];
        fprintf('Matching input component to template size...\n');
        to_match = match_sizes(to_match_name,'ch2.nii');

        if flip_val == 0
            flip = 'n';
        else
            flip = 'y';
        end

        [~,name] = fileparts(curr_img);

        if upper(flip) == 'Y'
            [to_match,curr_img] = flip_img(to_match,filepath,[name '_resampled']);
            to_match_name = [filepath,curr_img];
        end

        fprintf('Classifying component: %s\n',to_match_name);
        fprintf('Loading exemplars...\n');
        % slices of interest
        slices = prepare_slices();
        % exemplar images

        script_loc = mfilename('fullpath');
        script_path = fileparts(script_loc);

        exemplars = dir([script_path filesep 'Exemplars' filesep '*.nii']);
        exemplar_imgs(size(exemplars,1)) = struct();
        for i = 1:size(exemplars,1)
            ex = [script_path filesep 'Exemplars' filesep exemplars(i).name];
            ex = load_nii(ex);
            ex = ex.img;
            exemplar_imgs(i).img = ex;
        end

        num_exemplars = size(exemplars,1);
        %% One version
        % Produce the average correlation between the incoming network and each 
        % exemplar using the slices associated with each exemplar 

        fprintf('Calculating correlations per exemplar pattern (v1)...\n');
        ex_slice_groups(size(exemplars,1)) = struct();
        for i = 1:size(ex_slice_groups,2)
            j = shift(i);
            ex_slice_groups(i).name = slices(j).N;
        end
        ex_slice_groups = initialize_struct(ex_slice_groups,size(exemplars,1),0);

        flip_ex = zeros(1,num_exemplars);
        flip_val_pos = zeros(1,num_exemplars);
        flip_val_neg = zeros(1,num_exemplars);
        for num_ex = 1:size(exemplars,1)
            fprintf('Calculating correlation to exemplar %s/%s\n',num2str(num_ex),num2str(size(exemplars,1)));
            ex = exemplar_imgs(num_ex).img;
            j = shift(num_ex);

            ex_C_slices = slices(j).C;
            ex_A_slices = slices(j).A;

            [ex_slices,tm_slices] = get_all_slices(1,ex_C_slices,ex_A_slices,ex,to_match);

            orig_match = perform_correlation(1,ex_slice_groups,ex_slices,tm_slices,num_ex);
            flipped_match = perform_correlation(1,ex_slice_groups,ex_slices,tm_slices*(-1),num_ex);
            ex_slice_groups = orig_match;
            if num_ex ~= 8 && num_ex ~= 9
                if abs(flipped_match(num_ex).pos) > abs(orig_match(num_ex).pos)
                    flip_ex(num_ex) = 1;
                    flip_val_pos(num_ex) = flipped_match(num_ex).pos;
                end
            else
                if abs(flipped_match(num_ex).neg) > abs(orig_match(num_ex).neg)
                    flip_ex(num_ex) = 1;
                    flip_val_neg(num_ex) = flipped_match(num_ex).neg;
                end        
            end
        end

        A = 1;

        v1_summary = write_summary(1,ex_slice_groups);

        %% Another version
        % Produce the average correlation between each important pattern (slice group)
        % in the incoming network and the corresponding slices in each exemplar

        fprintf('Calculating correlations over all patterns (v3)...\n');

        C_slice_groups = [];
        A_slice_groups = [];
        for i = 1:size(slices,2)
            C_slice_groups = [C_slice_groups,slices(i).C];
            A_slice_groups = [A_slice_groups,slices(i).A];
        end

        num_slice_groups = size(C_slice_groups,2) + size(A_slice_groups,2);
        groups_ex(num_slice_groups) = struct();

        C_names = [];
        A_names = [];
        for i = 1:size(slices,2)
            names = slices(i).L;
            C_names = [C_names,names(1:size(slices(i).C,2))];
            A_names = [A_names,names(size(slices(i).C,2)+1:end)];
        end

        for i = 1:size(C_names,2)
            groups_ex(i).name = C_names(i);
        end

        for i = 1:size(A_names,2)
            offset = size(C_names,2)+i;
            groups_ex(offset).name = A_names(i);
        end

        groups_ex = initialize_struct(groups_ex,num_slice_groups,zeros(size(exemplars,1),1));

        for group = 1:size(C_slice_groups,2)
            fprintf('Calculating correlations for pattern %s/%s...\n',num2str(group),num2str(num_slice_groups));
            for num_ex = 1:size(exemplars,1)
                ex = exemplar_imgs(num_ex).img;

                C_group = C_slice_groups{group};
                C_slices = [];
                for i = 1:size(C_group,2)
                    C_slices = [C_slices,C_group{i}];
                end

                [ex_slices,tm_slices] = get_all_slices(2,C_slices,[],ex,to_match);
                groups_ex = perform_correlation(3,groups_ex,ex_slices,tm_slices,group,num_ex);
            end
        end

        for group = 1:size(A_slice_groups,2)
            offset = size(C_slice_groups,2) + group;
            fprintf('Calculating correlations for pattern %s/%s...\n',num2str(offset),num2str(num_slice_groups));

            for num_ex = 1:size(exemplars,1)
                ex = exemplar_imgs(num_ex).img;

                A_group = A_slice_groups{group};
                A_slices = [];
                for i = 1:size(A_group,2)
                    A_slices = [A_slices,A_group{i}];
                end

                [ex_slices,tm_slices] = get_all_slices(2,[],A_slices,ex,to_match);
                groups_ex = perform_correlation(3,groups_ex,ex_slices,tm_slices,offset,num_ex);
            end
        end

        v3_summary = write_summary(3,groups_ex);

        comp_summary = get_comp_summary(to_match);
        write_output(v1_summary,0,v3_summary,filepath,curr_img,comp_summary,num_exemplars,flip_ex,flip_val_pos,flip_val_neg,name,flip);
    end
end
    
%% Helper functions
% Refactoring any repeating code sections

%% Initializes the structures used to accumulate results
function output = initialize_struct(input, entries, val)
    for i = 1:entries
        input(i).pos = val;
        input(i).neg = val;
        input(i).tot = val;
        input(i).pc_pos = val;
        input(i).pc_neg = val;
        input(i).pc_tot = val;
    end
    output = input;
    clear input
end

%% Adds zeros to each struct item to dynamically extend it
function output = extend(input,idx)
    output = input;
    output(idx).pos = [output(idx).pos,0];
    output(idx).neg = [output(idx).neg,0];
    output(idx).tot = [output(idx).tot,0];
    output(idx).pc_pos = [output(idx).pc_pos,0];
    output(idx).pc_neg = [output(idx).pc_neg,0];
    output(idx).pc_tot = [output(idx).pc_tot,0];
end

%% Extracts all C and A slices from slices struct

function [Cs_pos,Cs_neg,As_pos,As_neg] = open_slices(slices)
    C_slice_groups_pos = [];
    A_slice_groups_pos = [];
    C_slice_groups_neg = [];
    A_slice_groups_neg = [];
    for i = 1:size(slices,2)
        if i ~= 7
            C_slice_groups_pos = [C_slice_groups_pos,slices(i).C];
            A_slice_groups_pos = [A_slice_groups_pos,slices(i).A];
        else
            C_slice_groups_neg = [C_slice_groups_neg,slices(i).C];
            A_slice_groups_neg = [A_slice_groups_neg,slices(i).A];
        end
    end

    C_open_pos = [];
    for i = 1:size(C_slice_groups_pos,2)
        C_open_pos = [C_open_pos,C_slice_groups_pos{i}];
    end
    A_open_pos = [];
    for i = 1:size(A_slice_groups_pos,2)
        A_open_pos = [A_open_pos,A_slice_groups_pos{i}];
    end
    C_open_neg = [];
    for i = 1:size(C_slice_groups_neg,2)
        C_open_neg = [C_open_neg,C_slice_groups_neg{i}];
    end
    A_open_neg = [];
    for i = 1:size(A_slice_groups_neg,2)
        A_open_neg = [A_open_neg,A_slice_groups_neg{i}];
    end
    
    Cs_pos = [];
    As_pos = [];
    Cs_neg = [];
    As_neg = [];
    for i = 1:size(C_open_pos,2)
        Cs_pos = [Cs_pos,C_open_pos{i}];
    end
    for i = 1:size(A_open_pos,2)
        As_pos = [As_pos,A_open_pos{i}];
    end
    for i = 1:size(C_open_neg,2)
        Cs_neg = [Cs_neg,C_open_neg{i}];
    end
    for i = 1:size(A_open_neg,2)
        As_neg = [As_neg,A_open_neg{i}];
    end
end


%% Gets all the exemplar and to_match slices

function [ex_all_slices,tm_all_slices] = get_all_slices(v,C_slices,A_slices,ex,to_match)
    ex_all_slices = [];
    tm_all_slices = [];
    
    if v == 1
        for sg = 1:size(C_slices,2)
            group = C_slices{sg};
            for s = 1:size(group,2)
                ex_all_slices = [ex_all_slices,squeeze(ex(:,group{s},:))];
                tm_all_slices = [tm_all_slices,squeeze(to_match(:,group{s},:))];
            end
        end
        for sg = 1:size(A_slices,2)
            group = A_slices{sg};
            for s = 1:size(group,2)
                ex_all_slices = [ex_all_slices,squeeze(ex(:,:,group{s}))];
                tm_all_slices = [tm_all_slices,squeeze(to_match(:,:,group{s}))];
            end        
        end
    elseif v == 2
        for s = 1:size(C_slices,2)
            ex_all_slices = [ex_all_slices,squeeze(ex(:,C_slices(1,s),:))];
            tm_all_slices = [tm_all_slices,squeeze(to_match(:,C_slices(1,s),:))];
        end
        for s = 1:size(A_slices,2)
            ex_all_slices = [ex_all_slices,squeeze(ex(:,:,A_slices(1,s)))];
            tm_all_slices = [tm_all_slices,squeeze(to_match(:,:,A_slices(1,s)))];
        end
    elseif v == 3
        if ~isempty(C_slices)
            for s = 1:size(C_slices,2)
                ex_all_slices = [ex_all_slices,squeeze(ex(:,C_slices(1,s),:))];
                tm_all_slices = [tm_all_slices,squeeze(to_match(:,C_slices(1,s),:))];
            end
        end
        if ~isempty(A_slices)
            for s = 1:size(A_slices,2)
                ex_all_slices = [ex_all_slices,squeeze(ex(:,:,A_slices(1,s)))];
                tm_all_slices = [tm_all_slices,squeeze(to_match(:,:,A_slices(1,s)))];
            end
        end
    end
end


%% Performs correlation between slices and updates the variable

function output = perform_correlation(v,input,ex_slices,tm_slices,idx1,idx2)
    output = input;
    [xsize,ysize] = size(ex_slices);
    rows = xsize*ysize;
    
    ex_slices = reshape(ex_slices,rows,1);
    tm_slices = reshape(tm_slices,rows,1);
    ex_nonzeros = ex_slices > 0.0001 | ex_slices < -0.0001;
    tm_nonzeros = tm_slices > 0.0001 | tm_slices < -0.0001;
    either_nonzero = ex_nonzeros & tm_nonzeros;
    non_zero = find(either_nonzero);
    ex_slices = ex_slices(non_zero);
    tm_slices = tm_slices(non_zero);
    
    if v == 1 || v == 2
        tms_pos = tm_slices.*(tm_slices>0);
        tms_neg = tm_slices.*(tm_slices<0);
        exs_pos = ex_slices.*(ex_slices>0);
        exs_neg = ex_slices.*(ex_slices<0);
        % use corr2 to correlate and sum up
        if (size(nonzeros(tms_pos),1) > 0) && (size(nonzeros(exs_pos),1) > 0)
            cp = corr2(tms_pos,exs_pos);
            output(idx1).pos = cp;
            if cp > 0
                output(idx1).pc_pos = cp;
            end
        end
        if (size(nonzeros(tms_neg),1) > 0) && (size(nonzeros(exs_neg),1) > 0)
            cn = corr2(tms_neg,exs_neg);
            output(idx1).neg = cn;
            if cn > 0
                output(idx1).pc_neg = cn;
            end
        end
        if (size(nonzeros(tm_slices),1) > 0) && (size(nonzeros(ex_slices),1) > 0)
            ct = corr2(tm_slices,ex_slices);
            output(idx1).tot = ct;
            if ct > 0
                output(idx1).pc_tot = ct;
            end
        end
    elseif v == 3
        tms_pos = tm_slices.*(tm_slices>0);
        tms_neg = tm_slices.*(tm_slices<0);
        exs_pos = ex_slices.*(ex_slices>0);
        exs_neg = ex_slices.*(ex_slices<0);
        % use corr2 to correlate and sum up
        if (size(nonzeros(tms_pos),1) > 0) && (size(nonzeros(exs_pos),1) > 0)
            cp = corr2(tms_pos,exs_pos);
            output(idx1).pos(idx2) = cp;
            if cp > 0
                output(idx1).pc_pos(idx2) = cp;
            end
        end
        if (size(nonzeros(tms_neg),1) > 0) && (size(nonzeros(exs_neg),1) > 0)
            cn = corr2(tms_neg,exs_neg);
            output(idx1).neg(idx2) = cn;
            if cn > 0
                output(idx1).pc_neg(idx2) = cn;
            end
        end
        if (size(nonzeros(tm_slices),1) > 0) && (size(nonzeros(ex_slices),1) > 0)
            ct = corr2(tm_slices,ex_slices);
            output(idx1).tot(idx2) = ct;
            if ct > 0
                output(idx1).pc_tot(idx2) = ct;
            end
        end
    end
end

%% Shifts index to match 8 slice types to 10 exemplars (1h vs 2h response, DMNs share slices)
function j = shift(i)
    if i > 8
        j = i-2;
    elseif i > 6
        j = i-1;
    else
        j = i;
    end
end

%% Writes summary into given file
function summary = write_summary(version,input_struct)

    c(size(input_struct,2)) = struct();

    for i = 1:size(input_struct,2)
        if version == 1 || version == 2
            if i < 6 || i > 9
                c(i).name = input_struct(i).name;
            else
                c(i).name = getname(i);
            end
        elseif version == 3
            c(i).name = input_struct(i).name;
        end
        
        c(i).pos = input_struct(i).pos;
        c(i).neg = input_struct(i).neg;
        c(i).tot = input_struct(i).tot;
        c(i).pc_pos = input_struct(i).pc_pos;
        c(i).pc_neg = input_struct(i).pc_neg;
        c(i).pc_tot = input_struct(i).pc_tot;
    end
    
    summary = c;
end

%% Writes copy-paste-able table into files
% one "sheet" per pos/neg combination
% copy-paste these rows into the real table
function write_output(v1,v2,v3,path,file,comp,num_exemplars,flip_ex,flip_val_pos,flip_val_neg,name,flip)

    [~,component,~] = fileparts(file);

    % converting R to Z scores
    v1 = z_convert(1,v1);
    v3 = z_convert(3,v3);
    
    flip_val_pos = z_convert(0,flip_val_pos);
    flip_val_neg = z_convert(0,flip_val_neg);
    
    perc_pos = num2str(round(comp.numpos/(comp.numpos+comp.numneg)*100,2));
    perc_neg = num2str(round(comp.numneg/(comp.numpos+comp.numneg)*100,2));
    
    if upper(flip) == 'N'
        f = fopen([path name '_match_pos.txt'],'w');
        fprintf(f,'Number of positive voxels in component:, %s (%s%%)\n',addComma(comp.numpos),perc_pos);
        fprintf(f,'Range of positive loadings in component:, [%s,%s]\n',num2str(comp.minpos),num2str(comp.maxpos));
        fprintf(f,'Number of negative voxels in component:, %s (%s%%)\n',addComma(comp.numneg),perc_neg);
        fprintf(f,'Range of negative loadings in component:, [%s,%s]\n\n',num2str(comp.minneg),num2str(comp.maxneg));
    elseif upper(flip) == 'Y'
        f = fopen([path name '_match_pp_negl.txt'],'w');
        fprintf(f,'Number of positive voxels in component:, %s (%s%%)\n',addComma(comp.numneg),perc_neg);
        fprintf(f,'Range of positive loadings in component:, [%s,%s]\n',num2str(-comp.maxneg),num2str(-comp.minneg));
        fprintf(f,'Number of negative voxels in component:, %s (%s%%)\n',addComma(comp.numpos),perc_pos);
        fprintf(f,'Range of negative loadings in component:, [%s,%s]\n\n',num2str(-comp.maxpos),num2str(-comp.minpos));
    end
    fprintf(f,'Component file: %s\n',[path,file]);
    order = order_groups('pos');

    v1_pos = vertcat(v1.pos);
    v1_name = vertcat(v1.name);
    v1_pos = [v1_pos(1:7);v1_pos(10:num_exemplars)];
    v1_name = [v1_name(1:7);v1_name(10:num_exemplars)];
    [v1_pos_sorted,pos_idx] = sort(v1_pos,'descend');
    v1_name_sorted = v1_name(pos_idx);
    fprintf(f,'\nFisher Z based on exemplar groupings\n\n'); %(v1)
    [v1_best_pos,max_ex] = max(v1_pos);
    fprintf(f,'Best match,%s\n',char(v1_name(max_ex)));
    v1_best_match_pos = char(v1_name(max_ex));
    for i = 1:size(vertcat(v1.name),1)
        if i < 8
            fprintf(f,'%s,%.2f\n',char(v1_name_sorted(i)),v1_pos_sorted(i));
        elseif i > 9
            fprintf(f,'%s,%.2f\n',char(v1_name_sorted(i-2)),v1_pos_sorted(i-2));
        end
    end

    % write all pattern matches
    fprintf(f,'\nFisher Z based on pattern groupings\n\n'); %(v3)
    ex_line(f,'pos');
    v3_name = vertcat(v3.name);
    for i = 1:size(order,2)
        v3_pos = v3(order(i)).pos;
        v3_pos = [v3_pos(1:7);v3_pos(10:num_exemplars)];
        listout = sprintf('%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f',v3_pos);
        fprintf(f,'%s,%s\n',char(v3_name(order(i))),listout);
    end
    fclose(f);
    
    if upper(flip) == 'N'
        f = fopen([path name '_match_neg.txt'],'w');
        fprintf(f,'Number of positive voxels in component:, %s (%s%%)\n',addComma(comp.numpos),perc_pos);
        fprintf(f,'Range of positive loadings in component:, [%s,%s]\n',num2str(comp.minpos),num2str(comp.maxpos));
        fprintf(f,'Number of negative voxels in component:, %s (%s%%)\n',addComma(comp.numneg),perc_neg);
        fprintf(f,'Range of negative loadings in component:, [%s,%s]\n\n',num2str(comp.minneg),num2str(comp.maxneg));
    elseif upper(flip) == 'Y'
        f = fopen([path name '_match_np_posl.txt'],'w');
        fprintf(f,'Number of positive voxels in component:, %s (%s%%)\n',addComma(comp.numneg),perc_neg);
        fprintf(f,'Range of positive loadings in component:, [%s,%s]\n',num2str(-comp.maxneg),num2str(-comp.minneg));
        fprintf(f,'Number of negative voxels in component:, %s (%s%%)\n',addComma(comp.numpos),perc_pos);
        fprintf(f,'Range of negative loadings in component:, [%s,%s]\n\n',num2str(-comp.maxpos),num2str(-comp.minpos));
    end
    fprintf(f,'Component file: %s\n',[path,file]);
    order = order_groups('neg');

    v1_neg = vertcat(v1.neg);
    v1_name = vertcat(v1.name);
    v1_neg = v1_neg(8:9);
    v1_name = v1_name(8:9);
    [v1_neg_sorted,v1_neg_idx] = sort(v1_neg,'descend');
    v1_name_sorted = v1_name(v1_neg_idx);
    fprintf(f,'\nFisher Z based on exemplar groupings\n\n'); %(v1)
    [v1_best_neg,max_ex] = max(v1_neg);
    v1_best_match_neg = char(v1_name(max_ex));
    fprintf(f,'Best match,%s\n',char(v1_name(max_ex)));
    for i = 1:size(v1_neg,1)
        fprintf(f,'%s,%.2f\n',char(v1_name_sorted(i)),v1_neg_sorted(i));
    end

    % write all pattern matches
    fprintf(f,'\nFisher Z based on pattern groupings\n\n'); %(v3)
    ex_line(f,'neg');
    for i = 1:size(order,2)
        v3_neg = v3(order(i)).neg;
        v3_neg = v3_neg(8:9);
        listout = sprintf('%.2f,%.2f',v3_neg);
        fprintf(f,'%s,%s\n',char(v3_name(order(i))),listout);
    end
    fclose(f);
    
    merged = vertcat(v1_pos_sorted,v1_neg_sorted);
    ymax = max(merged)+0.2;
    ymin = min(merged)-0.2;
    
    fprintf('Component classification complete.\n\n');
    
    fprintf('Classifying component: %s\n',file);
    fprintf('\nNumber of positive voxels in component: %s (%s%%)\n', addComma(comp.numpos),perc_pos);
    fprintf('Best match: %s, with Z = %s\n', v1_best_match_pos, num2str(v1_best_pos)); 
    fprintf('\nNumber of negative voxels in component: %s (%s%%)\n', addComma(comp.numneg),perc_neg);
    fprintf('Best match: %s, with Z = %s\n', v1_best_match_neg, num2str(v1_best_neg)); 

    fprintf('Done\n\n');
end

%% Convert R score to Z score

function z_val = z_convert(v,v_vals)
    [x,y] = size(v_vals);
    z_val = v_vals;
    
    if v == 1 || v == 2
        for i = 1:y
            z_val(i).pos = 0.5*(log(1+v_vals(i).pos) - log(1-v_vals(i).pos));
            z_val(i).neg = 0.5*(log(1+v_vals(i).neg) - log(1-v_vals(i).neg));
            z_val(i).tot = 0.5*(log(1+v_vals(i).tot) - log(1-v_vals(i).tot));
            z_val(i).pc_pos = 0.5*(log(1+v_vals(i).pc_pos) - log(1-v_vals(i).pc_pos));
            z_val(i).pc_neg = 0.5*(log(1+v_vals(i).pc_neg) - log(1-v_vals(i).pc_neg));
            z_val(i).pc_tot = 0.5*(log(1+v_vals(i).pc_tot) - log(1-v_vals(i).pc_tot));
        end
    elseif v == 3
        for i = 1:y
            z_val(i).pos = 0.5*(log(1+v_vals(i).pos) - log(1-v_vals(i).pos));
            z_val(i).neg = 0.5*(log(1+v_vals(i).neg) - log(1-v_vals(i).neg));
            z_val(i).tot = 0.5*(log(1+v_vals(i).tot) - log(1-v_vals(i).tot));
            z_val(i).pc_pos = 0.5*(log(1+v_vals(i).pc_pos) - log(1-v_vals(i).pc_pos));
            z_val(i).pc_neg = 0.5*(log(1+v_vals(i).pc_neg) - log(1-v_vals(i).pc_neg));
            z_val(i).pc_tot = 0.5*(log(1+v_vals(i).pc_tot) - log(1-v_vals(i).pc_tot));
        end
    elseif v == 0
        for i = 1:y
            z_val(i) = 0.5*(log(1+v_vals(i)) - log(1-v_vals(i)));
        end
    end
end

%% Adds a line with list of exemplar names to file

function ex_line(f,type)
    if strcmp(type,'pos')
            fprintf(f,',CE,LANG,INT,EXT,INIT,RESP1,RESP2,AUD,AAR,FVF\n');
    elseif strcmp(type,'neg')
        fprintf(f,',DMNT, DMNN\n');
    end
end

%% Gets name of exemplar

function name = getname(i)
    if i == 6
        name = {'One-handed Response'};
    elseif i == 7
        name = {'Two-handed Response'};
    elseif i == 8
        name = {'Traditional DMN'};
    elseif i == 9
        name = {'Novel DMN'};
    end
end

%% Get min, max values and num pos, neg voxels

function comp_summary = get_comp_summary(component)

    comp_summary = struct();
    cshape = size(component);
    component = reshape(component,[cshape(1)*cshape(2)*cshape(3),1]);
    c_positives = component > 0.0001;
    c_negatives = component < -0.0001;
    comp_summary.total = cshape(1)*cshape(2)*cshape(3);
    comp_summary.maxpos = max(component);
    comp_summary.minpos = min(component(c_positives));
    comp_summary.maxneg = max(component(c_negatives));
    comp_summary.minneg = min(component);
    comp_summary.numpos = sum(c_positives);
    comp_summary.numneg = sum(c_negatives);
    
end

%% Adding commas to numbers

function numOut = addComma(numIn)
        jf=java.text.DecimalFormat; % comma for thousands, three decimal places
        numOut=char(jf.format(numIn)); % omit "char" if you want a string out
end

