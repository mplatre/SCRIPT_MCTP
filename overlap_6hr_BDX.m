%Clears and deletes all variable and plots from MATLAB before running
clear;
close all;

thresh_val  = 0.52;
rm_small = 8000;
rm_large = 150000;

%Prompt user to select folders for experiment and immediately error checks
parent_folder = uigetdir('', 'Select parent folder of genotype and condtion folder');
if parent_folder == 1
    return
end

output_path = uigetdir('', 'Select folder for data output');
if output_path == 0
    return
end

output_name = input('Type in a name for the output movie and images:\n','s');
output_name = strrep(output_name, " ", "_");

tic;

%Child 1 is the name of the folder that contains all of the experiment data
%formatted with underscores and hyphens
child_1 = dir(parent_folder);
dir_child_1 = [child_1.isdir] & ~strcmp({child_1.name},'.') & ~strcmp({child_1.name},'..');
child_1_sub = child_1(dir_child_1);

%Child 2 is the name of the folder that is created in Keyence
child_2 = dir(strcat(parent_folder, '\', child_1_sub.name));
dir_child_2 = [child_2.isdir] & ~strcmp({child_2.name},'.') & ~strcmp({child_2.name},'..');
child_2_sub = child_2(dir_child_2);

%Child 3 is the list of XY01,XY02,XY03,XY04,XY05,XY06
child_3 = dir(strcat(parent_folder, '\', child_1_sub.name, '\', child_2_sub.name));
dir_child_3 = [child_3.isdir] & ~strcmp({child_3.name},'.') & ~strcmp({child_3.name},'..');
child_3_sub = child_3(dir_child_3);

%Child 4 is the list of 145 timesteps for the 12 hour analysis
child_4 = dir(strcat(parent_folder, '\', child_1_sub.name, '\', child_2_sub.name, '\', child_3_sub(1).name));
dir_child_4 = [child_4.isdir] & ~strcmp({child_4.name},'.') & ~strcmp({child_4.name},'..');
child_4_sub = child_4(dir_child_4);

%Save part of the file path to shorten for further use
pre_x = strcat(parent_folder, '\', child_1_sub.name, '\', child_2_sub.name, '\');

%Next we will split up and process the folder names into data for export
parent_split_1 = strsplit(parent_folder, '\');
parent_split_2 = strsplit(parent_split_1{end}, '_');
date = str2num(parent_split_2{1});
experiment_num = parent_split_2{2};
gene_cond = (strsplit(child_1_sub.name, '_')');
split_gene_cond = cellfun( @(gene_cond) strsplit(gene_cond, '-' ), gene_cond, 'UniformOutput', false );

%Initializing variables for storing root data and images
root_vals = {};
t_step = {};
t_step_proc = {};
xfold_t_step = {};
xfold_t_step_proc = {};

for n = 1:length(child_3_sub) %loops through each XY
    
    for i = 1:length(child_4_sub) %loops through each T
        
        collage = [];
        collage_loc = strcat(pre_x, child_3_sub(n).name, '\', child_4_sub(i).name);
        files = dir(strcat(collage_loc,'\**\*.tif'));
        
        %Make the collage of all of the images for each timestep
        for j = 1:length(files)
            fprintf('XY: %d T: %d  F: %d\n', n, i, j);
            image_dir = strcat(collage_loc, '\',  files(j).name);
            imdata = imread(image_dir);
            %Static 30% overlap based on a width of 1920
            %0.3*1920 = 576
            collage = [collage imdata(:,577:1920)];
        end
        
        %Save the collage for later video making
        t_step{i} = collage;
        
        %Perform analysis for each time before moving to next timestep
        
        %Thresholding the image
        BW = imbinarize(collage,'adaptive','ForegroundPolarity','dark','Sensitivity', thresh_val);
        BW_remove = imcomplement(bwareafilt(imcomplement(BW),[rm_small rm_large]));
        
        %Image smoothing to correct for jagged root edges
        N = 21;
        kernel = ones(N, N, N) / N^3;
        blurryImage = convn(double(BW_remove), kernel, 'same');
        sBW = imcomplement(blurryImage < 0.01);
        
        row_pixel = cleanedges(sBW); %calls a function that cleans the edges
        
        [h,w] = size(row_pixel);
        
        %We save the black and white photo for later video making
        t_step_proc{i} = row_pixel;
        
        %Calls a function to convert from binary to a cell array containing
        %connected pixel midpoints
        image_root_pos = con_midpoint(row_pixel, 1);
        
        
        anchor = image_root_pos{1};
        root_vals{i} = anchor;
        
        if i > 1
            %If a root is added or dropped, ignore and use previous starting
            %anchor and update in dataset
            if length(anchor) ~= length(root_vals(i-1))
                anchor = root_vals{i-1};
                root_vals(i) = root_vals(i-1);
            end
        end
        
        for k = 1:length(anchor) %loops for each root
            
            prev_root_center = anchor(k);
            cur_row = 2;
            
            for l = 2:h %loops through all rows of image after first row
                cur_root_center = image_root_pos{l};
                
                %If it reaches the end of all roots, stop
                if isempty(cur_root_center)
                    break
                end
                
                %Loops through all detected segments in row
                for m = 1:length(cur_root_center)
                    delta_center = abs(cur_root_center(m)-prev_root_center);
                    
                    %Detects the end of the specific root in case there
                    %is another object below it
                    if cur_row < l
                        break
                    end
                    
                    %If there is a segment with a center within 20 pixels
                    %of the previous center of the same root
                    if delta_center < 20
                        
                        if delta_center > 1 %special case for large deltaX
                            factor = 2;
                            
                            if cur_root_center(m) < prev_root_center
                                root_data{l-1,k,i} = prev_root_center - delta_center/factor;
                                prev_root_center = cur_root_center(m) - delta_center/factor;
                                cur_row = cur_row + 1;
                            elseif cur_root_center(m) > prev_root_center
                                root_data{l-1,k,i} = prev_root_center + delta_center/factor;
                                prev_root_center = cur_root_center(m) + delta_center/factor;
                                cur_row = cur_row + 1;
                            end
                            break
                        end
                        
                        %Normal case math and stepping through rows for center line
                        root_data{l-1,k, i} = prev_root_center;
                        prev_root_center = cur_root_center(m);
                        cur_row = cur_row + 1;
                        
                        break
                    end
                end
            end
        end
        
        root_data{h,1,1} = []; %ensures consistent height for dataset
        
    end
    
    root_data{h,k,i} = []; %ensures consistent depth for dataset
    
    xfold_t_step{n} = t_step;
    xfold_t_step_proc{n} = t_step_proc;
    xfold_all_root_data{n} = root_data;
    
    clear root_data;
    
end


xfold_all_root_rate = {};
xfold_all_root_length = {};
xfold_all_root_dist = {};
xfold_all_root_overwrite = {};


for n = 1:length(child_3_sub)
    
    root_overwrite = xfold_all_root_data{n};
    [h,root_num,time] = size(xfold_all_root_data{n});
    
    for i = 1:root_num
        for j = 2:time
            prev = root_overwrite(:,i,j-1);
            next = root_overwrite(:,i,j);
            index = find(~cellfun('isempty', prev),1,'last');
            
            updated = vertcat(prev(1:index), next(index+1:end));
            
            if isempty(updated) %When the root is empty, repost it as empty
                updated = prev;
            end
            root_overwrite(:,i,j) = updated;
        end
    end
    
    %Next we need to convert that to a 2-D array of each root (x) and it's
    %distance at each timepoint (y)
    
    all_root_dist = {};
    
    %Next we do the math to calculate the distance based on the center
    for i = 1:root_num
        for j = 1:time
            c_dist = 0;
            for k = 2:h
                if isempty(root_overwrite{k,i,j})
                    all_root_dist{j,i} = 0;
                    break
                end
                c_dist = c_dist + (sqrt(((root_overwrite{k,i,j} - root_overwrite{(k-1),i,j})^2)+1));
            end
            all_root_dist{j,i} = c_dist;
        end
    end
    
    fprintf('\nXY REGION: %d\n', n);
    %Next we remove any roots that jump in length due to algorithm errors
    to_remove = [];
    for i = 1:root_num
        if (sum(diff(cell2mat(all_root_dist(:,i))) > 500) > 0)
            %checks if there is a jump greater than 500
            %save i value for root to remove
            to_remove = [to_remove i];
            fprintf('JUMP: %d\n', i);
        elseif (nnz(~diff(cell2mat(all_root_dist(:,i)))) > 360)
            %checks if there are more than 360 timesteps where the length
            %does not change and saves i for root to remove
            fprintf('STUCK: %d\n', i);
            to_remove = [to_remove i];
        end
    end
    
    %We then subset the cell array using the list previously created
    for i = length(to_remove):-1:1
        all_root_dist(:,to_remove(i)) = [];
        root_overwrite(:,to_remove(i),:) = [];
    end
    
    xfold_all_root_dist{n} = all_root_dist;
    xfold_all_root_overwrite{n} = root_overwrite;
    
    final_dim = size(all_root_dist);
    final_root_num = final_dim(2);
    
    all_root_rate = {};
    all_root_length = {};
    
    for i = 1:final_root_num
        
        for j = 2:time
            all_root_length{j,i} = all_root_dist{j,i}-all_root_dist{1,i};
            all_root_rate{j,i} = all_root_dist{j,i}-all_root_dist{j-1,i};
        end
        
    end
    
    xfold_all_root_rate{n} = all_root_rate;
    xfold_all_root_length{n} = all_root_length;
    
    clear all_root_dist root_overwrite all_root_rate all_root_length
    
end

%Change to specified directory to output files and videos
cd(output_path);

for h = 1:length(xfold_t_step_proc)
    %Outer loop for each XY folder
    
    %Get variables to use in inner loops from the dataset
    [hh,root_num,time] = size(xfold_all_root_data{h});
    
    raw_frames = xfold_t_step{h};
    proc_frames = xfold_t_step_proc{h};
    
    %Initialize a datastructure to store each frame of the video
    F(length(proc_frames)) = struct('cdata',[],'colormap',[]);
    
    for i = 1:length(proc_frames)
        
        %Display the images
        imshow([raw_frames{i}; uint8(255*proc_frames{i})]);
        hold on;
        
        root_overwrite = xfold_all_root_overwrite{h};
        root_data = xfold_all_root_data{h};
        
        all_root_dist = xfold_all_root_dist{h};
        final_dim = size(all_root_dist);
        final_root_num = final_dim(2);
        
        %Plot the lines on the images
        for k = 1:final_root_num
            x = [];
            x = cell2mat(root_overwrite(:,k,i));
            plot(x, 1+hh:(length(x)+hh), 'LineWidth', 1.5);
            text((x(end)+50), (length(x)+hh), ['R', num2str(k)]);
        end
        
        for k = 1:root_num
            x = [];
            x = cell2mat(root_data(:,k,i));
            if isempty(x)
                continue
            end
            plot(x, 1:(length(x)), 'LineWidth', 1.5);
            text((x(end)+50), (length(x)), ['R', num2str(k)]);
        end
        
        %Saves the image with lines as a single frame
        F(i) = getframe(gcf);
        hold off;
    end
    
    close all;
    
    v = VideoWriter(strcat(char(output_name), '_', num2str(h), '.avi'));
    open(v)
    writeVideo(v,F)
    close(v)
    
end

%Plot the length values and save jpg
figure;
hold on
for i = 1:length(xfold_all_root_length)
    all_root_length = xfold_all_root_length{i};
    all_root_dist = xfold_all_root_dist{i};
    final_dim = size(all_root_dist);
    final_root_num = final_dim(2);
    
    for j = 1:final_root_num
        plot(cell2mat(all_root_length(:,j)));
    end
end
title('Root Growth from T=0');
xlabel('Time Step');
ylabel('Center Line Distance in Pixels');
saveas(gcf,strcat(output_name, '.jpg'))
hold off;
close all;


%Plot the average length values and save jpg
figure;
hold on
for i = 1:length(xfold_all_root_length)
    all_root_length = xfold_all_root_length{i};
    all_root_dist = xfold_all_root_dist{i};
    final_dim = size(all_root_dist);
    final_root_num = final_dim(2);
    
    %plot(mean(cell2mat(all_root_length(:,:))'), 'DisplayName', strcat('XY', num2str(i)));
    plot(mean(cell2mat(all_root_length(:,:))'), 'DisplayName', gene_cond{i});

end
legend('Location', 'northwest');
title('Average Root Growth from T=0');
xlabel('Time Step');
ylabel('Center Line Distance in Pixels');
saveas(gcf,strcat(output_name, '_average', '.jpg'))
hold off;
close all;


export_matrix = {};
for h = 1:length(xfold_all_root_length)
    all_root_length = xfold_all_root_length{h};
    all_root_dist = xfold_all_root_dist{h};
    final_dim = size(all_root_dist);
    final_root_num = final_dim(2);
    
    sub_matrix = {};
    
    for i = 1:final_root_num
        for j = 1:length(t_step_proc)
            
            %Date in column 1
            sub_matrix{((i-1)*length(t_step_proc))+j,1} = date;
            %Experiment Number in column 2
            sub_matrix{(((i-1)*length(t_step_proc))+j),2} = experiment_num;
            %Genotype in column 3
            sub_matrix{(((i-1)*length(t_step_proc))+j),3} = split_gene_cond{h}{1};
            %Condition in column 4
            sub_matrix{(((i-1)*length(t_step_proc))+j),4} = split_gene_cond{h}{2};
            %Root number in column 5
            sub_matrix{(((i-1)*length(t_step_proc))+j),5} = i;
            %Timestep in column 6
            sub_matrix{(((i-1)*length(t_step_proc))+j),6} = j*5;
            
            
            if isempty(all_root_length{j,i})
                sub_matrix{(((i-1)*length(t_step_proc))+j),7} = 0;
            else
                sub_matrix{(((i-1)*length(t_step_proc))+j),7} = all_root_length{j,i};
            end
        end
    end
    
    export_matrix = vertcat(export_matrix, sub_matrix);
    
end

T = array2table(export_matrix);

if ~isempty(T) %Only tries to save the table if it exists - prevents error cases from crashing the entire code
    T.Properties.VariableNames = {'Date' 'EXP' 'Genotype' 'Condition' 'Root_Num' 'Timestep_min' 'Measure'};
    writetable(T, strcat(output_name, '.xlsx'), 'Sheet', 1);
else
    fprintf('\nAlgorithm failed to detect any roots successfully.\n');
end


toc;