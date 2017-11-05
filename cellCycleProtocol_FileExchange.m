% Nuclear fluorescent stains such as DAPI can produce a fluorescent 
% signal that is proportional to concentration nuclear DNA.  Therefore, the 
% integrated intensity of a nuclei's fluorescent signal can suggest whether a cell is
% in the G1, S, or G2 stage of the cell cycle.  In this script I provide a
% workflow to process images of fluorescently stained nuclei and 
% perform analysis of nuclei integrated intensity to assign cells to 
% specific stages of the cells cycle. The idea for this script is inspired 
% by the publication discussed in the documentation folder of this FileExchange post.

% I would like to acknoldege the use of the expandAxes.m and togglefig.m functions developed 
% by Brett Shoelson, that are used here and provided with the FileExchange download. 

%% STEP 1: Use target image(s) to optimize image processing parameters prior to automated analysis 

% PARAMETERS
% set parameters used in the following functions to segment fluorescent
% nuclei for a target image

rolling_ball = 80;  % used with strel('disk', rolling_ball)
salt_size = 50;     % used with bwareopen(BW, salt_size)
lowhighArea = [50, 5000];   % used with bwpropfilt(bw, 'Area', lowhighArea)
% group names should match the folder names that contain images for each
% treatment condition
groupNames = {'Untreated', 'Aphidicolin', 'Camptothecin', 'Etoposide'};

% IMPORT TARGET
[filename, path] = uigetfile({'*.*'}, 'FileSelector');
I = imread(fullfile(path, filename));

% CORRECT SHADING 

I1 = im2double(I);
se = strel('disk', rolling_ball);
background = imopen(I1, se);
I2 = I1-background;
I3 = imadjust(I2, [min(I2(:)) max(I2(:))], []);

% THRESHOLD

% graythresh provides an initial target value to enter as a second argument 
% to the im2bw function
graylevel = graythresh(I3);  
BW = im2bw(I3, graylevel);  
%BW = imbinarize(I3); after R2016a
% Fill holes
BW = imfill(BW, 'holes');
% remove isolated regions that are not nuclei
BW = bwareaopen(BW, salt_size);
% expand the boundaries slightly by dialating BW
BW = imdilate(BW, strel('disk', 1));
% Clear borders
BW = imclearborder(BW);

% IMAGE SEGMENTATION

% perform distance transform on the binary image
D = -bwdist(~BW);

% Prevent oversegmentation by reducing the number of catchment basins.
mask = imextendedmin(D,2);
D2 = imimposemin(D,mask);
Ld2 = watershed(D2);
BW(Ld2 == 0) = 0;

% ROIs

% Remove regions that are outside of the defined nuclei area
BW = bwpropfilt(BW, 'Area', lowhighArea);

boundaries = bwboundaries(BW, 'noholes');
% Find the boundary for each foreground nuclei
numberOfBoundaries = size(boundaries, 1);
togglefig('Segmented Nuclei')
imshow(I3)
hold on
for k = 1 : numberOfBoundaries
    thisBoundary = boundaries{k};
    plot(thisBoundary(:,2), thisBoundary(:,1), 'g', 'LineWidth', 2);
end
hold off;

%% STEP 2: View detected nuclei for every image to be evaluated (Optional)
% The expandAxes function written by Brett Shoelson allows you to view if settings are
% appropriate for all images in your data set. Change the groupNames 
% variable to match the groups you wish to evaluate.  The
% output is a single figure window that contains all images associated with 
% the group being evaluated. that contains individual images for each
% group.  Outlines of the identified nuclei are
% superimposed on the image.  Clicking on an individual image will
% generate an enlarged figure to make it easier to view
% the segmentation result.

for k = 1:size(groupNames, 2)

    filename = dir(fullfile(pwd, groupNames{k}, '*.tif'));
    fileNumber = size(filename, 1);
    figure;
        
    for i = 1:fileNumber

        I = imread(fullfile(groupNames{k}, filename(i).name));
        
        % CORRECT SHADING 
        I1 = im2double(I);
        se = strel('disk', rolling_ball);
        background = imopen(I1, se);
        I2 = I1-background;
        I3 = imadjust(I2, [min(I2(:)) max(I2(:))], []);
        % create expandAxes to evaluate images
        ax(i) = subplot(ceil(sqrt(fileNumber)), ceil(sqrt(fileNumber)), i);

        % THRESHOLD IMAGE 
        %BW = imbinarize(I3); after R2016a
        graylevel = graythresh(I3);
        BW = im2bw(I3, graylevel);
        % Fill holes
        BW = imfill(BW, 'holes');
        % Remove isolated regions that are not nuclei
        BW = bwareaopen(BW, salt_size);
        % Expand the boundaries slightly by dialating BW
        BW = imdilate(BW, strel('disk', 1));


        % IMAGE SEGMENTATION
        % perform distance transform on the binary image
        D = -bwdist(~BW);

        % Prevent oversegementation
        mask = imextendedmin(D,2);
        D2 = imimposemin(D,mask);
        Ld2 = watershed(D2);
        BW(Ld2 == 0) = 0;
        
        % Clear borders
        BW = imclearborder(BW);
        
        % ROIs
        % Remove regions that are outside of the defined nuclei area
        BW = bwpropfilt(BW, 'Area', lowhighArea);

        imshow(I3)
        set(gcf, 'Name', groupNames{k})
        title(filename(i).name, 'interpreter', 'none', 'FontSize', 12); 
        axis image; 
        hold on;
        boundaries = bwboundaries(BW);
        % Find the boundary for each foreground nuclei
        numberOfBoundaries = size(boundaries, 1);

        % Validate the regions that have been identified by superimposing a
        % green outline around each nuclei
        for m = 1 : numberOfBoundaries
            thisBoundary = boundaries{m};
            plot(thisBoundary(:,2), thisBoundary(:,1), 'g', 'LineWidth', 2);
        end
        hold off;

    end
   
    expandAxes(ax)
    clear ax
end


%% STEP 3: Obtain integrated intensity values for nuclei
% A structure will be produced that contains the normalized integrated intensity
% values for all nuclei identified within each group of images.  

intIntensity = []; % short for Integrated Intensity

for k = 1:size(groupNames, 2)

    filename = dir(fullfile(pwd, groupNames{k}, '*.tif'));
    fileNumber = size(filename, 1);
       
    for i = 1:fileNumber

         
        % code contained in the 'i' for loop is the same as code used for
        % the target image
        I = imread(fullfile(groupNames{k}, filename(i).name));
        
        % CORRECT SHADING
        I1 = im2double(I);    
        I1 = im2double(I1);
        se = strel('disk', rolling_ball);
        background = imopen(I1, se);
        I2 = I1-background;
        I3 = imadjust(I2, [min(I2(:)) max(I2(:))], []);

        % THRESHOLD IMAGE
        graylevel = graythresh(I3);
        BW = im2bw(I3, graylevel);
        % Fill holes
        BW = imfill(BW, 'holes');
        % Remove isolated regions that are not nuclei
        BW = bwareaopen(BW, salt_size);
        % Expand the boundaries slightly by dialating BW
        BW = imdilate(BW, strel('disk', 1));

        % IMAGE SEGMENTATION
        % perform distance transform on the binary image
        D = -bwdist(~BW);

        % Prevent oversegmentation
        mask = imextendedmin(D,2);
        D2 = imimposemin(D,mask);
        Ld2 = watershed(D2);
        BW(Ld2 == 0) = 0;

        % Clear borders
        BW = imclearborder(BW);
        
        % ROIs
        % Remove regions that are outside of the defined nuclei area
        BW = bwpropfilt(BW, 'Area', lowhighArea);
        % Assign a numerical identifier to each forground region
        labeledImage = bwlabel(BW, 8); 
        % Obtain all region measurments for each nuclei that has been labeled
        nucleiMeasurements = regionprops(labeledImage, I3, {'PixelValue'}); 
        numberOfNuclei = size(nucleiMeasurements, 1);

        for m = 1:numberOfNuclei
            thisNucleiIntIntensity = nucleiMeasurements(m).PixelValues;
            intIntensity = [intIntensity, sum(thisNucleiIntIntensity)];
        end
    end
    
    % INTEGRATED INTENSITY NORMALIZATION
    % normalize integrated intensity to values between 0-1 for each group
    numerator = intIntensity - min(intIntensity);
    normalized(k).intIntensity = numerator./(max(intIntensity) - min(intIntensity));
    normalized(k).group = groupNames{k};
    intIntensity = [];

end

%% STEP 4: Generate a histogram illustrating integrated intensity distribution
% This section will generate histograms illustrating the distribution of 
% integrated intensity values for each group.  The histograms will be saved
% to the current directory as .png files.  The histograms can be referenced to determine
% appropriate thresholds used to catagorize cells according to the stage of
% the cell cycle.

binNumber = 50;

for k = 1:size(groupNames, 2)

    figure('Name', groupNames{k})
    % create histogram
    h(k) = histogram(normalized(k).intIntensity, binNumber);
    % alter histogram to display frequency rather than raw counts
    h(k).Normalization = 'probability';
    % apply labels
    ylabel('Fraction of cells');
    xlabel('Normalized DAPI Integrated Intensity');
    title(groupNames{k});
    set(gca, 'fontsize', 12);
    % alter figure dimensions
    fig = gcf;
    fig.Units = 'inches';
    fig.Position = [0 0 6 4.5];
    % save
    print(groupNames{k}, '-dpng');
 
    
end

%% STEP 5: Set normalized intensity values for gating G1, S and G2 cells
% After evaluating the distribution of integrated intensity for cells within 
% the control group you can set specific thresholds to catagorize cells in 
% the G1, S and G2 stages of the cell cycle.

g1 = [0.1, 0.26];
s = [0.26, 0.36];
g2 = [0.36, 0.5];

percent = [];

for k = 1:size(groupNames, 2)
    
    num_cells = length(normalized(k).intIntensity);
    % index into structure to identify G1 cells
    idx = normalized(k).intIntensity > g1(1) & ...
        normalized(k).intIntensity < g1(2);
    G1(k).values = normalized(k).intIntensity(idx);
    % index into structure to identify S cells
    idx = normalized(k).intIntensity >= s(1) & ...
        normalized(k).intIntensity < s(2);
    S(k).values = normalized(k).intIntensity(idx);
    % index into structure to identify G2 cells
    idx = normalized(k).intIntensity >= g2(1) & ...
        normalized(k).intIntensity < g2(2);
    G2(k).values = normalized(k).intIntensity(idx);
    
    percent = [percent; length(G1(k).values)/num_cells, length(S(k).values)/num_cells,...
        length(G2(k).values)/num_cells];
end

%% STEP 6: Generate a table and graph illustrating the proportion of cells in G1, S and G2 stages of the cell cycle

stage_percent = array2table(percent);
stage_percent.Properties.VariableNames = {'G1', 'S', 'G2'};
stage_percent.Properties.RowNames = groupNames;
stage_percent
writetable(stage_percent, 'cellCyclePercentages.csv', 'WriteRowNames', true)

figure;
hb = bar(percent');
set(gca, 'xticklabels', {'G1', 'S', 'G2'})
ylabel('Proportion')
legend(groupNames, 'Location', 'bestoutside')