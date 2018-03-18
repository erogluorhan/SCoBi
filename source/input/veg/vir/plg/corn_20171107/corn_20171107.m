classdef corn_20171107 < handle
    %CORNGENERATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        stage;
        
        % Dates of measurement which starts/finishes a stage
        cornStageDates =  ['6/11/12'; '6/26/12'; '7/17/12'; '8/15/12'; ...
                            '9/19/12'; '9/25/12'];
                        
        cornStagesStruct = struct('s1', 'V1-V9', ...
                                    's2', 'V10-VT', ...
                                    's3', 'R1-R4', ...
                                    's4', 'R5', ...
                                    's5', 'R6') ;
                                
        % TO-DO: Below values are to be taken by a probability distribution 
        leafAzimuthAngleDegMin = 0;        % The min leaf azimuth angle in degrees
        leafAzimuthAngleDegMax = 360;      % The max leaf azimuth angle in degrees
        cobBeginAngleDegMin = 20;         % The min cob angle from the stem (zenith) in degrees
        cobBeginAngleDegMax = 40;         % The max cob angle from the stem (zenith) in degrees
        cobAzimuthAngleDegMin = 0;        % The min cob azimuth angle in degrees
        cobAzimuthAngleDegMax = 360;      % The max cob azimuth angle in degrees
       
        pdfStemHeight;
        pdfStemBottomDiameter;
        pdfStemTopDiameter;
        pdfStemDielRe;
        pdfStemDielIm;
        
        pdfLeafNumbers; 
        pdfLeafLength;
        pdfLeafWidth;
        pdfLeafThickness;
        pdfLeafBeginningAngle; 
        pdfLeafDielRe; 
        pdfLeafDielIm; 
        
        pdfCobNumbers;
        pdfCobLength; 
        pdfCobDiameter; 
        pdfCobDielRe; 
        pdfCobDielIm;

    end
    
    methods
        
        % The initializor function that is directly called from the
        % simulator
        function initialize(obj)
            
            obj.stage = obj.cornStagesStruct.s1;
            
            % Set Vegetation (Corn) Generator 
            obj.setCornGenerator;

        end    
        
        function isTheSame = isTheSame(obj, objOther)
        
            
%         leafAzimuthAngleDegMax = 360;      % The max leaf azimuth angle in degrees
%         cobBeginAngleDegMin = 20;         % The min cob angle from the stem (zenith) in degrees
%         cobBeginAngleDegMax = 40;         % The max cob angle from the stem (zenith) in degrees
%         cobAzimuthAngleDegMin = 0;        % The min cob azimuth angle in degrees
%         cobAzimuthAngleDegMax = 360;
            
            isTheSame = 0;
            
            if strcmp( obj.stage, objOther.stage )
                if obj.leafAzimuthAngleDegMin == objOther.leafAzimuthAngleDegMin
                    if obj.leafAzimuthAngleDegMin == objOther.leafAzimuthAngleDegMin
                        if obj.leafAzimuthAngleDegMin == objOther.leafAzimuthAngleDegMin
                            if obj.leafAzimuthAngleDegMin == objOther.leafAzimuthAngleDegMin
                                if obj.leafAzimuthAngleDegMin == objOther.leafAzimuthAngleDegMin
                                    if obj.leafAzimuthAngleDegMin == objOther.leafAzimuthAngleDegMin
                                        isTheSame = 1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        
        end
        
        function out = get.stage(obj)
            out = obj.stage;
        end          
        
        % The plant generator function that is directly called by the
        % simulator
        function [plantHeight, particleData] = generatePlant( obj, TYPES, position )
            
            % Set Vegetation (Corn) Generator 
            [plantHeight, particleData] = obj.generateCorn( TYPES, position );
            
        end
        
        
        function [plantHeight, particleData] = generateCorn( obj, TYPES, position)
            
            particleData = [];
            
            plantHeight = 0;

            % TO-DO: Work on this. Finalize & Justify vegetation algorithms
            [stemLength, stemBottomDiam, stemTopDiam, epsrStemRe, epsrStemIm] = obj.getRandomStemParams();

            stemCOG = obj.calculateCenterOfGravity(position, 0, 0, stemLength/2);
            
            particleData = [particleData; TYPES.T, 1, stemCOG', 0, 0, ...
                stemBottomDiam/2, stemTopDiam/2, stemLength, epsrStemRe, epsrStemIm];

            % Update plant height w.r.t stem heights
            if stemLength > plantHeight
                plantHeight = stemLength;
            end

            numberLeaves = obj.getRandomLeafNumbers();

            for ii = 1 : numberLeaves

                % TO-DO: Determine realistic azimuth angle interval
                leafAzimuthAngle = obj.leafAzimuthAngleDegMin + ( obj.leafAzimuthAngleDegMax - obj.leafAzimuthAngleDegMin ) * rand;

                % TO-DO: Work on this. Finalize & Justify vegetation algorithms
                [leafLength, leafWidth, leafThickness, leafBeginAngle, epsrLeafRe, epsrLeafIm ] = obj.getRandomLeafParams();

                % TO-DO: Determine a realistic method
                % Leaf z-position is taken as its beginning point (Connection point to the stem)
                if ii <= 2  % Any corn has at least 2 leaves at the top
                    position(3) = stemLength;
                else
                    position(3) = 0.3 * stemLength + ( 0.9 - 0.3 ) * stemLength * rand;  % Leaves can be between [0.3,0.9] of the stem's height
                end

                leafCOG = obj.calculateCenterOfGravity(position, leafBeginAngle, leafAzimuthAngle, leafLength/2);

                % Update top layer height w.r.t leaf heights
                leafTopHeight = position(3) + leafLength * cos(degtorad(leafBeginAngle)); % the top height of the leaf

                if leafTopHeight > plantHeight
                    plantHeight = leafTopHeight; 
                end

                particleData = [particleData; TYPES.L, 1, leafCOG', leafBeginAngle, ...
                    leafAzimuthAngle, leafLength/2, leafWidth/2, leafThickness, ...
                    epsrLeafRe, epsrLeafIm];
            end

            % TO-DO: Check cob positioning
            numberCobs = obj.getRandomCobNumbers();

            for ii = 1 : numberCobs
                cobBeginAngle = obj.cobBeginAngleDegMin + ( obj.cobBeginAngleDegMax - obj.cobBeginAngleDegMin ) * rand;
                cobAzimuthAngle = obj.cobAzimuthAngleDegMin + ( obj.cobAzimuthAngleDegMax - obj.cobBeginAngleDegMin ) * rand;

                [ cobLength, cobWidth, epsrCobRe, epsrCobIm ] = obj.getRandomCobParams();

                % TO-DO: Magic numbers below!!!
                % Cob z-position is taken as its beginning point (Connection point to the stem)
                position(3) = 0.35 * stemLength + ( 0.75 - 0.35 ) * stemLength * rand; % Cobs can be between [0.35,0.75] of the stem's height
                
                cobCOG = obj.calculateCenterOfGravity(position, cobBeginAngle, cobAzimuthAngle, cobLength/2);

                % Update top layer height w.r.t cob heights
                cobTopHeight = position(3) + cobLength * cos(degtorad(cobBeginAngle)); % the top height of the cob

                if cobTopHeight > plantHeight
                    plantHeight = cobTopHeight; 
                end

                particleData = [particleData; TYPES.B, 1, cobCOG', cobBeginAngle, ...
                    cobAzimuthAngle, cobWidth/2, cobWidth/2, cobLength, ...
                    epsrCobRe, epsrCobIm];
            end
        end
        
        
        function setCornGenerator(obj)


            if strcmp( obj.stage, obj.cornStagesStruct.s1 )
                stageStartIndex = 1;
            elseif strcmp( obj.stage, obj.cornStagesStruct.s2 )
                stageStartIndex = 2;
            elseif strcmp( obj.stage, obj.cornStagesStruct.s3 )
                stageStartIndex = 3;
            elseif strcmp( obj.stage, obj.cornStagesStruct.s4 )
                stageStartIndex = 4;
            elseif strcmp( obj.stage, obj.cornStagesStruct.s5 )
                stageStartIndex = 5;
            end

            stageStart = obj.cornStageDates( stageStartIndex, : );
            stageFinish = obj.cornStageDates( stageStartIndex + 1, : );

            % Take all the dates the of full corn growth
            pathname = fileparts(mfilename('fullpath'));
            filename = strcat( pathname,'\corn_data\size\Dates.txt');
            fid = fopen(filename);
            dates = textscan(fid, '%s');
            fclose(fid);

            % Take the date interval for the corn growth stage of interest
            dates = dates{1,1};
            [N, ~] = size(dates);
            index1 = strfind(dates, stageStart);
            index1 = find(not(cellfun('isempty', index1)));
            index2 = strfind(dates, stageFinish);
            index2 = find(not(cellfun('isempty', index2)));

            % Determine indices for the whole data and that of interest
            indicesAll = (1 : N)';
            indices = (index1:1:index2 )';
            indicesAll_poly = ( min(indicesAll) : 0.1 : max( indicesAll ) )';
            
            obj.pdfStemHeight = obj.fitdistStemHeight( indices, indicesAll, indicesAll_poly );
            obj.pdfStemBottomDiameter = obj.fitdistStemBottomDiameter(indices, indicesAll, indicesAll_poly );
            obj.pdfStemTopDiameter = obj.fitdistStemTopDiameter(indices, indicesAll, indicesAll_poly );
            obj.pdfStemDielRe = obj.fitdistStemDielRe(indices, indicesAll, indicesAll_poly );
            obj.pdfStemDielIm = obj.fitdistStemDielIm(indices, indicesAll, indicesAll_poly );
            obj.pdfLeafNumbers = obj.fitdistLeafNumbers(indices);
            obj.pdfLeafLength = obj.fitdistLeafLength(indices, indicesAll, indicesAll_poly );
            obj.pdfLeafWidth = obj.fitdistLeafWidth(indices, indicesAll, indicesAll_poly );
            obj.pdfLeafThickness = obj.fitdistLeafThickness(indices, indicesAll, indicesAll_poly );
            obj.pdfLeafBeginningAngle = obj.fitdistLeafBeginningAngle(indices, indicesAll, indicesAll_poly );
            obj.pdfLeafDielRe = obj.fitdistLeafDielRe(indices, indicesAll, indicesAll_poly );
            obj.pdfLeafDielIm = obj.fitdistLeafDielIm(indices, indicesAll, indicesAll_poly );
            obj.pdfCobNumbers = obj.fitdistCobNumbers(indices);
            obj.pdfCobLength = obj.fitdistCobLength(indices, indicesAll, indicesAll_poly );
            obj.pdfCobDiameter = obj.fitdistCobDiameter(indices, indicesAll, indicesAll_poly );
            obj.pdfCobDielRe = obj.fitdistCobDielRe(indices, indicesAll, indicesAll_poly );
            obj.pdfCobDielIm = obj.fitdistCobDielIm(indices, indicesAll, indicesAll_poly );

            end


            function [pdStemHeight] = fitdistStemHeight(obj, indices, indicesAll, indicesAll_poly )

            %% STEM HEIGHT
            % Take the stem heights w.r.t. the date interval
            pathname = fileparts(mfilename('fullpath'));
            filename = strcat( pathname,'\corn_data\size\StemHeight.txt');
            fid = fopen(filename);
            stemHeightsAll = fscanf(fid,'%f');
            fclose(fid);

            stemHeightsAll = stemHeightsAll * 0.01; % To convert from cm to meters

            % Fit the 3rd-degree polynomial into the stem height data
            poly = polyfit( indicesAll, stemHeightsAll, 3);
            stemHeightsAll_poly = poly(1)*indicesAll_poly.^3 + poly(2)*indicesAll_poly.^2 + ...
                poly(3)*indicesAll_poly + poly(4);
            stemHeights_poly = stemHeightsAll_poly( find(indicesAll_poly == indices(1)) : find(indicesAll_poly == indices(end)), 1 );

            % Fit Gaussian Distribution and draw it
            pdStemHeight = fitdist( stemHeights_poly, 'normal');
            
            end


            function [pdNormalStemBottomDiameter] = fitdistStemBottomDiameter(obj, indices, indicesAll, indicesAll_poly )

            %% STEM BOTTOM DIAMETERS
            % Take the stem bottom diameters w.r.t. the date interval
            pathname = fileparts(mfilename('fullpath'));
            filename = strcat( pathname,'\corn_data\size\StemDiameter1.txt');
            fid = fopen(filename);
            stemBottomDiamsAll = fscanf(fid,'%f');
            fclose(fid);

            stemBottomDiamsAll = stemBottomDiamsAll * 0.001; % To convert from mm to meters

            % Fit the 3rd-degree polynomial into the stem bottom diameters data
            poly = polyfit(indicesAll, stemBottomDiamsAll, 3);
            stemBottomDiamsAll_poly = poly(1)*indicesAll_poly.^3 + poly(2)*indicesAll_poly.^2 + ...
                                poly(3)*indicesAll_poly + poly(4);
            stemBottomDiams_poly = stemBottomDiamsAll_poly( find(indicesAll_poly == indices(1)) : find(indicesAll_poly == indices(end)), 1 );

            % Fit Gaussian Distribution and draw it
            pdNormalStemBottomDiameter = fitdist(stemBottomDiams_poly,'Normal');


            end


            function [pdNormalStemTopDiameter] = fitdistStemTopDiameter(obj, indices, indicesAll, indicesAll_poly )

            %% STEM TOP DIAMETERS
            % Take the stem top diameters w.r.t. the date interval
            pathname = fileparts(mfilename('fullpath'));
            filename = strcat( pathname,'\corn_data\size\StemDiameter2.txt');
            fid = fopen(filename);
            stemTopDiamsAll = fscanf(fid,'%f');
            fclose(fid);

            stemTopDiamsAll = stemTopDiamsAll * 0.001; % To convert from mm to meters

            % Fit the 3rd-degree polynomial into the stem Top diameters data
            poly = polyfit( indicesAll, stemTopDiamsAll, 3);
            stemTopDiamsAll_poly = poly(1)*indicesAll_poly.^3 + poly(2)*indicesAll_poly.^2 + ...
                                poly(3)*indicesAll_poly + poly(4);
            stemTopDiams_poly = stemTopDiamsAll_poly( find(indicesAll_poly == indices(1)) : find(indicesAll_poly == indices(end)), 1 );

            % Fit Gaussian Distribution and draw it
            pdNormalStemTopDiameter = fitdist( stemTopDiams_poly, 'Normal' );

            end


            function [pdNormalStemDielRe] = fitdistStemDielRe( obj, indices, indicesAll, indicesAll_poly )

            %% STEM DIELECTRIC CONSTANT REAL PART
            % Take the stem dielectric constant real part w.r.t. the date interval
            pathname = fileparts(mfilename('fullpath'));
            filename = strcat( pathname,'\corn_data\diel\Stems_RE.txt');
            fid = fopen(filename);
            stemDielReAll = fscanf(fid,'%f');
            fclose(fid);

            % Fit the 3rd-degree polynomial into the stem dielectric constants real data
            poly = polyfit( indicesAll, stemDielReAll, 3);
            stemDielReAll_poly = poly(1)*indicesAll_poly.^3 + poly(2)*indicesAll_poly.^2 + ...
                                poly(3)*indicesAll_poly + poly(4);
            stemDielRe_poly = stemDielReAll_poly( find(indicesAll_poly == indices(1)) : find(indicesAll_poly == indices(end)), 1 );

            % Fit Gaussian Distribution and draw it
            pdNormalStemDielRe = fitdist( stemDielRe_poly, 'Normal' );

            end


            function [pdNormalStemDielIm] = fitdistStemDielIm( obj, indices, indicesAll, indicesAll_poly )

            %% STEM DIELECTRIC CONSTANT IMAGINARY PART
            % Take the stem dielectric constant imaginary part w.r.t. the date interval
            pathname = fileparts(mfilename('fullpath'));
            filename = strcat( pathname,'\corn_data\diel\Stems_IM.txt');
            fid = fopen(filename);
            stemDielImAll = fscanf(fid,'%f');
            fclose(fid);

            % Fit the 3rd-degree polynomial into the stem dielectric constants imaginary data
            poly = polyfit( indicesAll, stemDielImAll, 3);
            stemDielImAll_poly = poly(1)*indicesAll_poly.^3 + poly(2)*indicesAll_poly.^2 + ...
                                poly(3)*indicesAll_poly + poly(4);
            stemDielIm_poly = stemDielImAll_poly( find(indicesAll_poly == indices(1)) : find(indicesAll_poly == indices(end)), 1 );

            % Fit Gaussian Distribution and draw it
            pdNormalStemDielIm = fitdist( stemDielIm_poly, 'Normal' );

            end


            function [pdNormalLeafNumber] = fitdistLeafNumbers( obj, indices )

            %% NUMBER OF LEAVES
            % Take the number of leaves w.r.t. the date interval
            pathname = fileparts(mfilename('fullpath'));
            filename = strcat( pathname,'\corn_data\size\LeavesNumber.txt');
            fid = fopen(filename);
            leafNumbersAll = fscanf(fid,'%f');
            leafNumbers = leafNumbersAll(indices, 1);
            fclose(fid);

            % Fit Gaussian Distribution and draw it
            pdNormalLeafNumber = fitdist(leafNumbers,'Normal');

            end


            function [pdNormalLeafLength] = fitdistLeafLength(obj, indices, indicesAll, indicesAll_poly )

            %% LEAF LENGTH
            % Take the leaf lengths w.r.t. the date interval
            pathname = fileparts(mfilename('fullpath'));
            filename = strcat( pathname,'\corn_data\size\Length.txt');
            fid = fopen(filename);
            leafLengthsAll = fscanf(fid,'%f');
            fclose(fid);

            leafLengthsAll = leafLengthsAll * 0.01; % To convert from cm to meters

            % Fit the 3rd-degree polynomial into the stem height data
            poly = polyfit( indicesAll, leafLengthsAll, 3);
            leafLengthsAll_poly = poly(1)*indicesAll_poly.^3 + poly(2)*indicesAll_poly.^2 + ...
                                poly(3)*indicesAll_poly + poly(4);
            leafLengths_poly = leafLengthsAll_poly( find(indicesAll_poly == indices(1)) : find(indicesAll_poly == indices(end)), 1 );

            % Fit Gaussian Distribution and draw it
            pdNormalLeafLength = fitdist( leafLengths_poly, 'Normal' );

            end


            function [pdNormalLeafWidth] = fitdistLeafWidth(obj, indices, indicesAll, indicesAll_poly )

            %% LEAF WIDTH
            % Take the leaf widths w.r.t. the date interval
            pathname = fileparts(mfilename('fullpath'));
            filename = strcat( pathname,'\corn_data\size\Width.txt');
            fid = fopen(filename);
            leafleafWidthsAll = fscanf(fid,'%f');
            fclose(fid);

            leafleafWidthsAll = leafleafWidthsAll * 0.01; % To convert from cm to meters

            % Fit the 3rd-degree polynomial into the stem height data
            poly = polyfit( indicesAll, leafleafWidthsAll, 3);
            leafleafWidthsAll_poly = poly(1)*indicesAll_poly.^3 + poly(2)*indicesAll_poly.^2 + ...
                poly(3)*indicesAll_poly + poly(4);
            leafleafWidths_poly = leafleafWidthsAll_poly( find(indicesAll_poly == indices(1)) : find(indicesAll_poly == indices(end)), 1 );

            % Fit Gaussian Distribution and draw it
            pdNormalLeafWidth = fitdist( leafleafWidths_poly,'Normal');

            end


            function [pdNormalLeafThickness] = fitdistLeafThickness(obj, indices, indicesAll, indicesAll_poly )

            %% LEAF THICKNESS
            % Take the leaf lengths w.r.t. the date interval
            pathname = fileparts(mfilename('fullpath'));
            filename = strcat( pathname,'\corn_data\size\Thickness.txt');
            fid = fopen(filename);
            leafThicknessAll = fscanf(fid,'%f');
            fclose(fid);

            leafThicknessAll = leafThicknessAll * 0.001; % To convert from mm to meters

            % Fit the 3rd-degree polynomial into the stem height data
            poly = polyfit( indicesAll, leafThicknessAll, 3);
            leafThicknessAll_poly = poly(1)*indicesAll_poly.^3 + poly(2)*indicesAll_poly.^2 + ...
                                poly(3)*indicesAll_poly + poly(4);
            leafThickness_poly = leafThicknessAll_poly( find(indicesAll_poly == indices(1)) : find(indicesAll_poly == indices(end)), 1 );

            % Fit Gaussian Distribution and draw it
            pdNormalLeafThickness = fitdist( leafThickness_poly, 'Normal' );

            end


            function [pdNormalLeafBeginningAngle] = fitdistLeafBeginningAngle(obj, indices, indicesAll, indicesAll_poly )

            %% LEAF BEGIN ANGLE
            % Take the leaf beginning angles w.r.t. the date interval
            pathname = fileparts(mfilename('fullpath'));
            filename = strcat( pathname,'\corn_data\size\BegAngle.txt');
            fid = fopen(filename);
            leafBeginAnglesAll = fscanf(fid,'%f');
            fclose(fid);

            % Fit the 3rd-degree polynomial into the stem height data
            poly = polyfit( indicesAll, leafBeginAnglesAll, 3);
            leafBeginAnglesAll_poly = poly(1)*indicesAll_poly.^3 + poly(2)*indicesAll_poly.^2 + ...
                poly(3)*indicesAll_poly + poly(4);
            leafBeginAngles_poly = leafBeginAnglesAll_poly( find(indicesAll_poly == indices(1)) : find(indicesAll_poly == indices(end)), 1 );

            % Fit Gaussian Distribution and draw it
            pdNormalLeafBeginningAngle = fitdist( leafBeginAngles_poly,'Normal');

            end


            function [pdNormalLeafDielRe] = fitdistLeafDielRe( obj, indices, indicesAll, indicesAll_poly )

            %% LEAF DIELECTRIC CONSTANT REAL PART
            % Take the Leaf dielectric constant real part w.r.t. the date interval
            pathname = fileparts(mfilename('fullpath'));
            filename = strcat( pathname,'\corn_data\diel\Leaves_RE.txt');
            fid = fopen(filename);
            leafDielReAll = fscanf(fid,'%f');
            fclose(fid);

            % Fit the 3rd-degree polynomial into the leaf dielectric constants real data
            poly = polyfit( indicesAll, leafDielReAll, 3);
            leafDielReAll_poly = poly(1)*indicesAll_poly.^3 + poly(2)*indicesAll_poly.^2 + ...
                                poly(3)*indicesAll_poly + poly(4);
            leafDielRe_poly = leafDielReAll_poly( find(indicesAll_poly == indices(1)) : find(indicesAll_poly == indices(end)), 1 );

            % Fit Gaussian Distribution and draw it
            pdNormalLeafDielRe = fitdist( leafDielRe_poly, 'Normal' );

            end


            function [pdNormalLeafDielIm] = fitdistLeafDielIm( obj, indices, indicesAll, indicesAll_poly )

            %% LEAF DIELECTRIC CONSTANT IMAGINARY PART
            % Take the Leaf dielectric constant imaginary part w.r.t. the date interval
            pathname = fileparts(mfilename('fullpath'));
            filename = strcat( pathname,'\corn_data\diel\Leaves_IM.txt');
            fid = fopen(filename);
            leafDielImAll = fscanf(fid,'%f');
            fclose(fid);

            % Fit the 3rd-degree polynomial into the leaf dielectric constants imaginary data
            poly = polyfit( indicesAll, leafDielImAll, 3);
            leafDielImAll_poly = poly(1)*indicesAll_poly.^3 + poly(2)*indicesAll_poly.^2 + ...
                                poly(3)*indicesAll_poly + poly(4);
            leafDielIm_poly = leafDielImAll_poly( find(indicesAll_poly == indices(1)) : find(indicesAll_poly == indices(end)), 1 );

            % Fit Gaussian Distribution and draw it
            pdNormalLeafDielIm = fitdist( leafDielIm_poly, 'Normal' );

            end


            function [pdNormalCobNumber] = fitdistCobNumbers( obj, indices )

            %% NUMBER OF COBS
            % Take the number of cobs w.r.t. the date interval
            pathname = fileparts(mfilename('fullpath'));
            filename = strcat( pathname,'\corn_data\size\CobsNumber.txt');
            fid = fopen(filename);
            cobNumbersAll = fscanf(fid,'%f');
            cobNumbers = cobNumbersAll(indices, 1);
            fclose(fid);

            % Fit Gaussian Distribution and draw it
            pdNormalCobNumber = fitdist( cobNumbers, 'Normal' );

            end


            function [pdNormalCobLength] = fitdistCobLength(obj, indices, indicesAll, indicesAll_poly )

            %% COBS LENGTH
            % Take the Cob Lengths w.r.t. the date interval
            pathname = fileparts(mfilename('fullpath'));
            filename = strcat( pathname,'\corn_data\size\CobsLength.txt');
            fid = fopen(filename);
            cobLengthsAll = fscanf(fid,'%f');
            fclose(fid);

            cobLengthsAll = cobLengthsAll * 0.01; % To convert from cm to meters

            ind0 = find(cobLengthsAll == 0);
            indicesAllCob = ( max(ind0) : indicesAll(end) )';
            indicesAllCob_poly = ( min(indicesAllCob) : 0.1 : max( indicesAllCob ) )';
            cobLengthsAll = cobLengthsAll(indicesAllCob);

            % Fit the 3rd-degree polynomial into the cob length data
            poly = polyfit( indicesAllCob, cobLengthsAll, 3);
            cobLengthsAll_poly = poly(1)*indicesAllCob_poly.^3 + poly(2)*indicesAllCob_poly.^2 + ...
                                poly(3)*indicesAllCob_poly + poly(4);
            cobLengths_poly = cobLengthsAll_poly( find(indicesAllCob_poly == indices(1)) : find(indicesAllCob_poly == indices(end)), 1 );

            % Fit Gaussian Distribution and draw it
            if isempty( cobLengths_poly )
                cobLengths_poly = zeros(2,1);
            end
            pdNormalCobLength = fitdist( cobLengths_poly, 'Normal' );

            end


            function [pdNormalCobDiameter] = fitdistCobDiameter(obj, indices, indicesAll, indicesAll_poly )

            %% COBS DIAMETER
            % Take the Cob Diameters w.r.t. the date interval
            pathname = fileparts(mfilename('fullpath'));
            filename = strcat( pathname,'\corn_data\size\CobsDiameter.txt');
            fid = fopen(filename);
            cobDiametersAll = fscanf(fid,'%f');
            fclose(fid);

            cobDiametersAll = cobDiametersAll * 0.001; % To convert from mm to meters

            ind0 = find(cobDiametersAll == 0);
            indicesAllCob = ( max(ind0) : indicesAll(end) )';
            indicesAllCob_poly = ( min(indicesAllCob) : 0.1 : max( indicesAllCob ) )';
            cobDiametersAll = cobDiametersAll(indicesAllCob);

            % Fit the 3rd-degree polynomial into the cob diameter data
            poly = polyfit( indicesAllCob, cobDiametersAll, 3);
            cobDiametersAll_poly = poly(1)*indicesAllCob_poly.^3 + poly(2)*indicesAllCob_poly.^2 + ...
                                poly(3)*indicesAllCob_poly + poly(4);
            cobDiameters_poly = cobDiametersAll_poly( find(indicesAllCob_poly == indices(1)) : find(indicesAllCob_poly == indices(end)), 1 );

            % Fit Gaussian Distribution and draw it
            if isempty( cobDiameters_poly )
                cobDiameters_poly = zeros(2,1);
            end

            pdNormalCobDiameter = fitdist( cobDiameters_poly, 'Normal' );

            end


            function [pdNormalCobDielRe] = fitdistCobDielRe( obj, indices, indicesAll, indicesAll_poly )

            %% COB DIELECTRIC CONSTANT REAL PART
            % Take the Cob dielectric constant real part w.r.t. the date interval
            pathname = fileparts(mfilename('fullpath'));
            filename = strcat( pathname,'\corn_data\diel\Cobs_RE.txt');
            fid = fopen(filename);
            cobDielReAll = fscanf(fid,'%f');
            fclose(fid);

            ind0 = find(cobDielReAll == 0);
            indicesAllCob = ( max(ind0) : indicesAll(end) )';
            indicesAllCob_poly = ( min(indicesAllCob) : 0.1 : max( indicesAllCob ) )';
            cobDielReAll = cobDielReAll(indicesAllCob);

            % Fit the 3rd-degree polynomial into the cob dielectric constants real data
            poly = polyfit( indicesAllCob, cobDielReAll, 3);
            cobDielReAll_poly = poly(1)*indicesAllCob_poly.^3 + poly(2)*indicesAllCob_poly.^2 + ...
                                poly(3)*indicesAllCob_poly + poly(4);
            cobDielRe_poly = cobDielReAll_poly( find(indicesAllCob_poly == indices(1)) : find(indicesAllCob_poly == indices(end)), 1 );

            % Fit Gaussian Distribution and draw it
            if isempty( cobDielRe_poly )
                cobDielRe_poly = zeros(2,1);
            end

            % Fit Gaussian Distribution and draw it
            pdNormalCobDielRe = fitdist( cobDielRe_poly, 'Normal' );

            end


            function [pdNormalCobDielIm] = fitdistCobDielIm( obj, indices, indicesAll, indicesAll_poly )

            %% COB DIELECTRIC CONSTANT IMAGINARY PART
            % Take the Cob dielectric constant imaginary part w.r.t. the date interval
            pathname = fileparts(mfilename('fullpath'));
            filename = strcat( pathname,'\corn_data\diel\Cobs_IM.txt');
            fid = fopen(filename);
            cobDielImAll = fscanf(fid,'%f');
            fclose(fid);

            ind0 = find(cobDielImAll == 0);
            indicesAllCob = ( max(ind0) : indicesAll(end) )';
            indicesAllCob_poly = ( min(indicesAllCob) : 0.1 : max( indicesAllCob ) )';
            cobDielImAll = cobDielImAll(indicesAllCob);

            % Fit the 3rd-degree polynomial into the cob dielectric constants imaginary data
            poly = polyfit( indicesAllCob, cobDielImAll, 3);
            cobDielImAll_poly = poly(1)*indicesAllCob_poly.^3 + poly(2)*indicesAllCob_poly.^2 + ...
                                poly(3)*indicesAllCob_poly + poly(4);
            cobDielIm_poly = cobDielImAll_poly( find(indicesAllCob_poly == indices(1)) : find(indicesAllCob_poly == indices(end)), 1 );

            % Fit Gaussian Distribution and draw it
            if isempty( cobDielIm_poly )
                cobDielIm_poly = zeros(2,1);
            end

            % Fit Gaussian Distribution and draw it
            pdNormalCobDielIm = fitdist( cobDielIm_poly, 'Normal' );

            end
            

            function [stemHeight, stemBottomRadius, stemTopRadius, stemDielRe stemDielIm] = getRandomStemParams(obj)

            stemHeight = random(obj.pdfStemHeight);

            while stemHeight < 0.2 * obj.pdfStemHeight.mean() || stemHeight > 2 * obj.pdfStemHeight.mean()
                stemHeight = random(obj.pdfStemHeight);
            end


            stemBottomRadius = random(obj.pdfStemBottomDiameter);

            while stemBottomRadius < 0.2 * obj.pdfStemBottomDiameter.mean() || stemBottomRadius > 2 * obj.pdfStemBottomDiameter.mean()
                stemBottomRadius = random(obj.pdfStemBottomDiameter);
            end


            stemTopRadius = random(obj.pdfStemTopDiameter);

            while stemTopRadius < 0.2 * obj.pdfStemTopDiameter.mean() || stemTopRadius > 2 * obj.pdfStemTopDiameter.mean()
                stemTopRadius = random(obj.pdfStemTopDiameter);
            end


            stemTopRadius = random(obj.pdfStemTopDiameter);

            while stemTopRadius < 0.2 * obj.pdfStemTopDiameter.mean() || stemTopRadius > 2 * obj.pdfStemTopDiameter.mean()
                stemTopRadius = random(obj.pdfStemTopDiameter);
            end


            stemDielRe = random(obj.pdfStemDielRe);

            while stemDielRe < 0.2 * obj.pdfStemDielRe.mean() || stemDielRe > 2 * obj.pdfStemDielRe.mean()
                stemDielRe = random(obj.pdfStemDielRe);
            end


            stemDielIm = random(obj.pdfStemDielIm);

            while stemDielIm < 0.2 * obj.pdfStemDielIm.mean() || stemDielIm > 2 * obj.pdfStemDielIm.mean()
                stemDielIm = random(obj.pdfStemDielIm);
            end

            end


            function [leafNumbers] = getRandomLeafNumbers(obj)
                
            leafNumbers = random(obj.pdfLeafNumbers);

            while leafNumbers < 0.2 * obj.pdfLeafNumbers.mean() || leafNumbers > 2 * obj.pdfLeafNumbers.mean()
                leafNumbers = random(obj.pdfLeafNumbers);
            end

            leafNumbers = floor(leafNumbers + 0.5);

            end


            function [leafLength, leafWidth, leafThickness, leafBeginAngle, epsrLeafRe, epsrLeafIm] = getRandomLeafParams(obj)

            leafLength = random(obj.pdfLeafLength);

            while leafLength < 0.2 * obj.pdfLeafLength.mean() || leafLength > 2 * obj.pdfLeafLength.mean()
                leafLength = random(obj.pdfLeafLength);
            end


            leafWidth = random(obj.pdfLeafWidth);

            while leafWidth < 0.2 * obj.pdfLeafWidth.mean() || leafWidth > 2 * obj.pdfLeafWidth.mean()
                leafWidth = random(obj.pdfLeafWidth);
            end


            leafThickness = random(obj.pdfLeafThickness);

            while leafThickness < 0.2 * obj.pdfLeafThickness.mean() || leafThickness > 2 * obj.pdfLeafThickness.mean()
                leafThickness = random(obj.pdfLeafThickness);
            end


            leafBeginAngle = convertAngleTo360Range( random(obj.pdfLeafBeginningAngle) );

            epsrLeafRe = random(obj.pdfLeafDielRe);

            while epsrLeafRe < 0.2 * obj.pdfLeafDielRe.mean() || epsrLeafRe > 2 * obj.pdfLeafDielRe.mean()
                epsrLeafRe = random(obj.pdfLeafDielRe);
            end


            epsrLeafIm = random(obj.pdfLeafDielIm);

            while epsrLeafIm < 0.2 * obj.pdfLeafDielIm.mean() || epsrLeafIm > 2 * obj.pdfLeafDielIm.mean()
                epsrLeafIm = random(obj.pdfLeafDielIm);
            end

            end


            function [cobNumbers] = getRandomCobNumbers(obj)

            cobNumbers = random( obj.pdfCobNumbers );

            while cobNumbers < 0 || cobNumbers > obj.pdfCobNumbers.mean() + 2 * obj.pdfCobNumbers.sigma()
                cobNumbers = random(obj.pdfCobNumbers);
            end

            cobNumbers = floor( cobNumbers + 0.5);

            end


            function [cobLength, cobWidth, epsrCobRe, epsrCobIm] = getRandomCobParams(obj)

            cobLength = random(obj.pdfCobLength);

            while cobLength < 0.2 * obj.pdfCobLength.mean() || cobLength > 2 * obj.pdfCobLength.mean()
                cobLength = random(obj.pdfCobLength);
            end


            cobWidth = random(obj.pdfCobDiameter);

            while cobWidth < 0.2 * obj.pdfCobDiameter.mean() || cobWidth > 2 * obj.pdfCobDiameter.mean()
                cobWidth = random(obj.pdfCobDiameter);
            end


            epsrCobRe = random(obj.pdfCobDielRe);

            while epsrCobRe < 0.2 * obj.pdfCobDielRe.mean() || epsrCobRe > 2 * obj.pdfCobDielRe.mean()
                epsrCobRe = random(obj.pdfCobDielRe);
            end


            epsrCobIm = random(obj.pdfCobDielIm);

            while epsrCobIm < 0.2 * obj.pdfCobDielIm.mean() || epsrCobIm > 2 * obj.pdfCobDielIm.mean()
                epsrCobIm = random(obj.pdfCobDielIm);
            end

            end        
        
    end
    
    methods (Access = private)

        function [ posCOG ] = calculateCenterOfGravity( obj, position, ...
                particleDownAngle, particleAzimuthAngle, particleHalfLength )

        posCOG = position;

        particleDownAngle = degtorad( particleDownAngle );
        particleAzimuthAngle = degtorad( particleAzimuthAngle );

        
        posCOG(1,1) = posCOG(1,1) + particleHalfLength * sin( particleDownAngle ) * cos( particleAzimuthAngle ) ; % the x-value of center of gravity of the particle
        posCOG(2,1) = posCOG(2,1) + particleHalfLength * sin( particleDownAngle ) * sin( particleAzimuthAngle ); % the y-value  of center of gravity of the particle
        posCOG(3,1) = posCOG(3,1) + particleHalfLength * cos( particleDownAngle ); % the height of center of gravity of the particle

        end
    end
    
end

