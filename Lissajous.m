% Clear the workspace, command window
clear all; clc;

% Input parameters
x0 = 0;         % x initial position (10^-6 m)
y0 = 0;         % y initial position (10^-6 m)
A = 1;          % x peak to peak (10^-6 m)
B = 1;          % y peak to peak (10^-6 m)
fxMin = 1;      % minimum integer x-frequency to test (Hz)
fxMax = 10;     % maximum integer x-frequency to test (Hz)
fyMin = 1;      % minimum integer y-frequency to test (Hz)
fyMax = 10;     % maximum integer y-frequency to test (Hz)
dt = 0.001;     % time step size (s)
tMin = 0.000;   % start time (s)
tMax = 1.000;   % end time (s)
p = 6;          % number of decimal places relevant to calculation

% Generate data points for t, x, y, vx, vy, and VMag for each fx, fy pair
% using nested for-loops.
for fx = fxMin:fxMax
    for fy = fyMin:fyMax

        % Setup input vectors
        xi0 = [x0 y0];              % x,y initial position vector
        Ai = [A B];                 % x,y amplitudes vector
        fi = [fx fy];               % fx,fy frequencies vector
        tMax = lcm(fx,fy)/(fx*fy);  % calculate scan time
        t = tMin:dt:tMax;           % time step vector
        
        % Calculate data points for specified fx and fy.
        xi = Xi(xi0,Ai,fi,t);       % x,y-position points
        vi = Vi(Ai,fi,t);           % x,y-velocity component points
        v  = VMag(vi);              % speed
        
        % Calculate intersection points
        xInter = intersections(xi,p);

        % Plot scan and intersection points of scan curve.
        plot(xi(1,:),xi(2,:))
        xlabel(fx)
        ylabel(fy)
        hold on
        if sum(size(xInter))>0
            scatter(xInter(1,:),xInter(2,:))
        end
        pause;
        hold off
    end
end


% Function for rectangular coordinate position (xi, i = 1,2,...).
%
% xi-position (m)
% This function returns a matrix with the i-th position at each point 
% in time. The rows are the i-th position and the columns are different 
% points along the time interval. The inputs are i-th amplitude row vector 
% Ai (m), frequency fi (Hz) and row vector time t (s).
function xi = Xi(xi0,Ai,fi,t)
    
    % Check inputs for errors. If there are input errors, attempt
    % corrections and warn the user, or throw an error and stop. Inputs
    % xi0, Ai, and fi should all be row vectors of the same size and t is 
    % a row vector.
    %
    % Create flags that check for row vectors and matching vectors lengths.
    isRowFlag = isrow(xi0) && isrow(Ai) && isrow(fi) && isrow(t);
    lengthFlag = length(xi0)==length(Ai) && length(Ai)==length(fi);
    if isRowFlag && lengthFlag
        xi = xi0' + (Ai'/2).*sin(2*pi*fi'.*t);    % calculate positions
    else
        % Check if 'xi0' is a row vector. If not, assume column vector
        % input, transpose vector, and warn user.
        if ~isrow(xi0)
            xi0 = xi0';
            warning('Position function input xi0 should be a row vector');
        end
        
        % Check if 'Ai' is a row vector. If not, assume column vector
        % input, transpose vector, and warn user.
        if ~isrow(Ai)
            Ai = Ai';
            warning('Position function input Ai should be a row vector');
        end
        
        % Check if 'fi' is a row vector. If not, assume column vector
        % input, transpose vector, and warn user.
        if ~isrow(fi)
            fi = fi';
            warning('Position function input fi should be a row vector');
        end
        
        % Check if 't' is a row vector. If not, assume column vector
        % input, transpose vector, and warn user.
        if ~isrow(t)
            t = t';
            warning('Position function input t should be a row vector');
        end
        
        % Check if 'xi0', 'Ai', and 't' are now row vectors. If not, throw
        % error and stop, otherwise calculate positions.
        %
        % Update flags that check for row vectors and matching vectors
        % lengths.
        isRowFlag = isrow(xi0) && isrow(Ai) && isrow(fi) && isrow(t);
        lengthFlag = length(xi0)==length(Ai) && length(Ai)==length(fi);
        if isRowFlag && lengthFlag
            xi = xi0' + (Ai'/2).*sin(2*pi*fi'.*t); % calculate positions
        else
            message = "Error: Invalid input for position function. ";
            message = message + "Check the sizes of input vectors.";
            error(message);
        end
    end
end

% Functions for rectangular coordinate velocity (vi, i = 1,2,...) and speed.
%
% vi-velocity (m/s)
% This function returns a matrix with the i-th velocity at each point 
% in time. The rows are the i-th velocity and the columns are different 
% points along the time interval. The inputs are i-th amplitude row vector 
% Ai (m), frequency fi (Hz) and row vector time t (s).
function vi = Vi(Ai,fi,t)
    
    % Check inputs for errors. If there are input errors, attempt
    % corrections and warn the user, or throw an error and stop. Inputs
    % Ai and fi should both be row vectors of the same size and t is a row
    % vector.
    if  isrow(Ai) && isrow(fi) && isrow(t) && length(Ai)==length(fi)
        vi = pi*fi*Ai'.*cos(2*pi*fi'.*t); % calculate velocities
    else
        % Check if 'Ai' is a row vector. If not, assume column vector
        % input, transpose vector, and warn user.
        if ~isrow(Ai)
            Ai = Ai';
            warning('Velocity function input Ai should be a row vector');
        end
                
        % Check if 'fi' is a row vector. If not, assume column vector
        % input, transpose vector, and warn user.
        if ~isrow(fi)
            fi = fi';
            warning('Position function input fi should be a row vector');
        end
        
        % Check if 't' is a row vector. If not, assume column vector
        % input, transpose vector, and warn user.
        if ~isrow(t)
            t = t';
            warning('Velocity function input t should be a row vector');
        end
        
        % Check if 'Ai' and 't' are now row vectors and check to see if  
        % the size of 'Ai' and 'fi' are equal. If not, throw error
        % and stop, otherwise calculate positions.
        if isrow(Ai) && isrow(fi) && isrow(t) && length(Ai)==length(fi)
            vi = pi*fi*Ai'.*cos(2*pi*fi'.*t); % calculate velocities
        else
            message = "Error: Invalid input for velocity function. ";
            message = message + "Check the sizes of input vectors.";
            error(message);
        end
    end
end

% speed (m/s)
% This function returns a row vector with the speed at each point in time.
% Input v is a matrix where the rows are vx, vy, vz, etc. and the columns
% are different points along the time interval (e.g. v = [vx;vy]).
function vMag = VMag(v)
    vMag = sqrt(sum(v.^2));
end


% A function to find and return a matrix 'xInter' of (x,y) intersections 
% for the input curve. The first row represents the x's and the second row
% represents the y's. The columns are the i-th coordinate pair such that
% we have [x1 y1; x2 y2; ... ; xn yn]. The input is a 2xN matrix 'xi1' of 
% the data points that form the curve. In 'xi1' the first row also 
% represents the x's and the second row represents the y's. The columns 
% are the i-th coordinate pair of the data such that we have 
% [x1 y1; x2 y2; ... ; xn yn].
function xInter = intersections(xi1,p)
    
    xi1 = round(xi1,p);

    % Intialize an empty return matrix
    xInter =[];
    
    % Determine the size of the data set to calculate the number of line
    % segments and validate input.
    [r c] = size(xi1);
    
    % Validate input matrix
    if r ~=2 || c < 4
        message = 'Invalid input for the intersections() function.';
        message = message + 'Check size of input matrix.';
        error(message)
    end
    
    % Validate precision value.
    if ~isscalar(p)
        message = "Invalid input 'p' for the intersections() function.";
        message = message + "Input 'p' should be a scalar.";
        error(message)
    end
    
    % Check to see if the curve has overlapping sections.  If there is
    % overlap, select a subset that is representative of the entire curve,
    % but does not have overlap. Perform the intersection analysis on this
    % subset.
    %
    % Store the indices where the curve starts to double over itself.
    turn_indices = [];
    for i=2:(c-1)
        % Check to see if the point before and after each step is equal to
        % see if the curve overlaps.
        if ~round(xi1(1,i-1) - xi1(1,i+1), p) && ~round(xi1(2,i-1) - xi1(2,i+1), p)
            turn_indices = [turn_indices, i];
        end
        
    end
    
    % If the curve overlaps, select a representative subset and reset the
    % size of the data set.
    if length(turn_indices) > 1
        xi1 = xi1(:, turn_indices(1):turn_indices(2));
        [r c] = size(xi1);
    end
   
    % Copy input matrix and shift first data point the last. The curves are
    % supposed to be periodic, so the last line segment connects to the
    % first. This new matrix contains the line segment endpoints.
    xi2 = [xi1(:,2:c),xi1(:,1)];
    
    % Simulatenously calculate all of the line segment slopes.
    dx = xi1(1,:) - xi2(1,:);
    dy = xi1(2,:) - xi2(2,:);
    m = round(dy./dx, p);
    
    % Simulateneously calculate all of the line segment y-intercept points.
    b = round(xi1(2,:) - m.*xi1(1,:), p);
    
    % Use a double for-loop to check if line segments intersect. The first
    % for-loop runs through every line segment.
    for i = 1:(c-2)
        
        % Determine the upper and lower x-boundaries for intersection from
        % the 1st line segment.
        if xi1(1,i) >= xi1(1,i+1)
            x_upper1 = xi1(1,i);
            x_lower1 = xi1(1,i+1);
        elseif xi1(1,i) < xi1(1,i+1)
            x_upper1 = xi1(1,i+1);
            x_lower1 = xi1(1,i);
        end
        
        % Determine the upper and lower y-boundaries for intersection from
        % the 1st line segment.
        if xi1(2,i) >= xi1(2,i+1)
            y_upper1 = xi1(2,i);
            y_lower1 = xi1(2,i+1);
        elseif xi1(2,i) < xi1(2,i+1)
            y_upper1 = xi1(2,i+1);
            y_lower1 = xi1(2,i);
        end
        
        % The second for-loop, which runs through every line segment not
        % previously checked by the first for-loop.  Note that the starting
        % point for 'j' is 'i+2' because the the 'i' and 'i+1' line
        % segments will aways intersect at their shared end point. This is
        % not an intersection of the curve that the line segments
        % represent.
        for j = i+2:(c-1)
            
            % Determine the upper and lower x-boundaries for intersection
            % from the 2nd line segment.
            if xi1(1,j) >= xi1(1,j+1)
                x_upper2 = xi1(1,j);
                x_lower2 = xi1(1,j+1);
            elseif xi1(1,j) < xi1(1,j+1)
                x_upper2 = xi1(1,j+1);
                x_lower2 = xi1(1,j);
            end

            % Determine the upper and lower y-boundaries for intersection
            % from the 2nd line segment.
            if xi1(2,j) >= xi1(2,j+1)
                y_upper2 = xi1(2,j);
                y_lower2 = xi1(2,j+1);
            elseif xi1(2,j) < xi1(2,j+1)
                y_upper2 = xi1(2,j+1);
                y_lower2 = xi1(2,j);
            end
            
            % Check to see if the lines are parallel, if not find their
            % intersection point.  If the difference between the slopes of
            % the lines is zero or both slopes are infinite, then they're 
            % parallel.
            delta_m = round(m(j) - m(i),p); % Difference of line slopes
            
            % Scenario 1:
            % Check to see if both line segments are not verticle or
            % parallel to each other.  If not, determine the intercept.
            if abs(m(i)) ~= Inf && abs(m(j)) ~= Inf && delta_m ~= 0
                x_inter = -(b(j) - b(i))/delta_m;   % x-intercept point
                y_inter = m(i)*x_inter + b(i);      % y-intercept point
            %
            % Scenario 2:
            % Check to see if second line segment is verticle and the first
            % line is not verticle. If so, determine the intercept point.
            elseif abs(m(i)) ~= Inf && abs(m(j)) == Inf
                x_inter = x_upper2;                 % x-intercept point
                y_inter = m(i)*x_inter + b(i);      % y-intercept point
            %
            % Scenario 3:
            % Check to see if first line segment is verticle and the second
            % line is not verticle. If so, determine the intercept point.
            elseif abs(m(i)) == Inf && abs(m(j)) ~= Inf
                x_inter = x_upper1;                 % x-intercept point
                y_inter = m(j)*x_inter + b(j);      % y-intercept point
            %
            % Scenario 4:
            % The lines are parallel and do not intersect. Continue on to
            % next set of line segments.
            else
                continue;
            end
            
            % Round off any excess digits for the x-y intercept.
            x_inter = round(x_inter,p);
            y_inter = round(y_inter,p);
            
            % Check if the x,y-intersection point is inside both of 
            % the line segment intervals. If it is, add it to the 
            % 'xInter' matrix.
            x_condition1 = (x_inter <= x_upper1) && (x_inter >= x_lower1);
            y_condition1 = (y_inter <= y_upper1) && (y_inter >= y_lower1);
            x_condition2 = (x_inter <= x_upper2) && (x_inter >= x_lower2);
            y_condition2 = (y_inter <= y_upper2) && (y_inter >= y_lower2);
            if x_condition1 && y_condition1 && x_condition2 && y_condition2
                xInter = [xInter, [x_inter; y_inter]];
            end
        end
    end
end
