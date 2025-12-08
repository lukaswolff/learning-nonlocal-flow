classdef Macro_model < handle
    % MACRO_MODEL
    % Each instance of this class models the "macro model" on a given
    % geometry. Parallel computation and comparison of different sets of
    % parameters or models is possible.
    % To come:
    % -Input of different geometries
    % -Tools for parameter  estimation



    properties
        %geometry
        alpha   % Angle of obstacle
        a  % coordinate where obstacle starts

        % Conveyor belt
        dt    % Step size: Shutter time of camera is integer multiple
        dx   % space step-size x
        dy    % space step-size y 
        lambda_x    % gridconstant
        lambda_y   % gridconstant
        T       % end of time horizon
        sup1      % region for convolution

        Nx                % Number of discrete points in x
        Ny               % Number of discrete points in y
        NT                % Number of discrete points in t

        vx
        vy
        Vx
        Vy

        X
        Y

        x
        y
        t
        t_data
        tk_data
        t_sol_init
        sparsityfactor_time_data
        relation_time_data
        
        obstacle
        input_obstacle

        % other
        rho_sol   % Solution Matrix of the Density
        rho_init
        rho_net
        results
        velocities
        velocities_NN
        NNinput
        % Parameters
        eps
        sigma_1
        gamma

        % Data
        MD_centers
        X_MD
        v_T
        rho_data_60
        jsonfile
        datafile
        times
        data
        cargosize
        conveyor_z_coord

    end

    methods
        function self = Macro_model(alpha,options)
            %Construct an instance for a given setup of the conveyor belt.
            arguments
                alpha = 45 %angle of regulator
                options.obstacle = {} %Obstacle position and size {x,y,r}, Default: no obstacle
                options.parameters = [1.5,2000,200] %Numerical parameters eps, sigma, gamma
            end

            % Geometry
            self.alpha = alpha;
            self.a = 0.51;
            if alpha < 55
                self.a = 0.38;
            end


            self.input_obstacle = options.obstacle;


            load('default_setup.mat') %#ok<LOAD>
            self.MD_centers = MD_centers;
            self.X_MD = X_MD;
            % load('rho_data_60.mat') %#ok<LOAD>
            % self.rho_data_60 = rho_data_60;


            self.initialize_conveyor_belt();
            %self.initialize_singularizer()

            % Parameters

            self.eps = options.parameters(1); %2.1
            self.sigma_1 = options.parameters(2);%10000;
            self.gamma = options.parameters(3);%500; %Coefficient for Heaviside Approximation


            self.fetch_initial_data();

        end

        function recompute(self)
            % reset model and computation
            self = Macro_model(self.alpha);
            self.compute();
            self.plot();

        end

        function initialize_conveyor_belt(self)
            % Initialize conveyor belt
            % Independent of eps, sigma_1, gamma !!!

            %barrier_factor_singularizer = 1;            % controls the strength of the vector field behind the singularizer
            barrier_factor_boundary = 1;                % controls the strength of the vector field outside of the conveyor belt
            num_bordercells = 1; %5

            self.dt = 0.02;%0.0156 / 3;           % Step size: Shutter time of camera is integer multiple
            self.dx = 0.01;                 % space step-size x
            self.dy = 0.01;                 % space step-size y
            self.lambda_x = self.dt/self.dx;           % gridconstant
            self.lambda_y = self.dt/self.dy;           % gridconstant
            self.T = 12;                      % end of time horizon
            self.sup1 = 0.1 ;                % region for convolution
            self.x = 0:self.dx:0.7;              % discrete space vector  // width of the conveyorbelt
            self.y = 0:self.dy:0.8;                 % discrete space vector  // length of the conveyorbelt
            self.t = 0:self.dt:self.T;                    % discrete time vector
            self.v_T  = 0.1;%0.37;%0.335;                 % Conveyorbelt velocity
            v_T = self.v_T; %#ok<*PROP>

            self.Nx = length(self.x);                     % Number of discrete points in x
            self.Ny = length(self.y);                     % Number of discrete points in y
            self.NT = length(self.t);                     % Number of discrete points in t


            %rho_s = zeros(self.Nx,self.Ny);               % Temporary Solution Matrix of the Density // Just used for the Dimension Splitting
            self.obstacle = zeros(self.Nx,self.Ny);

            % Basic Velocity field of the conveyor belt
            self.vx = 0*ones(self.Nx,self.Ny);
            self.vy = v_T*ones(self.Nx,self.Ny);
            self.Vx = 0*ones(self.Nx,self.Ny);
            self.Vy = v_T*ones(self.Nx,self.Ny);

            [X, Y] = meshgrid(self.y,self.x); % change of direction in coordinates x-->
            self.X = X;self.Y = Y;


            % Compute region C
            ox = [self.a,0];
            h1 = 0.44;  % Length of singularizer (Hypotenuse)
            h2 = 0.04;%0.040; % Thickness of singularizer (Region C)
            h3 = 0.5 ;  % Size of the region B1

            alpha = self.alpha;
            for i=1:self.Nx
                for j=1:self.Ny
                    if dot([cos(pi*alpha/180) sin(pi*alpha/180)],ox+h1*[cos(pi*alpha/180) sin(pi*alpha/180)]-[X(i,j) Y(i,j)])>0 && dot([cos(pi*alpha/180) sin(pi*alpha/180)],ox-h1*[cos(pi*alpha/180) sin(pi*alpha/180)]-[X(i,j) Y(i,j)])<0      &&     dot([sin(pi*alpha/180) -cos(pi*alpha/180)],ox+h2*[sin(pi*alpha/180) -cos(pi*alpha/180)]-[X(i,j) Y(i,j)])>0 && dot([sin(pi*alpha/180) -cos(pi*alpha/180)],ox-h2*[sin(pi*alpha/180) -cos(pi*alpha/180)]-[X(i,j) Y(i,j)])<0
                        self.Vx(i,j) =  0;%v_T*sin(pi*alpha/180)*cos(pi*alpha/180);
                        self.Vy(i,j) =  0;%v_T*cos(pi*alpha/180)*cos(pi*alpha/180);
                        self.obstacle(i,j) = 1;
                    end
                end
            end


            % Compute region B1
            ox = [self.a,0] + h3*[sin(pi*alpha/180) -cos(pi*alpha/180)];

            for i=1:self.Nx
                for j=1:self.Ny
                    if dot([cos(pi*alpha/180) sin(pi*alpha/180)],ox+h1*[cos(pi*alpha/180) sin(pi*alpha/180)]-[X(i,j) Y(i,j)])>0 && dot([cos(pi*alpha/180) sin(pi*alpha/180)],ox-h1*[cos(pi*alpha/180) sin(pi*alpha/180)]-[X(i,j) Y(i,j)])<0      &&     dot([sin(pi*alpha/180) -cos(pi*alpha/180)],ox+h3*[sin(pi*alpha/180) -cos(pi*alpha/180)]-[X(i,j) Y(i,j)])>0 && dot([sin(pi*alpha/180) -cos(pi*alpha/180)],ox-h3*[sin(pi*alpha/180) -cos(pi*alpha/180)]-[X(i,j) Y(i,j)])<0
                        self.Vx(i,j) =  0;%barrier_factor_singularizer* v_T*cos(pi*alpha/180);
                        self.Vy(i,j) = 0;%- barrier_factor_singularizer* v_T*sin(pi*alpha/180);
                        self.obstacle(i,j) = 1;
                    end
                end
            end

            %additional obstacle
            if length(self.input_obstacle) > 0
                center = [self.input_obstacle{2} self.input_obstacle{1}];
                %radius = 0.05;
                for i=1:self.Nx
                    for j=1:self.Ny
                        if norm([X(i,j) Y(i,j)]-center)<self.input_obstacle{3}
                            self.Vx(i,j) =  0;%barrier_factor_singularizer* v_T*cos(pi*alpha/180);
                            self.Vy(i,j) = 0;%- barrier_factor_singularizer* v_T*sin(pi*alpha/180);
                            self.obstacle(i,j) = 1;
                        end
                    end
                end
            end


            % Define Vx, Vy at the border of the conveyor belt

            self.obstacle(1:num_bordercells, :) = 1;
            self.obstacle(end-num_bordercells+1:end, :) = 1;

            self.Vx(1:num_bordercells, :) = barrier_factor_boundary*v_T;
            self.Vy(1:num_bordercells, :) = 0;
            %
            %             % Smoothing of the velocity field by a gaussian blur
            %             PSF = fspecial('gaussian',1/self.dx/25,1/self.dx/50);
            %             %PSF = fspecial('gaussian',1/self.dx/5,1/self.dx/10);
            %             self.Vx = imfilter(self.Vx,PSF,'symmetric','conv');
            %             self.Vy = imfilter(self.Vy,PSF,'symmetric','conv');
        end

        function load_setup(self,jsonfile)
            %Input jsonfile: 0 for manual selection, 1 to load previous
            %data, filename for explicit filename

            if jsonfile == 0
                [file,path] = uigetfile('*.json', 'Select Multiple Files', 'MultiSelect', 'off' );
                jsonfile = fullfile(path,file);
            end

            if jsonfile == 1 %load previous data
                jsonfile = self.jsonfile;
            end
            
            %if jsonfile input is a filename:
            self.jsonfile = jsonfile;

            jsondata = jsondecode(fileread(jsonfile));
            if jsondata.conveyor.type ~= 'linear'
                error('Unknown Conveyor Type')
            end

            self.v_T  = jsondata.conveyor.velocity(2);                 % Conveyorbelt velocity
            v_T = self.v_T; %#ok<*PROP>

            self.dt = 0.002/v_T; %0.005;       %0.015    % Step size: Shutter time of camera is integer multiple
            self.dx = 0.01;                 % space step-size x
            self.dy = 0.01;                 % space step-size y
            self.lambda_x = self.dt/self.dx;           % gridconstant
            self.lambda_y = self.dt/self.dy;           % gridconstant
            self.T = 25;                      % end of time horizon
            self.x = 0:self.dx:jsondata.conveyor.size(1) ;              % discrete space vector  // width of the conveyorbelt
            self.y = 0:self.dy:jsondata.conveyor.size(2);                 % discrete space vector  // length of the conveyorbelt
            self.t = 0:self.dt:self.T;                    % discrete time vector
            
            X_length = jsondata.conveyor.size(2);
            Y_length = jsondata.conveyor.size(1);
            self.conveyor_z_coord = jsondata.conveyor.size(3) + 0.5*jsondata.cargo.size(2);
            
            self.cargosize = jsondata.cargo.size(1);

            self.Nx = length(self.x);                     % Number of discrete points in x
            self.Ny = length(self.y);                     % Number of discrete points in y
            self.NT = length(self.t);                     % Number of discrete points in t

            self.obstacle = zeros(self.Nx,self.Ny);

            % Basic Velocity field of the conveyor belt
            self.vx = 0*ones(self.Nx,self.Ny);
            self.vy = v_T*ones(self.Nx,self.Ny);
            self.Vx = 0*ones(self.Nx,self.Ny);
            self.Vy = v_T*ones(self.Nx,self.Ny);

            [X, Y] = meshgrid(self.y,self.x); % change of direction in coordinates x-->
            self.X = X;self.Y = Y;
            regulators = {};
            %regulatorstrings = "regulator" + (1:3);

            if isfield(jsondata,'obstacles')
                obstacles = jsondata.obstacles;
                for a = 1:length(obstacles)
                    if strcmp(obstacles(a).type,'regulatorRectangle')
                        regulator = obstacles(a);
                        regulator.type = "rectangle";
                        regulators{1} = regulator;
                    elseif obstacles(a).type == 'cylinder'
                        obstacle = obstacles(a);
                        center = [obstacle.position(2),obstacle.position(1)];
                        for i=1:self.Nx
                            for j=1:self.Ny
                                if norm([X(i,j) Y(i,j)]-center)<obstacle.size(1)
                                    self.Vx(i,j) =  0;%barrier_factor_singularizer* v_T*cos(pi*alpha/180);
                                    self.Vy(i,j) = 0;%- barrier_factor_singularizer* v_T*sin(pi*alpha/180);
                                    self.obstacle(i,j) = 1;
                                end
                            end
                        end

                    else
                        error('Unknown Conveyor Type')
                    end
                end

            end

            if isfield(jsondata,'regulator1')
                regulators{1} = jsondata.regulator1;
            end
            if isfield(jsondata,'regulator2')
                regulators{2} = jsondata.regulator2;
            end
            if isfield(jsondata,'regulator3')
                regulators{3} = jsondata.regulator3;
            end

            
            %regulators = {jsondata.regulator1,jsondata.regulator2,jsondata.regulator3};
            for k = 1:length(regulators)
                regulator = regulators{k};
                if regulator.type == "rectangle"
                    ox = [regulator.position(2),regulator.position(1)];%-[0.01,0]; %move one pixel to the left
                    h1 = regulator.size(1); % Length of singularizer (Hypotenuse)
                    h2 = regulator.size(2); % Thickness of singularizer (thickened (why?))
                    alpha = pi/180 * regulator.position(6);
                    for i=1:self.Nx
                        for j=1:self.Ny
                            if dot([cos(alpha) sin(alpha)],ox+h1*[cos(alpha) sin(alpha)]-[X(i,j) Y(i,j)])>0 && dot([cos(alpha) sin(alpha)],ox-h1*[cos(alpha) sin(alpha)]-[X(i,j) Y(i,j)])<0      &&     dot([sin(alpha) -cos(alpha)],ox+h2*[sin(alpha) -cos(alpha)]-[X(i,j) Y(i,j)])>0 && dot([sin(alpha) -cos(alpha)],ox-[X(i,j) Y(i,j)])<0
                                self.Vx(i,j) =  0;%v_T*sin(pi*alpha/180)*cos(pi*alpha/180);
                                self.Vy(i,j) =  0;%v_T*cos(pi*alpha/180)*cos(pi*alpha/180);
                                self.obstacle(i,j) = 1;
                            end
                        end
                    end
                elseif regulator.type(1:7) == "polynom"
                    %Kommentar Annika: Die Kurven entsprechen der Beschreibung c*a^b. In der Konfigurationsdatei ist dies beim Abweisertyp als „polynom_cA(xb)“ hinterlegt, bspw. entspricht „polynom_0.8AAA“ einer Kurve mit 0.8*a^3.  Die Kurven haben alle die Länge -0.5 < a < 0.5. Die Position in der Konfigurationsdatei beschreibt den Startpunkt der Kurve an der linken Außenkante. 
                    c = str2double(regulator.type(9:11));
                    A = regulator.type(12:end);
                    b = count(A,"A");
                    if b ~= length(A)
                        error('regulator type doesnt fit the format: non A character in AAA')
                    end
                    
                    width = 0.08;
                    x_anchor = regulator.position;
                    for i=1:self.Nx
                        for j=1:self.Ny
                            if abs(X(i,j) - (x_anchor + 0.5)) <= 0.5
                                a = X(i,j) - (x_anchor + 0.5);
                                y_polynomial = -c*(-0.5)^(b) + c*a^(b);
                                if Y(i,j) < (y_polynomial) && Y(i,j) > (y_polynomial-width)
                                    self.Vx(i,j) =  0;%v_T*sin(pi*alpha/180)*cos(pi*alpha/180);
                                    self.Vy(i,j) =  0;%v_T*cos(pi*alpha/180)*cos(pi*alpha/180);
                                    self.obstacle(i,j) = 1;
                                end
                            end
                        end
                    end




                end
              
            end

            
            if isfield(jsondata,'obstacle')
                obstacle = jsondata.obstacle;
                if obstacle.type == 'cylinder'
                    center = [obstacle.position(2),obstacle.position(1)];
                    for i=1:self.Nx
                        for j=1:self.Ny
                            if norm([X(i,j) Y(i,j)]-center)<obstacle.size(1)
                                self.Vx(i,j) =  0;%barrier_factor_singularizer* v_T*cos(pi*alpha/180);
                                self.Vy(i,j) = 0;%- barrier_factor_singularizer* v_T*sin(pi*alpha/180);
                                self.obstacle(i,j) = 1;
                            end
                        end
                    end
                end
            end
                    


            % % Compute region B1
            % ox = [self.a,0] + h3*[sin(pi*alpha/180) -cos(pi*alpha/180)];
            %
            % for i=1:self.Nx
            %     for j=1:self.Ny
            %         if dot([cos(pi*alpha/180) sin(pi*alpha/180)],ox+h1*[cos(pi*alpha/180) sin(pi*alpha/180)]-[X(i,j) Y(i,j)])>0 && dot([cos(pi*alpha/180) sin(pi*alpha/180)],ox-h1*[cos(pi*alpha/180) sin(pi*alpha/180)]-[X(i,j) Y(i,j)])<0      &&     dot([sin(pi*alpha/180) -cos(pi*alpha/180)],ox+h3*[sin(pi*alpha/180) -cos(pi*alpha/180)]-[X(i,j) Y(i,j)])>0 && dot([sin(pi*alpha/180) -cos(pi*alpha/180)],ox-h3*[sin(pi*alpha/180) -cos(pi*alpha/180)]-[X(i,j) Y(i,j)])<0
            %             self.Vx(i,j) =  0;%barrier_factor_singularizer* v_T*cos(pi*alpha/180);
            %             self.Vy(i,j) = 0;%- barrier_factor_singularizer* v_T*sin(pi*alpha/180);
            %             self.obstacle(i,j) = 1;
            %         end
            %     end
            % end
            %
            % %additional obstacle
            % if length(self.input_obstacle) > 0
            %     center = [self.input_obstacle{2} self.input_obstacle{1}];
            %     %radius = 0.05;
            %     for i=1:self.Nx
            %         for j=1:self.Ny
            %             if norm([X(i,j) Y(i,j)]-center)<self.input_obstacle{3}
            %                 self.Vx(i,j) =  0;%barrier_factor_singularizer* v_T*cos(pi*alpha/180);
            %                 self.Vy(i,j) = 0;%- barrier_factor_singularizer* v_T*sin(pi*alpha/180);
            %                 self.obstacle(i,j) = 1;
            %             end
            %         end
            %     end
            % end


            % Define Vx, Vy at the border of the conveyor belt

            num_bordercells = 1;
            barrier_factor_boundary = 1;

            self.obstacle(1:num_bordercells, :) = 1;
            self.obstacle(end-num_bordercells+1:end, :) = 1;

            self.Vx(1:num_bordercells, :) = barrier_factor_boundary*v_T;
            self.Vy(1:num_bordercells, :) = 0;

            %self.fetch_initial_data();
            %
            %             % Smoothing of the velocity field by a gaussian blur
            %             PSF = fspecial('gaussian',1/self.dx/25,1/self.dx/50);
            %             %PSF = fspecial('gaussian',1/self.dx/5,1/self.dx/10);
            %             self.Vx = imfilter(self.Vx,PSF,'symmetric','conv');
            %             self.Vy = imfilter(self.Vy,PSF,'symmetric','conv');

        end

        function load_data(self,datafile)

            % Takes input Data with n = 6 datapoints per piece (x,y,z,vx,vy,vz)
            %n = 3; %number of datapoints per entry (e.g. 3:(x,y,z); 6:(x,y,z,vx,vy,vz))

            if datafile == 0
                [file,path] = uigetfile('*.csv', 'Select Multiple Files', 'MultiSelect', 'off' );
                datafile = fullfile(path,file);
            end

            if datafile == 1 %load previous data
                datafile = self.datafile;
            end

            self.datafile = datafile;

            % Set up the Import Options and import the data
            opts = delimitedTextImportOptions("NumVariables", 193);

            % Specify range and delimiter
            opts.DataLines = [1, Inf];
            opts2 = opts;
            opts.Delimiter = [";","[","]",","];
            opts.MissingRule = 'omitvar';

            opts.ConsecutiveDelimitersRule = "join";
            temp(1,1:193) = "double";
            opts.VariableTypes = temp;


            % Import the data
            data = readmatrix(datafile,opts);
            opts2.Delimiter = [";"];
            temp1 = readtable(datafile, opts2);
            temp1 = table2cell(temp1);
            n = length(str2num(temp1{1,2}));



            % Convert to output type
            times = data(:,1);
            if isnan(data(1,end)) %Last line sometimes is NaN, change to end-1
                data = data(:,2:end-1); 
            else 
                data = data(:,2:end);
            end
            sz = size(data);
            
            data = permute(reshape(data',n,sz(2)/n,sz(1)),[3,2,1]); %data: (timestep,cargo,position)

            if n == 3
                for k = 1:size(data,1)-1
                    data(k,:,4:6) = (data(k+1,:,1:3) - data(k,:,1:3))/(times(k+1)-times(k));
                end
            end

            self.times = times;
            dt_data = times(2);
            T_end = times(end);
            if dt_data ~= 0.005
                error('Check timestep size')
            end
            self.data = data;

            rho_max = 2004;
            N = size(data,2);
            NT = size(data,1);

            %Smooth data
            sig0 = 0.012;            % parameter initial density smoothed
            %sig0 = 0.01;
            r = 5; %sparsity factor: only choose every r'th timestep
            K = floor(T_end/(r*self.dt));
            rho = single(zeros(self.Nx,self.Ny,K));              % Solution Matrix of the Density
            X = self.X;
            Y = self.Y;
            dx = self.dx;
            dy = self.dy;

            self.conveyor_z_coord = data(1,1,3);


            self.v_T = data(4,3,5); %read velocity of particle 3 at k = 4 

            % agr = size(self.MD_centers);
            % N = agr(1);   % N = 192;
            %figure;
            %hold on

            t_data = [];
            f = r*self.dt / dt_data;
            %Explanation of timescales: evaluate data at every dt*r seconds.  
            for k=1:K  %can use parfor
                t_data(k) = r*(k-1)*self.dt; %time in seconds at k'th position in rho_init
                tk_data(k) = floor(1+f*(k-1)); %position in data() correspondent to this time  
                t_sol_init(k) = 1 + r*(k-1);%position in rho_sol correspondent to this time
             
                %tic
                for i = 1:N

                    if abs(data(tk_data(k),i,3)-self.conveyor_z_coord) > 0.01
                        continue %ignore particles that fell off the end of the conveyor
                    end

                    part_x = data(tk_data(k),i,2);
                    part_y = data(tk_data(k),i,1);

                    % Trapez Rule for Integration

                    % rho(:,:,k) = rho(:,:,k) + 0.125*1/(rho_max*2*pi*sig0^2)*exp(-1/(2*sig0^2)*((X-part_x-0.0*dx).^2+(Y-part_y-0.5*dy).^2));
                    % rho(:,:,k) = rho(:,:,k) + 0.125*1/(rho_max*2*pi*sig0^2)*exp(-1/(2*sig0^2)*((X-part_x-0.0*dx).^2+(Y-part_y+0.5*dy).^2));
                    % rho(:,:,k) = rho(:,:,k) + 0.125*1/(rho_max*2*pi*sig0^2)*exp(-1/(2*sig0^2)*((X-part_x+0.5*dx).^2+(Y-part_y-0.0*dy).^2));
                    % rho(:,:,k) = rho(:,:,k) + 0.125*1/(rho_max*2*pi*sig0^2)*exp(-1/(2*sig0^2)*((X-part_x-0.5*dx).^2+(Y-part_y+0.0*dy).^2));
                    % rho(:,:,k) = rho(:,:,k) + 0.5*1/(rho_max*2*pi*sig0^2)*exp(-1/(2*sig0^2)*((X-part_x).^2+(Y-part_y).^2));

                    rho(:,:,k) = rho(:,:,k) + 1/(rho_max*2*pi*sig0^2)*exp(-1/(2*sig0^2)*((X-part_x).^2+(Y-part_y).^2));
                    %radius = 0.012;
                    %rho(:,:,k) = rho(:,:,k) + sqrt(6/pi)*radius^2*1/sig0*exp(-1/(2*sig0^2)*((X-part_x).^2+(Y-part_y).^2));

                end
                k;
            end
            self.t_data = t_data;
            self.tk_data = tk_data;
            self.rho_init = rho;
            self.t_sol_init = t_sol_init;
            self.relation_time_data = f;
            self.sparsityfactor_time_data = r;

        end

        function samples = extract_fluxdata(self)
            %rho = self.rho_init;
            Nk = size(self.rho_init,3);
            Ni = size(self.data,2);
            samples = {};
            for k = 1:Nk-1
                for i = 1:Ni

                    %Velocity
                    f = self.relation_time_data;
                    velocity = squeeze(self.data(round(1+f*(k-1)),i,4:5));
                    velocity_diff = round(velocity - [0;self.v_T],2);
                    if self.data(round(1+f*(k-1)),i,2) > 2 % || velocity_diff == [0;0]  instead filter for rho > 0.9 (below)
                        continue
                    end
                    if abs(self.data(round(1+f*(k-1)),i,3)-self.conveyor_z_coord) > 0.02 % pieces that fell off
                        continue
                    end
                    

                    %Density
                    size_x = size(self.rho_init,2);
                    size_y = size(self.rho_init,1);
                    x_particle = self.data(round(1+f*(k-1)),i,2);
                    y_particle = self.data(round(1+f*(k-1)),i,1);
                    x_index = round(x_particle/self.dx); 
                    y_index = round(y_particle/self.dy);

                    assert(y_index > 0, 'Positive Indizes');
                    

                    if self.rho_init(y_index,x_index,k) < 0.7 
                        continue
                    end

                    % Relevanter Ausschnitt der Daten

                    scopesize_x = 10;
                    scopesize_y = 10;

                    offset_x_l = max(1-(x_index - scopesize_x),0);
                    offset_x_r = max((x_index + scopesize_x)-size_x,0);

                    offset_y_l = max(1-(y_index - scopesize_y),0);
                    offset_y_r = max((y_index + scopesize_y)-size_y,0);

                    scope_x = x_index - scopesize_x + offset_x_l : x_index + scopesize_x - offset_x_r; %Maybe add Factor 2 to look ahead
                    scope_y = y_index - scopesize_y + offset_y_l : y_index + scopesize_y - offset_y_r;

                    sample.density = single(zeros(2*scopesize_y+1,2*scopesize_x+1));
                    sample.density(1+offset_y_l : 2*scopesize_y+1 - offset_y_r,1+offset_x_l : 2*scopesize_x+1 - offset_x_r) = self.rho_init(scope_y,scope_x,k);

                    sample.obstacle = boolean(ones(2*scopesize_y+1,2*scopesize_x+1));
                    sample.obstacle(1+offset_y_l : 2*scopesize_y+1 - offset_y_r,1+offset_x_l : 2*scopesize_x+1 - offset_x_r) = self.obstacle(scope_y,scope_x);

                    sample.velocity = squeeze(self.data(round(1+f*(k-1)),i,4:5)) - [0;self.v_T];
                    sample.velocity = sample.velocity/self.v_T;
                    sample.time = k;
                    samples{end+1} = sample;
                end
                k;
            end
            %samples = 0;
            % figure(100)
            % imagesc(density)

        end

        function plot_obstacle(self)
            %%% Plot
            figure(5);
            imagesc(self.y,self.x,-self.obstacle);
            colormap gray
            hold on
            title('Conveyor Belt Setup')
            ylabel('Width in m')
            xlabel('Direction of flow in m')
            %clim([0,1])
            %colorbar
            set(gca,'YDir','normal')
            hold off
        end
        
        function plot_velocityfield(self)
            % Plot the velocityfield of the given geometry
            figure(2);
            tiledlayout(1,2)
            nexttile(1)
            imagesc(self.y,self.x,self.Vx)
            hold on
            X_Regler = [685 475]/1000;
            Y_Regler = [395 25]/1000;
            X_Wand = [0 0.800];
            plot(X_Regler,Y_Regler,'k','LineWidth',2)
            plot(X_Wand,[0.6 0.6],'k','LineWidth',1)
            plot(X_Wand,[0.025 0.025],'k','LineWidth',1)
            title('x-Direction')
            clim([0,1])
            set(gca,'YDir','normal')
            colorbar
            hold off
            nexttile(2)
            imagesc(self.y,self.x,self.Vy)
            hold on
            plot(X_Regler,Y_Regler,'k','LineWidth',2)
            plot(X_Wand,[0.6 0.6],'k','LineWidth',1)
            plot(X_Wand,[0.025 0.025],'k','LineWidth',1)
            title('y-Direction')
            clim([0,1])
            set(gca,'YDir','normal')
            colorbar
            hold off
        end

        function fetch_initial_data(self)
            %legacy
            sig0 = 0.02;            % parameter initial density smoothed
            rho = zeros(self.Nx,self.Ny,self.NT);              % Solution Matrix of the Density
            X = self.X;
            Y = self.Y;
            dx = self.dx;
            dy = self.dy;
            rho(:,:,1) = 0;

            rho_max = 2004;
            k = 370;
            agr = size(self.MD_centers);
            N = agr(1);   % N = 192;
            %figure;
            %hold on

            for i=1:N

                part_x = 0.82-self.X_MD(k,6*(i-1)+3)*1.170;
                part_y = 0.615-self.X_MD(k,6*(i-1)+1)*1.045;

                % Trapez Rule for Integration
                if part_x ~= 0.82
                    rho(:,:,1) = rho(:,:,1) + 0.125*1/(rho_max*2*pi*sig0^2)*exp(-1/(2*sig0^2)*((X-part_x-0.0*dx).^2+(Y-part_y-0.5*dy).^2));
                    rho(:,:,1) = rho(:,:,1) + 0.125*1/(rho_max*2*pi*sig0^2)*exp(-1/(2*sig0^2)*((X-part_x-0.0*dx).^2+(Y-part_y+0.5*dy).^2));
                    rho(:,:,1) = rho(:,:,1) + 0.125*1/(rho_max*2*pi*sig0^2)*exp(-1/(2*sig0^2)*((X-part_x+0.5*dx).^2+(Y-part_y-0.0*dy).^2));
                    rho(:,:,1) = rho(:,:,1) + 0.125*1/(rho_max*2*pi*sig0^2)*exp(-1/(2*sig0^2)*((X-part_x-0.5*dx).^2+(Y-part_y+0.0*dy).^2));
                    rho(:,:,1) = rho(:,:,1) + 0.5*1/(rho_max*2*pi*sig0^2)*exp(-1/(2*sig0^2)*((X-part_x).^2+(Y-part_y).^2));
                end
                % plot(part_x,part_y,'o')

            end
            self.rho_init = rho;

        end

        function compute(self)
            self.rho_sol =self.compute_p(self.eps,self.sigma_1,self.gamma,self.rho_init(:,:,1),self.NT);
        end
               
        function rho = compute_p(self,eps,sigma_1,gamma,rho_0,NT)
            %Computes density over time. Depends on data saved in self:
            %initial data, parameters
            rho = zeros(self.Nx,self.Ny,NT);
            rho(:,:,1) = rho_0;
            %             rho_obstacle = rho +self.obstacle;
            rho_s = zeros(self.Nx,self.Ny);

            %sigma_1 = self.sigma_1;
            sup1 = self.sup1; %#ok<*PROPLC>
            eps = self.v_T*eps;
            dx = self.dx;
            dy = self.dy;
            H = @(u) (atan(gamma*(u-1)))/pi + 0.5;     % Approximation of Heaviside function


            % Computation ETA

            [X1, Y1] = meshgrid((-sup1+(dx/2)):dx:sup1,(-sup1+(dy)):dy:sup1-dy);
            [X2, Y2] = meshgrid((-sup1+dx):dx:sup1-dx,(-sup1+(dy/2)):dy:sup1);
            ETA = 1*((sigma_1)/(2*pi))* exp( -0.5*sigma_1*(X2.^2+Y2.^2));

            DETA_X = ((sigma_1)/(2*pi))* exp( -0.5*sigma_1*(X1.^2+Y1.^2)).*(-sigma_1.*X1);
            DETA_Y = ((sigma_1)/(2*pi))* exp( -0.5*sigma_1*(X2.^2+Y2.^2)).*(-sigma_1.*Y2);

            % DETA_X = min(sqrt(X2.^2+Y2.^2)-sup1,0).*(X2./sqrt(X2.^2+Y2.^2));
            % DETA_Y = max(sqrt(X2.^2+Y2.^2)-sup1,0).*(Y2./sqrt(X2.^2+Y2.^2));
            %if done like this, use eps ~= 10000

            self.velocities{1} = zeros(self.Nx,self.Ny,self.NT);
            self.velocities{2} = zeros(self.Nx,self.Ny,self.NT);


            N_eta1 = size(ETA,1)/2;
            N_eta2 = size(ETA,1)/2;

            % Computation Numerical Fluxes
            for k=1:NT-1
                k

                Dx = (dx*dy*conv2(rho(:,:,k)+self.obstacle,DETA_Y));
                Dy = (dx*dy*conv2(rho(:,:,k)+self.obstacle,DETA_X));



                Dx = Dx(N_eta1:end-N_eta1,N_eta2:end-N_eta2+1);
                Dy = Dy(N_eta1:end-N_eta1+1,N_eta2:end-N_eta2);
                %Since DETA_X/DETA_Y have even size, Dx(i) corresponds to
                %Dx(x_(i-1/2)). This property follows through for vx,vy

                vx = -eps.*Dx./sqrt(1+Dx.*Dx + Dy.*Dy);
                vy = -eps.*Dy./sqrt(1+Dx.*Dx + Dy.*Dy);

                self.velocities{1}(:,:,k) = vx;
                self.velocities{2}(:,:,k) = vy;
                sum = 0;

                % Compute F_plus and F_minus
                for j=2:self.Ny

                    for i=(self.Nx-1):-1:2

                        if self.obstacle(i,j) ~= 1
                            %%%F_plus: Flux between (i,j) and (i+1,j), only
                            %%%if both are not an obstacle
                            v = rho(i,j,k);
                            w = rho(i+1,j,k);
                            F_plus = 0;
                            F_minus = 0;

                            if self.obstacle(i+1,j) ~= 1
                                if vx(i+1,j)>0
                                    %if flux is positive, we have flow from
                                    %v to w
                                    F_plus = vx(i+1,j)*v*H(v);
                                else
                                    %if negative, w to v
                                    F_plus = vx(i+1,j)*w*H(w);
                                end
                                if self.Vx(i+1,j)>0
                                    F_plus = F_plus +  self.Vx(i+1,j)*v;
                                else
                                    F_plus = F_plus + self.Vx(i+1,j)*w;
                                end
                            end
                            %%% F_minus: Flux between (i,j) and (i-1,j)
                            v = rho(i-1,j,k);
                            w = rho(i,j,k);
                            if self.obstacle(i-1,j) ~= 1
                                if vx(i,j)>0
                                    F_minus = vx(i,j)*v*H(v);
                                else
                                    F_minus = vx(i,j)*w*H(w);
                                end

                                if self.Vx(i,j)>0
                                    F_minus = F_minus + self.Vx(i,j)*v;
                                else
                                    F_minus = F_minus + self.Vx(i,j)*w;
                                end
                            end
                            %%%
                            rho_s(i,j) = rho(i,j,k) - (self.lambda_x)*(F_plus - F_minus);
                            sum = sum - (self.lambda_x)*(F_plus - F_minus) ;
                        end
                    end
                end
                sum; %#ok<VUNUS>
                for i=2:self.Nx

                    for j=(self.Ny-1):-1:2

                        if self.obstacle(i,j) ~= 1
                            v = rho_s(i,j);
                            w = rho_s(i,j+1);

                            F_plus = 0;
                            F_minus = 0;

                            if self.obstacle(i,j+1) ~= 1
                                if vy(i,j+1) >0
                                    F_plus = vy(i,j+1)*v*H(v);
                                else
                                    F_plus = vy(i,j+1)*w*H(w);
                                end

                                if self.Vy(i,j+1) >0
                                    F_plus = F_plus +  self.Vy(i,j+1)*v;
                                else
                                    F_plus = F_plus + self.Vy(i,j+1)*w;
                                end
                            end

                            %%% F_minus
                            v = rho_s(i,j-1);
                            w = rho_s(i,j);

                            if self.obstacle(i,j-1) ~= 1
                                if vy(i,j) > 0
                                    F_minus = vy(i,j)*v*H(v);
                                else
                                    F_minus = vy(i,j)*w*H(w);
                                end

                                if self.Vy(i,j) > 0
                                    F_minus = F_minus + self.Vy(i,j)*v;
                                else
                                    F_minus = F_minus + self.Vy(i,j)*w;
                                end
                            end

                            %%%

                            rho(i,j,k+1) = rho_s(i,j) - (self.lambda_y)*(F_plus - F_minus);
                            sum = sum - (self.lambda_x)*(F_plus - F_minus) ;
                        end
                    end
                end
                sum; %#ok<VUNUS>
            end
        end

        function compare_parameters(self,bool_par,eps,sigma_1,gamma)
            %Inputs: Lists of length n

            n = length(eps);
            ['Number of simultations:', num2str(n)]
            results = cell(1,n);
            if bool_par
                parfor i = 1:n   %if parfor doesn't work, simply use for instead.
                    ['Simulation number:', num2str(i)]

                    results{i}.rho = self.compute_p(eps(i),sigma_1(i),gamma(i)); %#ok<PFBNS>
                    results{i}.parameters = [eps(i) sigma_1(i) gamma(i)];
                end
            else
                for i = 1:n   %if parfor doesn't work, simply use for instead.
                    ['Simulation number:', num2str(i)] %#ok<*NOPRT>

                    results{i}.rho = self.compute_p(eps(i),sigma_1(i),gamma(i));
                    results{i}.parameters = [eps(i) sigma_1(i) gamma(i)];
                end
            end
            self.results = results;
        end

        function plot_comparison(self,framerate)
            %framerate = 1 highest, choose integer for lower framerate
            fig = figure();%'visible','off');
            fig.Position = [10 10 1000 1000]; %FullHD Display
            n = length(self.results);
            tiledlayout(ceil(sqrt(n+1)),ceil(sqrt(n+1)));
            %nt = round(2/self.dt);
            for k=1:floor((self.NT-1)/(3*framerate))
                %Simulations
                for i = 1:n
                    nexttile(i)
                    imagesc(self.y,self.x,self.results{i}.rho(:,:,k*3*framerate));
                    hold on
                    title(sprintf('$\\varepsilon = %0.1f, \\sigma_1 = %d, \\gamma =%d $',self.results{i}.parameters),'FontSize', 6,'Interpreter','latex');
                    X_Regler = [685 475]/1000;
                    Y_Regler = [395 25]/1000;
                    X_Wand = [0 0.800];
                    %plot(X_Regler,Y_Regler,'k','LineWidth',1)
                    plot(X_Wand,[0.6 0.6],'k','LineWidth',1)
                    plot(X_Wand,[0.025 0.025],'k','LineWidth',1)
                    clim([0,1])
                    colorbar
                    set(gca,'YDir','normal')
                    hold off
                end

                %Data
                nexttile(n+1)
                imagesc(self.y,self.x,self.rho_init(:,:,k*framerate)); %Datapoint
                %only every 3 timesteps
                %imagesc(self.y,self.x,self.obstacle+self.velocities(:,:,k*3));
                hold on
                title('Data')
                X_Regler = [685 475]/1000;
                Y_Regler = [395 25]/1000;
                X_Wand = [0 0.800];
                %plot(X_Regler,Y_Regler,'k','LineWidth',2)
                plot(X_Wand,[0.6 0.6],'k','LineWidth',1)
                plot(X_Wand,[0.025 0.025],'k','LineWidth',1)
                clim([0,1])
                %colorbar
                set(gca,'YDir','normal')
                hold off
                %M(k) = getframe(fig);
                k;
                pause(0.0001)
            end
            %movie(M)
        end

        function plot_rho(self,rho)    
            figure(1);
            %tiledlayout(1,2)
            nt = round(2/self.dt);
            for k=1:(self.NT/3)
                imagesc(self.y,self.x,rho(:,:,3*k));
                hold on
                clim([0,1.5])
                colorbar
                set(gca,'YDir','normal')
                %surf(y,x,rho(:,:,k));
                pause(0.001)
                hold off
    
            end
        end

        function plot_net(self)
            self.plot_rho(self.rho_net)
        end

        function plot(self)
            %%% Plot
            figure(1);
            %tiledlayout(1,2)
            nt = round(2/self.dt);
            for k=1:1200
                %nexttile(1)
                k;
                imagesc(self.y,self.x,self.rho_sol(:,:,3*k));
                hold on
                % X_Regler = [685 475]/1000;
                % Y_Regler = [395 25]/1000;
                % X_Wand = [0 0.800];
                % plot(X_Regler,Y_Regler,'k','LineWidth',2)
                % plot(X_Wand,[0.6 0.6],'k','LineWidth',1)
                % plot(X_Wand,[0.025 0.025],'k','LineWidth',1)
                clim([0,1.5])
                colorbar
                set(gca,'YDir','normal')
                %surf(y,x,rho(:,:,k));
                pause(0.001)
                hold off
                % nexttile(2)
                % imagesc(y,x,rho_data_60(:,:,k));
                % hold on
                % X_Regler = [710 420]/1000;
                % Y_Regler = [305 25]/1000;
                % X_Wand = [0 0.800];
                % plot(X_Regler,Y_Regler,'k','LineWidth',2)
                % plot(X_Wand,[0.6 0.6],'k','LineWidth',1)
                % plot(X_Wand,[0.025 0.025],'k','LineWidth',1)
                % caxis([0,1])
                % colorbar
                % set(gca,'YDir','normal')
                % %surf(y,x,rho(:,:,k));
                % %pause(0.0001)
                % hold off
                % M2(k) = getframe(fig);
            end
        end
        
        function plot_data(self)
            %%% Plot
            figure(1);
            %tiledlayout(1,2)
            %nt = round(2/self.dt);
            for k=1:self.NT/self.sparsityfactor_time_data
                %nexttile(1)

                % nexttile(2)
                imagesc(self.y,self.x,self.rho_init(:,:,k));
                hold on
                % X_Regler = [685 475]/1000;
                % Y_Regler = [395 25]/1000;
                % X_Wand = [0 0.800];
                % plot(X_Regler,Y_Regler,'k','LineWidth',1)
                % plot(X_Wand,[0.6 0.6],'k','LineWidth',1)
                % plot(X_Wand,[0.025 0.025],'k','LineWidth',1)
                
                clim([0,1])
                colorbar
                set(gca,'YDir','normal')
                %surf(y,x,rho(:,:,k));
                pause(0.001)
                hold off
                % M2(k) = getframe(fig);
            end
        end

        function outflow_diagram(self, rhos, times, titles)

            index_out = floor(size(rhos{1},2)) -100 ;                     
            
            for j = 1:length(rhos)
                mass_init = sum(sum(rhos{j}(:,1:index_out,1)));
                T = size(rhos{j},3);
                for i=1:T
                    mass_act = sum(sum(rhos{j}(:,1:index_out,i)));  % Current mass before exit
                    outputs{j}(i) = mass_act/mass_init;
                end                
            end
            figure
            hold on
            for j = 1:length(rhos)
                plot(times{j},outputs{j})
            end
            legend(titles)

        end

        function compare(self)
            %%% Plot
            figure(1);
            tiledlayout(2,1)
            nt = round(2/self.dt);
            for k=1:self.NT-1
                nexttile(1)
                imagesc(self.y,self.x,self.rho_sol(:,:,k*3));
                hold on
                % X_Regler = [685 475]/1000;
                % Y_Regler = [395 25]/1000;
                % X_Wand = [0 0.800];
                % plot(X_Regler,Y_Regler,'k','LineWidth',2)
                % plot(X_Wand,[0.6 0.6],'k','LineWidth',1)
                % plot(X_Wand,[0.025 0.025],'k','LineWidth',1)
                clim([0,1])
                colorbar
                set(gca,'YDir','normal')
                %surf(y,x,rho(:,:,k));
                pause(0.001)
                hold off

                nexttile(2)
                imagesc(self.y,self.x,self.rho_init(:,:,k*3));
                hold on
                % X_Regler = [685 475]/1000;
                % Y_Regler = [395 25]/1000;
                % X_Wand = [0 0.800];
                % plot(X_Regler,Y_Regler,'k','LineWidth',2)
                % plot(X_Wand,[0.6 0.6],'k','LineWidth',1)
                % plot(X_Wand,[0.025 0.025],'k','LineWidth',1)
                clim([0,1])
                colorbar
                set(gca,'YDir','normal')
                %surf(y,x,rho(:,:,k));
                %pause(0.0001)
                hold off
                % M2(k) = getframe(fig);
            end
        end
        
        function fig = wasserstein_distance(self,rho_1,rho_2)
            %% Import necessary Python modules
            scipy_stats = py.importlib.import_module('scipy.stats');
            np = py.importlib.import_module('numpy');
            traceback = py.importlib.import_module('traceback');

            %% use densities from data for all timepoints
            wasserstein_over_time2 = [];
            times2 = [];
            wasserstein_over_time3 = [];
            times3 = [];


        

            rho1_times = rho_1;
            rho2_times = rho_2;
            %rho2_times = MM2.rho_sol(:,:,1:4:end);


            sum_of_mass_1 = sum(rho1_times(:,:,5),[1,2]);
            sum_of_mass_2 = sum(rho2_times(:,:,5),[1,2]);

            for t = 1:4:175
                t
                rho1 = rho1_times(:,:,t);
                rho2 = rho2_times(:,:,t);
                %rho3 = rho_net(:,:,t*5);

                rho1(40,end) = max(0,sum_of_mass_1 - sum(rho1,[1,2]));
                rho2(40,end) = max(0,sum_of_mass_2 - sum(rho2,[1,2]));
                %rho3(40,end) = max(0,sum_of_mass_sol - sum(rho3,[1,2]));

                factor = 0.1;
                rho1 = imresize(rho1,factor,"box");
                rho2 = imresize(rho2,factor,"box");
                %rho3 = imresize(rho3,factor,"box");
                %imagesc(rho2)

                sz = size(rho1);
                X = self.X(1:sz(1),1:sz(2))/factor; Y = self.Y(1:sz(1),1:sz(2))/factor;

                %imagesc(rho2)

                % Flatten the coordinates and densities
                %coords_flat = [MM.X(:),MM.Y(:)];

                coords_flat = [X(:),Y(:)];
                density1_flat = rho1(:);
                density2_flat = rho2(:);
                %density3_flat = rho3(:);

                %Normalize
                density1_flat = density1_flat / sum(density1_flat);
                density2_flat = density2_flat / sum(density2_flat);
                %density3_flat = density3_flat / sum(density3_flat);

                % Convert MATLAB arrays to numpy arrays
                coords_py = np.array(coords_flat);  % Grid coordinates
                density1_py = np.array(density1_flat);
                density2_py = np.array(density2_flat);
                %density3_py = np.array(density3_flat);

                % Compute Wasserstein distance using scipy.stats.wasserstein_distance_nd
                % For some reason, the algorithm fails for some instances. Impose
                % bounds in line 9690 of _stats_py.py to get an approximate solution.
                try
                    wasserstein_distance = scipy_stats.wasserstein_distance_nd(coords_py, coords_py, ...
                        density1_py, density2_py);
                    wasserstein_over_time2 = [wasserstein_over_time2 wasserstein_distance];
                    times2 = [times2 t];
                end
                % try
                % wasserstein_distance = scipy_stats.wasserstein_distance_nd(coords_py, coords_py, ...
                %     density1_py, density3_py);
                % wasserstein_over_time3 = [wasserstein_over_time3 wasserstein_distance];
                % times3 = [times3 t];
                % end
                %wasserstein_over_time = [wasserstein_over_time wasserstein_distance];
            end
            fig = figure(5);
            plot(5*self.dt*times2,wasserstein_over_time2)%,times3,wasserstein_over_time3)
            xlabel("Time in s");
            ylabel("Wasserstein Distance");
        end


        %% To be put in subclass
        function samples = compute_flux(self,samples)

            gamma = self.gamma;
            sigma = self.sigma_1;
            eps = self.eps;

            scopesize_y = 10;
            scopesize_x = 10;
            
            %eps = self.v_T*eps;
            dx = self.dx;
            dy = self.dy;

            sup1 = scopesize_x*dx; %#ok<*PROPLC>
            H = @(u) (atan(gamma*(u-1)))/pi + 0.5;     % Approximation of Heaviside function

            % Computation ETA
            %sup1 has to depend on the sample
            [X, Y] = meshgrid((-sup1):dx:sup1,(-sup1):dy:sup1);
            %[X2, Y2] = meshgrid((-sup1):dx:sup1,(-sup1):dy:sup1);
            %ETA = 1*((sigma)/(2*pi))* exp( -0.5*sigma*(X2.^2+Y2.^2));
            % ATTENTION: We adopted the mesh to represent
            % the flow on the lattice, not the interfaces.

            DETA_X = ((sigma)/(2*pi))* exp( -0.5*sigma*(X.^2+Y.^2)).*(-sigma.*X);
            DETA_Y = ((sigma)/(2*pi))* exp( -0.5*sigma*(X.^2+Y.^2)).*(-sigma.*Y);

            parfor i = 1:length(samples)

                %sample = samples

                Dx = (dx*dy*conv2(samples{i}.density+samples{i}.obstacle,DETA_Y,'valid'));
                Dy = (dx*dy*conv2(samples{i}.density+samples{i}.obstacle,DETA_X,'valid'));

                vx = Dx;
                vy = Dy;

                % vx =  - eps.*Dx./sqrt(1+Dx.*Dx + Dy.*Dy);
                % vy = -eps.*Dy./sqrt(1+Dx.*Dx + Dy.*Dy);

                samples{i}.velocity_computed =[vy;vx];
                
            end

        end
        
        function [vx,vy] = apply_net(self,rho,net)
            %net_x gives velocities orthogonal to the flow
            %net_y gives velocities in direction of the flow
            rho_padded = padarray(rho,[10,10],0);
            obst_padded = padarray(self.obstacle,[10,10],1);
            NNinput = cat(3,rho_padded,obst_padded);
            vx = zeros(size(rho));
            vy = zeros(size(rho));
            v_T = self.v_T; Nx = self.Nx-1; Ny =  self.Ny-1;
            Factor = 1.2;
            k=1;
            loc = false(self.Nx,self.Ny);
            Inputs = zeros(21,21,2,1000);
            for j=2:Ny
                for i=2:Nx
                    if max(rho(i-1:i+1,j-1:j+1),[],'all')>0.8 
                        % a = predict(net,NNinput(i:i+20,j:j+20,:));
                        % vx(i,j) = Factor*v_T*a(1);
                        % vy(i,j) = Factor*v_T*a(2);    
                        Inputs(:,:,:,k) = NNinput(i:i+20,j:j+20,:);  
                        k = k+1;
                        loc(i,j)=true;
                    end
                end
            end
            if k>1
                v = predict(net,Inputs(:,:,:,1:k-1));
                vx(loc==1) = Factor*v_T*v(:,1);
                vy(loc==1) = Factor*v_T*v(:,2);
            end

        end

        function rho = compute_continue_NN(self,net,rho_old,NT)
            %Continues computation of rho_old by NT time steps

            rho_continued = self.compute_Macro_NN_(net,rho_old(:,:,end),NT);
            rho = cat(3,rho_old,rho_continued(:,:,2:end));
        end
        
        function rho = compute_Macro_NN(self,net,NT)
            rho = self.compute_Macro_NN_(net,self.rho_init(:,:,1),NT);
        end

        function rho = compute_Macro_NN_(self,net,rho_0,NT)

            %Computes density over time. Depends on data saved in self:
            %initial data
            %net_x gives velocities in 
            rho = zeros(self.Nx,self.Ny,NT);
            rho(:,:,1) =rho_0 ;
            %             rho_obstacle = rho +self.obstacle;
            rho_s = zeros(self.Nx,self.Ny);
            %sigma_1 = self.sigma_1;
            %sup1 = self.sup1; %#ok<*PROPLC>
            %eps = self.v_T*eps;
            dx = self.dx;
            dy = self.dy;
            H = @(u) (atan(self.gamma*(u-1)))/pi + 0.5;     % Approximation of Heaviside function

            % DETA_X = min(sqrt(X2.^2+Y2.^2)-sup1,0).*(X2./sqrt(X2.^2+Y2.^2));
            % DETA_Y = max(sqrt(X2.^2+Y2.^2)-sup1,0).*(Y2./sqrt(X2.^2+Y2.^2));
            %if done like this, use eps ~= 10000

            self.velocities_NN{1} = zeros(self.Nx,self.Ny,self.NT);
            self.velocities_NN{2} = zeros(self.Nx,self.Ny,self.NT);

            % Computation Numerical Fluxes
            for k=1:NT
                k   
                % if k == 50
                %     a = 3;
                % end
                [vx,vy] = self.apply_net(min(1,rho(:,:,k)),net); %Velocities in cell center
                vx(2:end,:) = .5*(vx(1:end-1,:)+vx(2:end,:)); %vx(i,j) is approxiate velocity on interface i-1/2,j
                vy(:,2:end) = .5*(vy(:,1:end-1)+vy(:,2:end)); %vx(i,j) is approxiate velocity on interface i,j-1/2

                % vx = vx + self.Vx;
                % vy = vy + self.Vy;

                self.velocities_NN{1}(:,:,k) = vx; %orth to flow
                self.velocities_NN{2}(:,:,k) = vy; %direction of flow
                sum = 0;

                % Compute F_plus and F_minus
                for j=2:self.Ny

                    for i=(self.Nx-1):-1:2

                        if self.obstacle(i,j) ~= 1
                            %%%F_plus: Flux between (i,j) and (i+1,j), only
                            %%%if both are not an obstacle
                            v = rho(i,j,k);
                            w = rho(i+1,j,k);
                            F_plus = 0;
                            F_minus = 0;

                            if self.obstacle(i+1,j) ~= 1
                                if vx(i+1,j)>0
                                    %if flux is positive, we have flow from
                                    %v to w
                                    F_plus = vx(i+1,j)*v*H(v);
                                else
                                    %if negative, w to v
                                    F_plus = vx(i+1,j)*w*H(w);
                                end
                                if self.Vx(i+1,j)>0
                                    F_plus = F_plus +  self.Vx(i+1,j)*v;
                                else
                                    F_plus = F_plus + self.Vx(i+1,j)*w;
                                end
                            end
                            %%% F_minus: Flux between (i,j) and (i-1,j)
                            v = rho(i-1,j,k);
                            w = rho(i,j,k);
                            if self.obstacle(i-1,j) ~= 1
                                if vx(i,j)>0
                                    F_minus = vx(i,j)*v*H(v);
                                else
                                    F_minus = vx(i,j)*w*H(w);
                                end

                                if self.Vx(i,j)>0
                                    F_minus = F_minus + self.Vx(i,j)*v;
                                else
                                    F_minus = F_minus + self.Vx(i,j)*w;
                                end
                            end
                            %%%
                            rho_s(i,j) = rho(i,j,k) - (self.lambda_x)*(F_plus - F_minus);
                            sum = sum - (self.lambda_x)*(F_plus - F_minus) ;
                        end
                    end
                end
                sum; %#ok<VUNUS>
                for i=2:self.Nx

                    for j=(self.Ny-1):-1:2

                        if self.obstacle(i,j) ~= 1
                            v = rho_s(i,j);
                            w = rho_s(i,j+1);

                            F_plus = 0;
                            F_minus = 0;

                            if self.obstacle(i,j+1) ~= 1
                                if vy(i,j+1) >0
                                    F_plus = vy(i,j+1)*v*H(v);
                                else
                                    F_plus = vy(i,j+1)*w*H(w);
                                end

                                if self.Vy(i,j+1) >0
                                    F_plus = F_plus +  self.Vy(i,j+1)*v;
                                else
                                    F_plus = F_plus + self.Vy(i,j+1)*w;
                                end
                            end

                            %%% F_minus
                            v = rho_s(i,j-1);
                            w = rho_s(i,j);

                            if self.obstacle(i,j-1) ~= 1
                                if vy(i,j) > 0
                                    F_minus = vy(i,j)*v*H(v);
                                else
                                    F_minus = vy(i,j)*w*H(w);
                                end

                                if self.Vy(i,j) > 0
                                    F_minus = F_minus + self.Vy(i,j)*v;
                                else
                                    F_minus = F_minus + self.Vy(i,j)*w;
                                end
                            end

                            %%%

                            rho(i,j,k+1) = rho_s(i,j) - (self.lambda_y)*(F_plus - F_minus);
                            sum = sum - (self.lambda_x)*(F_plus - F_minus) ;
                        end
                    end
                end
                sum; %#ok<VUNUS>
            end
            self.rho_net = rho;
           
         
        end

    end
    

end

