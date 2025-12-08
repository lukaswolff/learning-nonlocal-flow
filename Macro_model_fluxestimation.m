classdef Macro_model_fluxestimation < Macro_model


    properties

    end


    methods
        function self = Macro_model_fluxestimation(alpha,options)
            % Call Macro_model constructor
            arguments
                alpha double
                options.obstacle = {}
                options.parameters = [1.5,2000,20]
            end

            self@Macro_model(alpha,obstacle = options.obstacle,parameters = options.parameters);
            self.load_setup(0);
        end

        function samples = compute_flux(samples)
            
            scopesize_y = 10;
            scopesize_x = 10;
            sup1 = self.sup1; %#ok<*PROPLC>
            %eps = self.v_T*eps;
            dx = self.dx;
            dy = self.dy;
            H = @(u) (atan(gamma*(u-1)))/pi + 0.5;     % Approximation of Heaviside function

            % Computation ETA
            %sup1 has to depend on the sample
            [X1, Y1] = meshgrid((-sup1+(dx/2)):dx:sup1,(-sup1+(dy)):dy:sup1-dy);
            [X2, Y2] = meshgrid((-sup1+dx):dx:sup1-dx,(-sup1+(dy/2)):dy:sup1);
            ETA = 1*((sigma_1)/(2*pi))* exp( -0.5*sigma_1*(X2.^2+Y2.^2));

            DETA_X = ((sigma_1)/(2*pi))* exp( -0.5*sigma_1*(X1.^2+Y1.^2)).*(-sigma_1.*X1);
            DETA_Y = ((sigma_1)/(2*pi))* exp( -0.5*sigma_1*(X2.^2+Y2.^2)).*(-sigma_1.*Y2);

            N_eta1 = size(ETA,1)/2;
            N_eta2 = size(ETA,1)/2;

            Dx = (dx*dy*conv2(rho(:,:,k)+self.obstacle,DETA_Y));
            Dy = (dx*dy*conv2(rho(:,:,k)+self.obstacle,DETA_X));



            Dx = Dx(N_eta1:end-N_eta1,N_eta2:end-N_eta2+1);
            Dy = Dy(N_eta1:end-N_eta1+1,N_eta2:end-N_eta2);
            %Since DETA_X/DETA_Y have even size, Dx(i) corresponds to
            %Dx(x_(i-1/2)). This property follows through for vx,vy

            vx = -eps.*Dx./sqrt(1+Dx.*Dx + Dy.*Dy);
            vy = -eps.*Dy./sqrt(1+Dx.*Dx + Dy.*Dy);

        end

        function apply_net(self,rho,net_x,net_y)
            rho_padded = padarray(rho,[10,10],1);
            vx = zeros(size(rho));
            vy = zeros(size(rho));
            for j=1:self.Ny
                for i=1:self.Nx
                    if rho(i,j)>0.9
                        vx(i,j) = predict(net_x,rho_padded(i:i+20,j,j+20));
                        vy(i,j) = predict(net_y,rho_padded(i:i+20,j,j+20));

                    end
                end
            end

        end

        function compute_Macro_NN(self,net_x,net_y)

                        %Computes density over time. Depends on data saved in self:
            %initial data, parameters
            rho = zeros(self.Nx,self.Ny,self.NT);
            rho(:,:,1) = self.rho_init(:,:,1);
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
            for k=1:self.NT-1
                %k;

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


    end

end
