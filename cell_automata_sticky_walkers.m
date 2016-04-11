%% Computational Modelling assignment #6 - cellular automata

function [ t_final, r_final, Domain ] = ...
    random_walkers_aggradation( L, N, T, plot )
%%-------------------------------------------------------------------------
% K. Kochanski, 11-Apr-2016
%--------------------------------------------------------------------------
%% DESCRIPTION 
% This code models a 'sticky' random walker.
%   The code shows a particle which performs a 2D random walk around the
% screen. The particle stops either: after time T has passed (particle
% 'dies'), or when it tries to move on top of a previous particle (particle
% is 'stuck'). Particles cannot move off the screen.
%
% Walkers can move in any of 8 directions, or be still, at each time step.
% The effective 'walker diffusivity', or mean cells moved per step, is 1.07
%
% Everything is dimensionless.
%
% PLOTTING
% To show the walkers walking (every time step)   : plot = 1
% To show only the final result (runs much faster): plot = 0
%
%-------------------------------------------------------------------------
%% Input variables
% -var------|--meaning----------|--type---------|--suggested value--------
% L         | lattice size      | integer       |  15
% N         | numner of walkers | integer       |  10
% T         | lifetime of walker| integer       |  210
% plot      | plot every step?  | boolean       |  1
%%------------------------------------------------------------------------

% Effective diffusion length scale: average distance travelled by a walker,
% normalized by the lattice size L
diff_length_scale = sqrt(1.07*T)/L;
disp(sprintf('Length scale of diffusion : %d', diff_length_scale));

%Create the domain. Every zero space is empty, every 1 space is full
Domain = zeros(L);
% Start with a single full space (the pole) in the center
Domain( round((L+1)./2), round((L+1)./2) ) = 1;

%Make lists to store useful details:
% the time each walker lasts for, t_final
t_final = zeros(1,N);
% and the radius each walker reaches, r_final
r_final = zeros(1,N);

%Create and walk N walkers
for ii = 1:N
    %Generate a random initial start point on the edge
    edge0 = randi(4); %Picks one of the 4 edges
    l0 = randi(L); %Picks a distance along the edge
    if edge0 == 1
        position = [l0 1];
    elseif edge0 == 2
        position = [l0 L];
    elseif edge0 == 3
        position = [1 l0];
    elseif edge0 == 4
        position = [L l0];
    end        
    
    % 'Domain' stores all the stuck particles in their final states
    % 'tempdomain' *also* stores the moving particle at each time step,
    % so i can use it in plotting.
    tempDomain = Domain; 
    
    t = 0; % Counts the number of steps
    stuck = 0; % Determines whether this walker has stuck yet
    %Walk until the walker sticks to a non-empty spot
    
    while stuck == 0;
        if plot == 1;
            %% plot walker position at each step
            % walker color changes with age
            walkercolor = 0.1 + 0.9*t/T;
            % Mark the new position on the plot
            tempDomain(position(1), position(2)) = walkercolor; 
            image(tempDomain, 'CDataMapping', 'scaled')
            axis off
            pause(0.01)
            % (and erase it from the plot, in preparation for next step)
            tempDomain(position(1), position(2)) = 0;
        end
        
        % Generate a 'movement' in x and y
        %by generating two random numbers which are either -1, 0 or +1
        % (Note: I am letting the walker move in any of 8 directions)
        new_position = position + ...
            [ randi(3)-2, randi(3)-2];
        
        % Check whether the new position is in the allowed domain
        % If it isn't, advance time without moving.
        if new_position(1) <= 0  ||  new_position(1) > L
            t = t + 1;
        elseif new_position(2) <= 0  ||  new_position(2) > L
            t = t + 1;
            
        % Then check whether the new positions is 'sticky'
        % If it is, make the current position sticky, flash its color as a 
        % visual indicator, and make kill the walker.
        elseif Domain(new_position(1), new_position(2)) == 1;
            t = t + 1;
            if plot == 1;
                tempDomain(new_position(1), new_position(2)) = walkercolor;
                tempDomain(position(1), position(2)) = walkercolor;
                image(tempDomain, 'CDataMapping','scaled')
                axis off
                pause(0.25)
            end
            Domain(position(1), position(2)) = 1;
            stuck = 1;
            
        % And if there is nothing exciting about the new spot, move walker
        else
            t = t + 1;
            position = new_position;
        end
        
        % Walker runs out of time after time T
        if t >= T
            Domain(position(1), position(2)) = 1;
            stuck = 1;
        end
       
    end
    
    % Store the final time and final radius reached
    t_final(ii) = t;
    r = (position(1)-50).^2 + (position(2)-50).^2;
    r_final(ii) = r;
       
    
    
end

%Make a plot of the final aggregation
image(Domain, 'CDataMapping','scaled')
axis off

end

