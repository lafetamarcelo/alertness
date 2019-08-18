function plotAnnealingResults(dt, fval)
global iter__  bVal__ dcurves__ ncurves__ ind__

iter__ = iter__ + 1; 

switch iter__
    case 1
        bVal__ = fval;
        plotnew = true;
    otherwise
        if fval < bVal__
         bVal__ = fval;
         plotnew = true;
        else
         plotnew = false; 
        end
end

% If is a new best value
if plotnew
  
  % get the size from 
  ns = length(dt.dayTime);
  
  if ind__ > 1
   for ki = 1 : ind__ - 1
    for kii = 1 : ns
      %set(gca,'Children',[dcurves__{ki,kii} ncurves__{ki,kii}])
      dcurves__{ki,kii}.Color(4) = 0.05;
      if kii < ns
        ncurves__{ki,kii}.Color(4) = 0.05;
      end
    end
   end 
  end
  
  h = figure(1); hold on;
  axis tight manual
  filename = 'realAlert.gif';
  for k = 1 : ns
    % Plot the day results
    dcurves__{ind__,k} = plot(dt.dayTime{k}, dt.dayOut{k}, 'r', 'LineWidth', 2.0);
    dcurves__{ind__,k}.Color(4) = 0.05;

    if k ~= ns
      % Plot the nigth results
      ncurves__{ind__,k} = plot(dt.nightTime{k}, dt.nightOut{k}, 'b', 'LineWidth', 2.0);
      ncurves__{ind__,k}.Color(4) = 0.05;  
    end
    
  end
  hold off;
  
  % Create the GIF
  
  % Capture the plot as an image 
  frame = getframe(h); 
  im = frame2im(frame); 
  [imind,cm] = rgb2ind(im,256); 
  % Write to the GIF File 
  if ind__ == 1 
    imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
  else 
    imwrite(imind,cm,filename,'gif','WriteMode','append'); 
  end 
  
  ind__ = ind__ + 1;
end

end