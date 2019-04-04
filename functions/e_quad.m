function perform = e_quad(dts,dte)

    Y = cell(length(dte.y),1);
    for w = 1 : length(dte.y)
    
        ind = zeros(length(dte.y{w}),1);
        for k = 1 : length(dte.y{w})
            error = abs(dte.t{w}(k) - dts.td{w});
            [~,ind(k)] = min(error);
        end
        
        Y{w} = dts.yd{w}(ind);
    end
    
    Y_1 = cell2mat(Y); Y_2 = cell2mat(dte.y);
    
    perform.quadratic = (Y_1 - Y_2)'*(Y_1 - Y_2)/length(Y_1);
end