%This function is a method of superclass SOLUTION.
%It gives back a rgb color code defined for plotting in contplot or solplot
%either defined by the user or automatically
%
%@obj:      Solution Superclass object
%@options:  options structure for contplot or solplot

function rgb = custom_color(obj,options)

    if isfield(options,'color') %gatekeeper ensured, that if the field is present, it has the right format

        if ischar(options.color) %choose color according to the Matlab standard colors
            if strcmpi(options.color,'r');  rgb =  obj.plot_color.red;        end
            if strcmpi(options.color,'g');  rgb =  obj.plot_color.green;      end
            if strcmpi(options.color,'b');  rgb =  obj.plot_color.blue;       end
            if strcmpi(options.color,'c');  rgb =  obj.plot_color.cyan;       end
            if strcmpi(options.color,'m');  rgb =  obj.plot_color.magenta;    end
            if strcmpi(options.color,'y');  rgb =  obj.plot_color.yellow;     end
            if strcmpi(options.color,'k');  rgb =  obj.plot_color.black;      end

        else  %An rgb code for the color was supplied
            rgb = options.color;
        end


    else %no color was supplied. Select a new color for every new plot in current figure from the predefined colors in the class

        fn = fieldnames(obj.plot_color);
        numplots =  mod(numel(get(gca,'Children')),numel(fn))+1;    %Gets the number of single plots present in the currently open figure (if there is one): 
                                                                   %This way, a new color gets selected for every new plot. mod(...) starts the cycle if there are more plots than colors.
        rgb = obj.plot_color.(fn{numplots,1});

    end

end