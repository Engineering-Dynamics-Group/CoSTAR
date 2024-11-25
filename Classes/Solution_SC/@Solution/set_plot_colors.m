 %This is a static method of the superclass Solution. It sets different standard colors to the plot_color property
 %of the Solution Superclass for plotting in contplot and solplot
 %
 %@obj: Solution subclass object
    
function set_plot_colors(obj)
    

    obj.plot_color.blue             = [ 46/255, 116/255, 181/255]; 		%blue
    obj.plot_color.green            = [112/255, 173/255,  71/255]; 		%green
    obj.plot_color.red              = [197/255,   0/255,  90/255]; 	    %red

    obj.plot_color.yellow           = [255/255, 192/255,   0/255];      %yellow
    obj.plot_color.cyan             = [ 40/255, 198/255, 198/255];      %cyan

    obj.plot_color.magenta          = [230/255,   0/255, 229/255];      %magenta
    obj.plot_color.orange           = [237/255, 125/255,  49/255];      %orange

    obj.plot_color.purple           = [112/255,  48/255, 160/255];  	%purple
    obj.plot_color.grey             = [128/255, 137/255, 146/255];      %grey

    obj.plot_color.black            = [   0   ,    0   ,    0   ];      %black

    %%%%%%%%%%%%%%%%%

    obj.plot_color.pale_blue        = [165/255, 200/255, 232/255];      %pale blue
    obj.plot_color.pale_green       = [198/255, 225/255, 181/255]; 	    %pale green
    obj.plot_color.pale_red 		= [255/255, 192/255, 221/255]; 	    %pale red

    obj.plot_color.pale_yellow      = [255/255, 221/255, 119/255];      %pale yellow
    obj.plot_color.pale_cyan        = [149/255, 234/255, 234/255];      %pale cyan

    obj.plot_color.pale_magenta     = [255/255, 143/255, 225/255];      %pale magenta
    obj.plot_color.pale_orange      = [244/255, 177/255, 131/255];      %pale orange

    obj.plot_color.pale_purple      = [199/255, 162/255, 227/255];  	%pale purple
    obj.plot_color.pale_grey        = [186/255, 191/255, 195/255];      %pale grey

    %%%%%%%%%%%%%%%%

    obj.plot_color.dark_blue        = [  0/255,  48/255,  90/255];      %dark blue
    obj.plot_color.dark_green       = [ 78/255, 121/255,  50/255]; 	    %dark green
    obj.plot_color.dark_red         = [139/255,  17/255,  65/255];      %dark red
    
    obj.plot_color.dark_yellow      = [179/255, 134/255,   0/255];      %dark yellow
    obj.plot_color.dark_cyan        = [24/255,  119/255, 119/255];      %dark cyan

    obj.plot_color.dark_magenta     = [160/255,   0/255, 161/255];      %dark_magenta
    obj.plot_color.dark_orange      = [184/255,  84/255,  16/255];      %dark orange
   
    obj.plot_color.dark_purple      = [190/255, 147/255, 222/255];  	%dark purple
    obj.plot_color.dark_grey        = [ 74/255,  83/255,  92/255];      %dark grey
    
end