 %This function displays a clickable hyperlink called name, which executes the function fcn
 %
 %@fcn:     executable function statement as string
 %@name:    string for name to be displayed
function s_disp_fcncall(fcn,name)


    fprintf(append('<a href="matlab:',fcn,'">',name,'</a>'));

end