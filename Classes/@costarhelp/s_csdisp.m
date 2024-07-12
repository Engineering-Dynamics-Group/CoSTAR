% This is a static method of class costarhelp
% 
% This wrapper function for the display function takes the string 
% given in "mystring", splits it at the Latex style line break 
% delimiter "\\" and displays the splitted string.
%
%@mystring:     String character containing "\\" for line break

function s_csdisp(mystring)
    
    splitted_string = split(mystring,["\\"]);   %Split the string at the Latex Style delimite "\\"
    
    
    for k = 1:length(splitted_string)           %splitted_string is now a cell. Display the single cell lines 
        disp(splitted_string{k,1})
    end


end