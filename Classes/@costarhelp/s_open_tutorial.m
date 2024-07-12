% Little help function to open linked tutorial or minimal example

function s_open_tutorial(file_name)

    help_struct.filename = char(file_name);

    edit(file_name)

end

