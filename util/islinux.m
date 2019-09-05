function result = islinux()

    if (~ismac && isunix)
        result = true;
    else
        result = false;
    end

end