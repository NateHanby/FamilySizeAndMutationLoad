# checks if a string of text can exist in a partially completed box that accetps only floats.  
def isfloat(input):
    if input == "" or input == ".":
        return True
    
    try:
        float(input)
    except ValueError:
        return False
    
    return True

# checks if a string of text can exist in a partially completed box that accepts only integers.  
def isint(input):
    if input == "":
        return True
    
    try:
        int(input)
    except ValueError:
        return False
    
    return True