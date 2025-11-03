import math

def magnitude (value):
    if (value == 0): return 0
    return int(math.floor(math.log10(abs(value))))