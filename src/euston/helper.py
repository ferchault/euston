import math

def human_readable_time(seconds):
    minutes = int(math.floor(seconds / 60.0))
    seconds -= minutes*60
    hours = int(math.floor(minutes / 60))
    minutes -= hours*60
    days = int(math.floor(hours / 24))

    return '%03d-%02d:%02d:%06.3f' % (days, hours, minutes, float(seconds))