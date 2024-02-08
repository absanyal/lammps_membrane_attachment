import numpy as np
import matplotlib.pyplot as plt

def filter_angle(theta):
    if 0 <= theta <= 90:
        return theta
    elif 90 < theta <= 180:
        return filter_angle(theta - 90)
    elif theta > 180:
        return filter_angle(theta - 180)
    elif -90 <= theta < 0:
        return filter_angle(theta + 90)
    elif -180 <= theta < -90:
        return filter_angle(theta + 180)
    elif theta < -180:
        return filter_angle(theta + 180)