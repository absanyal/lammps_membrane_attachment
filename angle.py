import numpy as np
import matplotlib.pyplot as plt

def filter_angle(theta):
    if 0 <= theta <= 90:
        return theta
    elif theta < 0:
        return filter_angle(theta + 90)
    elif theta > 90:
        return filter_angle(theta - 90)