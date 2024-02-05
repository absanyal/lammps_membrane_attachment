import numpy as np


def moving_average(array, window_size, padding='constant'):
    """
    Calculates the moving average of an array with padding.

    Args:
      array: The array to calculate the moving average of.
      window_size: The size of the moving window.
      padding: The type of padding to use.

    Returns:
      The moving average of the array.
    """

    if padding not in ['constant', 'reflect', 'edge']:
        raise ValueError(
            'Padding must be one of "constant", "reflect", or "edge".')

    if window_size > len(array):
        raise ValueError(
            'Window size must be less than or equal to the length of the array.')

    padding_width = window_size // 2
    if padding == 'constant':
        padded_array = np.pad(array, padding_width, 'constant')
    elif padding == 'reflect':
        padded_array = np.pad(array, padding_width, 'reflect')
    elif padding == 'edge':
        padded_array = np.pad(array, padding_width, 'edge')

    moving_average = np.convolve(padded_array, np.ones(
        window_size), 'valid') / window_size

    return moving_average


class cylinder():
    """
    Represents a cylinder in 3D space.

    Attributes:
        radius (float): The radius of the cylinder.
        rA (vector): The position vector of end-point A on the long axis if the cylinder.
        rB (vector): The position vector of end-point B on the long axis if the cylinder.
    """

    def __init__(self, radius, rA, rB):
        """
        Initialize the CylinderMath class.

        Args:
            radius (float): The radius of the cylinder.
            rA (float): Position of end-point A.
            rB (float): Position of end-point B.
        """
        self.radius = radius
        self.rA = rA
        self.rB = rB


def norm(A):
    """
    Calculate the Euclidean norm of a given vector.

    Parameters:
    A (list): The input vector.

    Returns:
    float: The Euclidean norm of the input vector.
    """
    sum = 0
    for i in range(len(A)):
        sum += A[i] * A[i]
    return np.sqrt(sum)


def distance_from_axis(cyl, rP):
    """
    Calculate the distance from the axis of a cylinder to a given point.

    Args:
        cyl: The cylinder object with attributes rA and rB representing two points
             defining the axis of the cylinder.
        rP: The point for which the distance from the cylinder axis is calculated.

    Returns:
        The distance from the cylinder axis to the given point.
    """
    e = cyl.rA - cyl.rB
    d = norm((np.cross(e, rP - cyl.rA)))/norm(e)
    return d


def distance_from_surface(cyl, rP):
    """
    Calculate the distance from the surface of a cylinder to a given point.

    Args:
        cyl: the cylinder object
        rP: the point from which to calculate the distance

    Returns:
        The distance from the surface of the cylinder to the given point
    """
    return cyl.radius - distance_from_axis(cyl, rP)
