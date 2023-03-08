# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 09:25:06 2020

@author: modified by esoria
"""
# -*- coding: utf-8 -*-
from __future__ import absolute_import

from __future__ import division

import numpy as np

from skimage.transform import warp


from warnings import warn

import math







def radon(image, theta=None, circle=None):

    """

    Calculates the radon transform of an image given specified

    projection angles.



    Parameters

    ----------

    image : array_like, dtype=float

        Input image. The rotation axis will be located in the pixel with

        indices ``(image.shape[0] // 2, image.shape[1] // 2)``.

    theta : array_like, dtype=float, optional (default np.arange(180))

        Projection angles (in degrees).

    circle : boolean, optional

        Assume image is zero outside the inscribed circle, making the

        width of each projection (the first dimension of the sinogram)

        equal to ``min(image.shape)``.

        The default behavior (None) is equivalent to False.



    Returns

    -------

    radon_image : ndarray

        Radon transform (sinogram).  The tomography rotation axis will lie

        at the pixel index ``radon_image.shape[0] // 2`` along the 0th

        dimension of ``radon_image``.

    xp: coordinates

    References

    ----------

    .. [1] AC Kak, M Slaney, "Principles of Computerized Tomographic

           Imaging", IEEE Press 1988.

    .. [2] B.R. Ramesh, N. Srinivasa, K. Rajgopal, "An Algorithm for Computing

           the Discrete Radon Transform With Some Applications", Proceedings of

           the Fourth IEEE Region 10 International Conference, TENCON '89, 1989



    Notes

    -----

    Based on code of Justin K. Romberg

    (http://www.clear.rice.edu/elec431/projects96/DSP/bpanalysis.html)



    """

    if image.ndim != 2:

        raise ValueError('The input image must be 2-D')

    if theta is None:

        theta = np.arange(180)

    if circle is None:

        warn('The default of `circle` in `skimage.transform.radon` '

             'will change to `True` in version 0.15.')

        circle = False



    if circle:

        radius = min(image.shape) // 2

        c0, c1 = np.ogrid[0:image.shape[0], 0:image.shape[1]]

        reconstruction_circle = ((c0 - image.shape[0] // 2) ** 2

                                 + (c1 - image.shape[1] // 2) ** 2)

        reconstruction_circle = reconstruction_circle <= radius ** 2

        if not np.all(reconstruction_circle | (image == 0)):

            warn('Radon transform: image must be zero outside the '

                 'reconstruction circle')

        # Crop image to make it square

        slices = []

        for d in (0, 1):

            if image.shape[d] > min(image.shape):

                excess = image.shape[d] - min(image.shape)

                slices.append(slice(int(np.ceil(excess / 2)),

                                    int(np.ceil(excess / 2)

                                        + min(image.shape))))

            else:

                slices.append(slice(None))

        slices = tuple(slices)

        padded_image = image[slices]

    else:

        diagonal = np.sqrt(2) * max(image.shape)

        pad = [int(np.ceil(diagonal - s)) for s in image.shape]

        new_center = [(s + p) // 2 for s, p in zip(image.shape, pad)]

        old_center = [s // 2 for s in image.shape]

        pad_before = [nc - oc for oc, nc in zip(old_center, new_center)]

        pad_width = [(pb, p - pb) for pb, p in zip(pad_before, pad)]

        padded_image = np.pad(image, pad_width, mode='constant',

                              constant_values=0)

    # padded_image is always square

    assert padded_image.shape[0] == padded_image.shape[1]

    radon_image = np.zeros((padded_image.shape[0], len(theta)))

    center = padded_image.shape[0] // 2
    
    b = np.ceil (math.sqrt (sum (np.array([np.size(padded_image,0),np.size(padded_image,1)])**2))/2 + 1)
    xp = np.arange(-b,b+1)


    shift0 = np.array([[1, 0, -center],

                       [0, 1, -center],

                       [0, 0, 1]])

    shift1 = np.array([[1, 0, center],

                       [0, 1, center],

                       [0, 0, 1]])



    def build_rotation(theta):

        T = np.deg2rad(theta)

        R = np.array([[np.cos(T), np.sin(T), 0],

                      [-np.sin(T), np.cos(T), 0],

                      [0, 0, 1]])

        return shift1.dot(R).dot(shift0)



    for i in range(len(theta)):

        rotated = warp(padded_image, build_rotation(theta[i]))

        radon_image[:, i] = rotated.sum(0)
        b = np.ceil (np.size(radon_image,0)/2 - 1)
        xp = np.arange(-b,b+1)

    return radon_image, xp


   
    
