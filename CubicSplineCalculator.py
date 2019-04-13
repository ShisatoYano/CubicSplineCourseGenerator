# -*- coding: utf-8 -*-

import numpy as np
import sys
import math
import matplotlib.pyplot as plt

class  CubicSplineCalculator(object):
    def __init__(self, ndarray_x, ndarray_y):
        self._sample_x = ndarray_x # ndarray
        self._sample_y = ndarray_y # ndarray
        self._calc_coefficient()
    
    def _calc_coefficient(self):
        point_num = len(self._sample_x)
        N         = point_num - 1
        matrix    = np.zeros([4*N, 4*N])
        y         = np.zeros([4*N])
        
        # simultaneous equations
        mat_idx = 0
        for i in range(N):
            for j in range(4):
                matrix[mat_idx, 4*i+j] = pow(self._sample_x[i], j)
            y[mat_idx] = self._sample_y[i]
            mat_idx += 1
        for i in range(N):
            for j in range(4):
                matrix[mat_idx, 4*i+j] = pow(self._sample_x[i+1], j)
            y[mat_idx] = self._sample_y[i+1]
            mat_idx += 1
        for i in range(N-1):
            for j in range(4):
                matrix[mat_idx, 4*i+j] = j*pow(self._sample_x[i+1], j-1)
                matrix[mat_idx, 4*(i+1)+j] = -j*pow(self._sample_x[i+1], j-1)
            mat_idx += 1
        for i in range(N-1):
            matrix[mat_idx, 4*i+3] = 3*self._sample_x[i+1]
            matrix[mat_idx, 4*i+2] = 1
            matrix[mat_idx, 4*(i+1)+3] = -3*self._sample_x[i+1]
            matrix[mat_idx, 4*(i+1)+2] = -1
            mat_idx += 1
        matrix[mat_idx,3] = 3*self._sample_x[0]
        matrix[mat_idx,2] = 1
        mat_idx += 1
        matrix[4*N-1,4*N-1] = 3*self._sample_x[N]
        matrix[4*N-1,4*N-2] = 1

        # calculate variables
        self.variables = np.dot(np.linalg.inv(matrix),y)

    def interpolation(self, t):
        point_num = len(self._sample_x)
        for index, j in enumerate(self._sample_x):
            if t < j:
                index -= 1
                break
        if index == -1:
            index += 1
        elif index == point_num-1:
            index -= 1
        a3 = self.variables[4*index + 3]
        a2 =  self.variables[4*index + 2]
        a1 = self.variables[4*index + 1]
        a0 = self.variables[4*index + 0]
        result = a3*pow(t,3) + a2*pow(t,2) + a1*t + a0
        return result