# -*- coding: utf-8 -*-

from matplotlib.widgets import Button
from matplotlib.backend_bases import MouseEvent
import matplotlib.pyplot as plt
import matplotlib.gridspec as gds
import numpy as np
import math as mt
import CubicSplineCalculator as csc
import pandas as pd

class CubicSplineCourseGenerator(object):
    def __init__(self):
        self._figure         = None
        self._axes_plot      = None
        self._axes_yaw       = None
        self._axes_button    = None
        self._axes_curvature = None
        self._plot_points    = None
        self._plot_line      = None
        self._plot_yaw       = None
        self._plot_curvature = None
        self._sample_x       = np.empty([0, 0])
        self._sample_y       = np.empty([0, 0])
        self._cubic_spline_x = np.empty([0, 0])
        self._cubic_spline_y = np.empty([0, 0])
        self._elap_time_s    = np.empty([0, 0])
        self._yaw_deg        = np.empty([0, 0])
        self._curvature      = np.empty([0, 0])
        self._dt_sec         = 0.05
        self._init_plot()
    
    def _init_plot(self):
        self._figure = plt.figure(figsize=(15, 6))
        gridPos = gds.GridSpec(2, 2)
        plt.subplots_adjust(wspace=0.4, hspace=0.7)
        # Generated course data
        axes_plot = plt.subplot(gridPos[:,0])
        axes_plot.set_xlim([-10, 90])
        axes_plot.set_ylim([-10, 90])
        axes_plot.set_xlabel('X [m]')
        axes_plot.set_ylabel('Y [m]')
        axes_plot.grid()
        self._axes_plot = axes_plot
        # Calculated Yaw Angle
        axes_yaw = plt.subplot(gridPos[0,1])
        axes_yaw.set_xlabel('Time [s]')
        axes_yaw.set_ylabel('Yaw Angle [deg]')
        axes_yaw.grid()
        self._axes_yaw = axes_yaw
        # Calculated Curvature
        axes_curvature = plt.subplot(gridPos[1,1])
        axes_curvature.set_xlabel('Time [s]')
        axes_curvature.set_ylabel('Curvature')
        axes_curvature.grid()
        self._axes_curvature = axes_curvature
        # Saving button as CSV
        axes_button = plt.axes([0.05, 0.7, 0.08, 0.05])
        button_to_csv = Button(axes_button, 'Save as CSV')
        button_to_csv.on_clicked(self._save_as_csv)
        self._figure.subplots_adjust(left=0.2, bottom=None, right=None, top=None)
        self._figure.canvas.mpl_connect('button_press_event', self._on_click)
        plt.show()

    def _find_neighber_point(self, event):
        distance_threshold = 3.0
        nearest_point      = None
        minimum_distance   = mt.sqrt(2*(100**2))
        for x, y in zip(self._sample_x, self._sample_y):
            distance = mt.hypot((event.xdata-x), (event.ydata-y))
            if distance < minimum_distance:
                minimum_distance = distance
                nearest_point    = (x, y)
        if minimum_distance < distance_threshold:
            return nearest_point
        return None

    def _add_point(self, event):
        if isinstance(event, MouseEvent):
            self._sample_x = np.append(self._sample_x, event.xdata)
            self._sample_y = np.append(self._sample_y, event.ydata)

    def _remove_point(self, point):
        if point[0] in self._sample_x and point[1] in self._sample_y:
            rem_idx_x = np.where(self._sample_x==point[0])
            self._sample_x = np.delete(self._sample_x, rem_idx_x)
            rem_idx_y = np.where(self._sample_y==point[1])
            self._sample_y = np.delete(self._sample_y, rem_idx_y)

    def _update_plot(self):
        if  self._sample_x.size == 0 and self._sample_y.size == 0:
            self._plot_points.set_data([], [])
            self._plot_line.set_data([], [])
        else:
            # add new plot
            if not self._plot_points:
                self._plot_points,    = self._axes_plot.plot(self._sample_x, self._sample_y, ".b", ms=10)
                self._plot_line,      = self._axes_plot.plot([], [], "b", ms=5)
                self._plot_yaw,       = self._axes_yaw.plot([], [], "b", ms=5)
                self._plot_curvature, = self._axes_curvature.plot([], [], "b", ms=5)
            # update current plot
            else:
                self._plot_points.set_data(self._sample_x, self._sample_y)
                # calculate cubic spline
                sample_point_index = np.arange(0, len(self._sample_x), 1, dtype=float)
                csc_obj_x = csc.CubicSplineCalculator(sample_point_index, self._sample_x)
                csc_obj_y = csc.CubicSplineCalculator(sample_point_index, self._sample_y)
                interp_index = np.arange(0, len(self._sample_x)-1, 0.01, dtype=float)
                self._cubic_spline_x = [csc_obj_x.interpolation(i) for i in interp_index]
                self._cubic_spline_y = [csc_obj_y.interpolation(i) for i in interp_index]
                self._elap_time_s    = [self._dt_sec*i for i in range(len(interp_index))]
                self._plot_line.set_data(self._cubic_spline_x, self._cubic_spline_y)
                # insert yaw
                self._calculate_yaw_angle()
                self._plot_yaw.set_data(self._elap_time_s, self._yaw_deg)
                self._axes_yaw.set_xlim([0, self._elap_time_s[-1]])
                self._axes_yaw.set_ylim([min(self._yaw_deg), max(self._yaw_deg)])
                # insert curvature
                self._calculate_curvature(csc_obj_x, csc_obj_y, interp_index)
                self._plot_curvature.set_data(self._elap_time_s, self._curvature)
                self._axes_curvature.set_xlim([0, self._elap_time_s[-1]])
                self._axes_curvature.set_ylim([min(self._curvature), max(self._curvature)])
        self._figure.canvas.draw()
    
    def _calculate_yaw_angle(self):
        yaw_index = range(len(self._cubic_spline_x))
        yaw_deg   = np.empty([0, 0])
        for i in yaw_index[1:-1]:
            x_forward  = self._cubic_spline_x[i+1]
            x_backward = self._cubic_spline_x[i-1]
            y_forward  = self._cubic_spline_y[i+1]
            y_backward = self._cubic_spline_y[i-1]
            yaw_rad    = mt.atan2((y_forward-y_backward), (x_forward-x_backward))
            yaw_deg = np.append(yaw_deg, np.rad2deg(yaw_rad))
        yaw_deg = np.insert(yaw_deg, 0, yaw_deg[0])
        yaw_deg = np.append(yaw_deg, yaw_deg[-1])
        self._yaw_deg = yaw_deg
    
    def _calculate_curvature(self, csc_obj_x, csc_obj_y, interp_index):
        curv_index = range(len(self._cubic_spline_x))
        curvature = np.empty([0, 0])
        for i in curv_index[1:-1]:
            delta_t = (i+1) - (i-1)
            d_x  = (self._cubic_spline_x[i+1]-self._cubic_spline_x[i-1])/delta_t
            dd_x = (self._cubic_spline_x[i+1]-2*self._cubic_spline_x[i]+self._cubic_spline_x[i-1])/pow(delta_t,2)
            d_y  = (self._cubic_spline_y[i+1]-self._cubic_spline_y[i-1])/delta_t
            dd_y = (self._cubic_spline_y[i+1]-2*self._cubic_spline_y[i]+self._cubic_spline_y[i-1])/pow(delta_t,2)
            k    = (dd_y*d_x - dd_x*d_y)/(pow(d_x,2) + pow(d_y,2))
            curvature = np.append(curvature, k)
        curvature = np.insert(curvature, 0, curvature[0])
        curvature = np.append(curvature, curvature[-1])
        self._curvature = curvature
    
    def _save_as_csv(self, event):
        columns_list = ['Time[s]','X[m]','Y[m]','Yaw[deg]','Curvature']
        data_frame = pd.DataFrame(index=[], columns=columns_list)
        data_frame['Time[s]']   = self._elap_time_s
        data_frame['X[m]']      = self._cubic_spline_x
        data_frame['Y[m]']      = self._cubic_spline_y
        data_frame['Yaw[deg]']  = self._yaw_deg
        data_frame['Curvature'] = self._curvature
        data_frame.to_csv('cubic_spline_course_data.csv',index=False)

    def _on_click(self, event):
        # left click
        if event.button == 1 and event.inaxes in [self._axes_plot]:
            point = self._find_neighber_point(event)
            if point:
                self._remove_point(point)
            else:
                self._add_point(event)
            self._update_plot()
        # right click
        elif event.button == 3 and event.inaxes in [self._axes_plot]:
            point = self._find_neighber_point(event)
            if point:
                self._remove_point(point)
                self._update_plot()

if __name__ == '__main__':

    course = CubicSplineCourseGenerator()