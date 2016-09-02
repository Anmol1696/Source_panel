"""
    This is a script that will be used to solve the source panel method
    This is strictly for an ellipse
"""

from math import sin, cos, pi

import matplotlib.pyplot    as plt
import numpy                as np

def get_panels(a, b, num_panels, plot):
    """
        Taking linear distributions of the angles
    """
    n = num_panels/4

    panels = []
    plot_curve = [[], []]

    theta = (90.0/n) - ((n - 1)/2.0)
    print theta 

    x1, y1 = 20.0, 0.0
    for i in range(n):
        sum_theta = ((i + 1)/2.0)*((2.0*theta) + i)
        x2 = a*cos(sum_theta*pi/180)
        y2 = b*sin(sum_theta*pi/180)
        
        panels.append([(x2, y2), (x1, y1)])

        x1, y1 = x2, y2

    for j in range(3):
        temp_panel = []
    
        for i in range(n):
            if j == 0:
                temp_panel.append([(panels[i][0][0] * (-1), panels[i][0][1]), (panels[i][1][0]*(-1), panels[i][1][1])])
            elif j == 1:
                temp_panel.append([(panels[i][0][0] * (-1), panels[i][0][1] * (-1)), (panels[i][1][0] * (-1), panels[i][1][1] * (-1))])
            elif j == 2:
                temp_panel.append([(panels[i][0][0], panels[i][0][1] * (-1)), (panels[i][1][0], panels[i][1][1] * (-1))])

        if not j == 1:
            temp_panel.reverse()
            for panel in panels:
                panel.reverse()

        panels += temp_panel

    if plot:
        x_points = [panel[0][0] for panel in panels]
        y_points = [panel[0][1] for panel in panels]
        plt.scatter(x_points, y_points)
        #plt.plot(x_points, y_points, '-o')
        plt.show()

    return panels

if __name__ == '__main__':
    panels = get_panels(20, 10, num_panels = 52, plot = True)

    print panels
    print 'Len of panels -> ', len(panels)
