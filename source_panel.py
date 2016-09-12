from math import sin, cos, pi

import matplotlib.pyplot as plt
import numpy as np

def get_panels(a, b, num_panels, plot):
    n = num_panels/4

    sub_panels  = []
    panels      = []
    plot_curve  = [[], []]

    theta = (90.0/n) - ((n - 1)/2.0)

    x1, y1 = 20.0, 0.0
    for i in range(n):
        sum_theta = ((i + 1)/2.0)*((2.0*theta) + i)
        x2 = round(a*cos(sum_theta*pi/180), 5)
        y2 = round(b*sin(sum_theta*pi/180), 5)
        
        sub_panels.append([(x2, y2), (x1, y1)])

        x1, y1 = x2, y2

    for j in range(4):
        for i in range(n):
            if j == 0:
                panels.append([(sub_panels[i][0][0]*(-1), sub_panels[i][0][1]), (sub_panels[i][1][0] * (-1), sub_panels[i][1][1])])
            elif j == 1:
                panels.append(sub_panels[::-1][i][::-1])
            elif j == 2:
                panels.append([(sub_panels[i][0][0], sub_panels[i][0][1] * (-1)), (sub_panels[i][1][0], sub_panels[i][1][1] * (-1))])
            elif j == 3:
                panels.append([(sub_panels[::-1][i][1][0] * (-1), sub_panels[::-1][i][1][1] * (-1)), (sub_panels[::-1][i][0][0] * (-1), sub_panels[::-1][i][0][1] * (-1))])

    if plot:
        x_points = [panel[0][0] for panel in panels]
        y_points = [panel[0][1] for panel in panels]
        plt.scatter(x_points, y_points)
        plt.plot(x_points, y_points, '-o')
        plt.show()

    return panels

if __name__ == '__main__':
    panels = get_panels(20, 10, num_panels = 24, plot = True)

    print panels
    print 'Len of panels -> ', len(panels)
