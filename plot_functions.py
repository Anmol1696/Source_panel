"""
    Script consisting of functions for ploting and others
"""

import matplotlib.pyplot    as plt
import numpy                as np

def plot_C_p_vs_x(C_pi_matrix, panels):
    """
        Plot the C_p vs the mid points of the panels
    """

    num_panels = len(panels)
    
    x_i_points = [(panel[0][0] + panel[1][0])/2.0 for panel in panels]
    print 'X_i points', x_i_points    
    
    plt.scatter(x_i_points[:num_panels/2], C_pi_matrix[:num_panels/2], color = 'red')
    plt.scatter(x_i_points[num_panels/2:], C_pi_matrix[num_panels/2:], color = 'blue')
    plt.show()
    
    return 0
