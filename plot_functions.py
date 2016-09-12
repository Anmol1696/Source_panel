import matplotlib.pyplot    as plt
import numpy                as np

def plot_C_p_vs_x(C_pi_matrix, panels):
    num_panels = len(panels)
    
    x_i_points = [(panel[0][0] + panel[1][0])/2.0 for panel in panels]
    
    plt.scatter(x_i_points[:num_panels/2], C_pi_matrix[:num_panels/2], color = 'green')
    plt.scatter(x_i_points[num_panels/2:], C_pi_matrix[num_panels/2:], color = 'red')
    plt.plot(x_i_points[:num_panels/2], C_pi_matrix[:num_panels/2], '-o')
    plt.plot(x_i_points[num_panels/2:], C_pi_matrix[num_panels/2:], '-o')
    plt.show()
    
    return 0
