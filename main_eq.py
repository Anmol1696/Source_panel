"""
    This module will calculate the I_ij for the panels
"""

from sympy  import Symbol
from math   import sin, cos, pi, atan
from math   import log as ln

import numpy as np

from source_panel import get_panels

def V_infinity(time):
    return 5 * abs(cos(time))

def get_I_ij(panel_i, panel_j, alpha):
    """
        This will contain the equation
        panel_i/j will have [(X_i+1, Y_i+1), (X_i, Y_i)] and same for j
        alpha should be in degrees
    """
    
    alpha = (alpha * pi) / 180.0

    phi_i = atan((panel_i[0][1] - panel_i[1][1])/(panel_i[0][0] - panel_i[1][0]))
    phi_j = atan((panel_j[0][1] - panel_j[1][1])/(panel_j[0][0] - panel_j[1][0]))

    x_i, y_i = (panel_i[0][0] - panel_i[1][0])/2, (panel_i[0][1] - panel_i[1][1])/2

    A   = -(x_i - panel_j[1][0])*cos(phi_j - alpha) - (y_i - panel_j[1][1])*sin(phi_j - alpha)
    B   = (x_i - panel_j[1][0])**2 + (y_i - panel_j[1][1])**2
    C   = sin(phi_i - phi_j)
    D   = (y_i - panel_j[1][1])*cos(phi_i - alpha) - (x_i - panel_j[1][0])*sin(phi_i - alpha)
    S_j = ((panel_j[0][0] - panel_j[1][0])**2 + (panel_j[0][1] - panel_j[1][1])**2)**0.5
    E   = (B - A**2)**0.5

    I_ij = (C/2)*ln((S_j**2 + (2.0*A*S_j) + B)/B) + ((D - (A*C))/E)*(atan((S_j + A)/E) - atan(A/E))

    return I_ij, phi_i

def get_matrixs(panels, alpha):
    """
        This will create a matrix for I and also the cos(beta_i) matrix
        alpha in degree
    """

    I_matrix    = [[0 for _ in range(len(panels))] for _ in range(len(panels))]
    cos_matrix  = [0 for _ in range(len(panels))]

    for i in range(len(panels)):
        for j in range(len(panels)):
            if i == j:
                I_matrix[i][j] = pi
            else:
                I_ij, phi_i = get_I_ij(panels[i], panels[j], alpha)
                I_matrix[i][j] = I_ij
                if cos_matrix[i] == 0:
                    cos_matrix[i] = cos(phi_i + (pi/2.0) - (alpha*pi/180.0))


    return np.array(I_matrix), np.array(cos_matrix)

def solve_for_lambda(panels, alpha):
    """
        alpha is in degree
        here x = lambda_i/2*pi*V_infinity
    """
    
    I_matrix, cos_matrix = get_matrixs(panels, alpha)

    x = np.linalg.solve(I_matrix, cos_matrix)

    return x

def get_sum_lambda_S(lambdas, panels):
    """
        Output will be Sigma(lambda_i*S_i)
    """
    sum_lambda_S = 0

    for i in range(len(panels)):
        S_i = ((panels[i][0][0] - panels[i][1][0])**2 + (panels[i][0][1] - panels[i][1][1])**2)**0.5
        sum_lambda_S += lambdas[i] * S_i * 2 * pi * V_infinity(0)

    return sum_lambda_S

def main():
    """
        Main function
    """
    panels = get_panels(a = 20, b = 10, num_panels = 52, plot = False)

    x = solve_for_lambda(panels, alpha = 0)

    print x
    print 'Len of x -> ', len(x)

    sum_lambda_S = get_sum_lambda_S(x, panels)

    print sum_lambda_S

if __name__ == "__main__":
    main()
