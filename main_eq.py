"""
    This module will calculate the I_ij for the panels
"""

from sympy  import Symbol
from math   import sin, cos, pi, atan2, atan
from math   import log as ln

import numpy as np

from source_panel   import get_panels
from plot_functions import plot_C_p_vs_x

def V_infinity(time):
    return 5 * abs(cos(time))

def get_I_ij_or_I_vs(panel_i, panel_j, alpha, which_I):
    """
        This will contain the equation
        panel_i/j will have [(X_i+1, Y_i+1), (X_i, Y_i)] and same for j
        alpha should be in degrees
        I_vs is the integral used for calculating V_s for a given pannel
        which_I is either I_ij or I_vs
    """
    
    alpha = (alpha * pi) / 180.0

    phi_i = atan2((panel_i[0][1] - panel_i[1][1]), (panel_i[0][0] - panel_i[1][0]))
    phi_j = atan2((panel_j[0][1] - panel_j[1][1]), (panel_j[0][0] - panel_j[1][0]))

    x_i, y_i = (panel_i[0][0] - panel_i[1][0])/2, (panel_i[0][1] - panel_i[1][1])/2

    A   = -(x_i - panel_j[1][0])*cos(phi_j) - (y_i - panel_j[1][1])*sin(phi_j)
    B   = (x_i - panel_j[1][0])**2 + (y_i - panel_j[1][1])**2
    C   = sin(phi_i - phi_j)
    D   = (y_i - panel_j[1][1])*cos(phi_i) - (x_i - panel_j[1][0])*sin(phi_i)
    S_j = ((panel_j[0][0] - panel_j[1][0])**2 + (panel_j[0][1] - panel_j[1][1])**2)**0.5
    E   = (B - (A**2))**0.5

    if which_I == 'I_ij':
        I_ij_cj = (C/2)*ln((S_j**2 + (2.0*A*S_j) + B)/B) + ((D - (A*C))/E)*(atan((S_j + A)/E) - atan(A / E))
        I_ij_mj = (C*S_j) + ((D - (2*A*C))/2)*ln((S_j**2 + (2*S_j*A) + B)/B) - ((B*C) + (A*(D - (2*A*C))))*(atan((S_j + A)/A) - atan(A/E))
    
        return I_ij_mj, I_ij_cj, phi_i
    elif which_I == 'I_vs':
        C_s = -cos(phi_i - phi_j) 
        D_s = ((y_i - panel_j[1][1])*sin(phi_i)) + ((x_i - panel_j[1][0])*cos(phi_i))
        P_s = D_s - (2*A*C_s)
        Q_s = -1 * B * C_s

        #I_vs_cj = (((D - (A*C))/(2*E))*ln((S_j**2 + (2*A*S_j) + B)/B)) - (C * (atan((S_j + A)/E) - atan(A/E)))
        I_vs_cj = ((C_s/2)*ln(((S_j**2 + (2*A*S_j) + B)/B))) + (((D_s - (A*C_s))/E)*(atan((S_j + A)/E) - atan(A/E)))
        I_vs_mj = (C_s*S_j) + (P_s/2)*ln((S_j**2 + (2*A*S_j) + B)/B) + ((Q_s - (A*P_s))/E)*(atan((S_j + A)/E) - atan(A/E))
        
        return I_vs_mj, I_vs_cj 
    else:
        print '[ERROR] NO I mentioned from I_ij and I_vs'
        
        return 0

def get_matrixs(panels, alpha):
    """
        This will create a matrix for I and also the cos(beta_i) matrix
        alpha in degree
    """

    I_matrix_m      = [[0 for _ in range(len(panels))] for _ in range(len(panels))]
    I_matrix_c      = [[0 for _ in range(len(panels))] for _ in range(len(panels))]
    I_vs_matrix_m   = [[0 for _ in range(len(panels))] for _ in range(len(panels))]
    I_vs_matrix_c   = [[0 for _ in range(len(panels))] for _ in range(len(panels))]
    S_i_matrix      = [0 for _ in range(len(panels))]
    beta_i_matrix   = [0 for _ in range(len(panels))]

    for i in range(len(panels)):
        S_i = ((panels[i][0][0] - panels[i][1][0])**2 + (panels[i][0][1] - panels[i][1][1])**2)**0.5
        S_i_matrix[i] = S_i
        
        for j in range(len(panels)):
            if i == j:
                I_matrix_m[i][j]    = S_i*pi/2
                I_matrix_c[i][j]    = pi
                #I_vs_matrix_m[i][j] = 2*pi*S_i*(ln(S_i/2)-1)
                I_vs_matrix_m[i][j] = 0
                I_vs_matrix_c[i][j] = 0
            else:
                I_matrix_m[i][j], I_matrix_c[i][j], phi_i = get_I_ij_or_I_vs(panels[i], panels[j], alpha, 'I_ij')
                I_vs_matrix_m[i][j], I_vs_matrix_c[i][j] = get_I_ij_or_I_vs(panels[i], panels[j], alpha, 'I_vs')
                
                if beta_i_matrix[i] == 0:
                    beta_i_matrix[i] = (phi_i + (pi/2.0) - ((alpha*pi)/180.0))
                    if beta_i_matrix[i] < 0:
                        beta_i_matrix[i] += 2*pi

    I_matrix = form_I_matrix(I_matrix_m, I_matrix_c, S_i_matrix)

    return np.array(I_matrix), np.array(beta_i_matrix), np.array(I_vs_matrix_m), np.array(I_vs_matrix_c)

def form_I_matrix(I_matrix_m, I_matrix_c, S_i_matrix):
    """
        Given the 3 matrix this will form a 2n * 2n matrix
    """
    num_panels = len(S_i_matrix)

    I_matrix = [[0 for _ in range(2*num_panels)] for _ in range(2*num_panels)]

    for i in range(num_panels):
        for j in range(num_panels):
            I_matrix[i][j] = I_matrix_m[i][j]
        for j in range(num_panels, 2*num_panels):
            I_matrix[i][j] = I_matrix_c[i][j - num_panels]
    for i in range(num_panels, 2*num_panels):
        for j in range(num_panels):
            if (i - num_panels) == j:
                I_matrix[i][j] = S_i_matrix[i - num_panels]
            else:
                I_matrix[i][j] = 0
        for j in range(num_panels, 2*num_panels):
            if i == j:
                I_matrix[i][j]   = 1
                if j == (2*num_panels) - 1:
                    I_matrix[i][0]  = -1
                else:
                    I_matrix[i][j+1] = -1

    return I_matrix
    
def solve_for_lambda(I_matrix, beta_i_matrix, v_infi):
    """
        alpha is in degree
        here x = lambda_i/2*pi*V_infinity
    """
    
    #I_matrix, cos_matrix = get_matrixs(panels, alpha)
    print "I_matrix -> ", I_matrix
    print 'len I_matrix', len(I_matrix)
    print "Beta i matrix -> ", beta_i_matrix

    b_matrix = map(lambda x: v_infi * cos(x) * (-1), beta_i_matrix) + [0]*len(beta_i_matrix)

    x_matrix = np.linalg.solve(map(lambda x: (x/(2*pi)), I_matrix), b_matrix)

    m_j_matrix = x_matrix[:len(beta_i_matrix)]
    c_j_matrix = x_matrix[len(beta_i_matrix):]

    return m_j_matrix, c_j_matrix

def get_sum_lambda_S(m_j_matrix, c_j_matrix, panels):
    """
        Output will be Sigma(lambda_i*S_i)
    """
    sum_lambda_S = 0

    for i in range(len(panels)):
        S_i = ((panels[i][0][0] - panels[i][1][0])**2 + (panels[i][0][1] - panels[i][1][1])**2)**0.5
        sum_lambda_S += (m_j_matrix[i]*S_i**2) + c_j_matrix[i]*S_i

    return sum_lambda_S

def get_V_i(I_vs_matrix_m, I_vs_matrix_c, m_j_matrix, c_j_matrix, beta_i_matrix, v_infi):
    """
        lambda_v_infi is the result of solve_for lambda, lambda_i/2*pi*V_infinity
        V_i_matrix element is V_i/V_infi
    """
    len_lambdas = len(beta_i_matrix)
    V_i_infi_matrix  = [0 for _ in range(len_lambdas)]
    
    for i in range(len_lambdas):
        temp_sum = 0
        for j in range(len_lambdas):
            temp_sum += ((m_j_matrix[j]*I_vs_matrix_m[i][j]) + (c_j_matrix[j]*I_vs_matrix_c[i][j]))/(2*pi)
        V_i_infi_matrix[i] = (v_infi * sin(beta_i_matrix[i])) + temp_sum

    return V_i_infi_matrix

def get_C_pi(V_i_infi_matrix, v_infi):
    """
        This will calculate C_pi = 1 - (V_i/V_infinity)**2
    """
    C_pi_matrix = map(lambda x: (1 - (x/v_infi)**2), V_i_infi_matrix)

    return C_pi_matrix

def get_C_l_and_C_d(a, b, C_pi_matrix, beta_i_matrix, panels):
    """
        C_l and C_d are defined by C_l = sum(S_j * C_pi * sin(beta_i))/(2*a)
        C_d = sum(S_j * C_pi * cos(beta_i))/(2*b)
    """
    C_l = 0
    C_d = 0

    for i in range(len(panels)):
        S_i = ((panels[i][0][0] - panels[i][1][0])**2 + (panels[i][0][1] - panels[i][1][1])**2)**0.5
        C_l += (sin(beta_i_matrix[i]) * C_pi_matrix[i] * S_i) / (2*a)
        C_d += (cos(beta_i_matrix[i]) * C_pi_matrix[i] * S_i) / (2*b)

    return C_l, C_d

def main():
    """
        Main function
        alpha in degree
        Lambda_norm is lambda/2*pi*V_infinity
    """
    alpha = 0
    a = 20.0
    b = 10.0
    v_infi = 20

    panels = get_panels(a, b, num_panels = 24, plot = False)
    I_matrix, beta_i_matrix, I_vs_matrix_m, I_vs_matrix_c = get_matrixs(panels, alpha)

    m_j_matrix, c_j_matrix = solve_for_lambda(I_matrix, beta_i_matrix, v_infi)
    
    print 'm_j_matrix -> ', m_j_matrix
    print 'c_j_matrix -> ', c_j_matrix
    
    sum_lambda_S = get_sum_lambda_S(m_j_matrix, c_j_matrix, panels)
    print 'Sum of lambda -> ', sum_lambda_S
    
    V_i_infi_matrix = get_V_i(I_vs_matrix_m, I_vs_matrix_c, m_j_matrix, c_j_matrix, beta_i_matrix, v_infi) 
    
    print 'V_i_infi_matrix ->', V_i_infi_matrix
    
    C_pi_matrix = get_C_pi(V_i_infi_matrix, v_infi)

    C_l, C_d = get_C_l_and_C_d(a, b, C_pi_matrix, beta_i_matrix, panels)
   
    print 'Panels -> ', panels
    print 'I_vs_matrix_m -> ', I_vs_matrix_m
    print 'I_vs_matrix_c -> ', I_vs_matrix_c
    #print 'Values of lambda -> ',lambda_norm
    #print 'Len of x -> ', len(lambda_norm)
    print 'C_pi_matrix -> ', C_pi_matrix

    #sum_lambda_S = get_sum_lambda_S(lambda_norm, panels)

    #print 'Sum of lambda -> ', sum_lambda_S
    print 'C_l -> ', C_l
    print 'C_d -> ', C_d

    plot_C_p_vs_x(C_pi_matrix, panels)
    
if __name__ == "__main__":
    main()
