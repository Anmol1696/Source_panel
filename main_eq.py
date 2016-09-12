from sympy  import Symbol
from math   import sin, cos, pi, atan2, atan
from math   import log as ln

import numpy as np

from source_panel   import get_panels
from plot_functions import plot_C_p_vs_x

def get_sum_lambda_S(m_j_matrix, c_j_matrix, panels):
    sum_lambda_S = 0

    for i in range(len(panels)):
        S_i = ((panels[i][0][0] - panels[i][1][0])**2 + (panels[i][0][1] - panels[i][1][1])**2)**0.5
        sum_lambda_S += (m_j_matrix[i]*S_i**2) + c_j_matrix[i]*S_i

    return sum_lambda_S

def get_matrixs(panels, alpha):
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
                I_vs_matrix_m[i][j] = 2*S_i
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


def get_I_ij_or_I_vs(panel_i, panel_j, alpha, which_I):
    phi_i = atan2((panel_i[0][1] - panel_i[1][1]), (panel_i[0][0] - panel_i[1][0]))
    phi_j = atan2((panel_j[0][1] - panel_j[1][1]), (panel_j[0][0] - panel_j[1][0]))

    x_i, y_i = (panel_i[0][0] + panel_i[1][0])/2, (panel_i[0][1] + panel_i[1][1])/2

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

        I_vs_cj = (C_s/2)*ln((S_j**2 + (2.0*A*S_j) + B)/B) + ((D_s - (A*C_s))/E)*(atan((S_j + A)/E) - atan(A / E))
        I_vs_mj = (C_s*S_j) + ((D - (2*A*C_s))/2)*ln((S_j**2 + (2*S_j*A) + B)/B) - ((B*C) + (A*(D_s - (2*A*C_s))))*(atan((S_j + A)/A) - atan(A/E))

        return I_vs_mj, I_vs_cj 
    else: 
        return 0

def form_I_matrix(I_matrix_m, I_matrix_c, S_i_matrix):
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
                    I_matrix[i][num_panels]  = -1
                else:
                    I_matrix[i][j + 1] = -1
    
    return I_matrix
    
def solve_for_lambda(I_matrix, beta_i_matrix, v_infi):
    b_matrix = map(lambda x: v_infi * cos(x) * (-1), beta_i_matrix) + [0]*len(beta_i_matrix)

    x_matrix = np.linalg.solve(map(lambda x: (x/(2.0*pi)), I_matrix), b_matrix)

    m_j_matrix = x_matrix[:len(beta_i_matrix)]
    c_j_matrix = x_matrix[len(beta_i_matrix):]

    return m_j_matrix, c_j_matrix

def get_C_l_and_C_d(a, b, C_pi_matrix, beta_i_matrix, panels):
    C_l = 0
    C_d = 0

    for i in range(len(panels)):
        S_i = ((panels[i][0][0] - panels[i][1][0])**2 + (panels[i][0][1] - panels[i][1][1])**2)**0.5
        C_l += (sin(beta_i_matrix[i]) * C_pi_matrix[i] * S_i) / (2*a)
        C_d += (cos(beta_i_matrix[i]) * C_pi_matrix[i] * S_i) / (2*b)

    return C_l, C_d

def get_V_i(I_vs_matrix_m, I_vs_matrix_c, m_j_matrix, c_j_matrix, beta_i_matrix, v_infi):
    len_lambdas = len(beta_i_matrix)
    V_i_infi_matrix  = [0 for _ in range(len_lambdas)]
    
    for i in range(len_lambdas):
        temp_sum = 0
        for j in range(len_lambdas):
            temp_sum += ((m_j_matrix[j]*I_vs_matrix_m[i][j]) + (c_j_matrix[j]*I_vs_matrix_c[i][j]))/(2*pi)
        V_i_infi_matrix[i] = (v_infi * sin(beta_i_matrix[i])) + temp_sum

    return V_i_infi_matrix

def get_C_pi(V_i_infi_matrix, v_infi):
    C_pi_matrix = map(lambda x: (1 - (x/v_infi)**2), V_i_infi_matrix)

    return C_pi_matrix

def main():
    alphas_list = [69, -21, 9]
    a = 20
    b = 10
    v_infi = 20
    
    panels = get_panels(a, b, num_panels = 24, plot = False)

    for alpha in alphas_list:
        print 'Alpha -> ', alpha
        I_matrix, beta_i_matrix, I_vs_matrix_m, I_vs_matrix_c = get_matrixs(panels, alpha)

        m_j_matrix, c_j_matrix = solve_for_lambda(I_matrix, beta_i_matrix, v_infi)
    
        sum_lambda_S = get_sum_lambda_S(m_j_matrix, c_j_matrix, panels)
    
        V_i_infi_matrix = get_V_i(I_vs_matrix_m, I_vs_matrix_c, m_j_matrix, c_j_matrix, beta_i_matrix, v_infi) 
     
        C_pi_matrix = get_C_pi(V_i_infi_matrix, v_infi)

        C_l, C_d = get_C_l_and_C_d(a, b, C_pi_matrix, beta_i_matrix, panels)
   
        print 'Sum of lambda -> ', sum_lambda_S
        print 'C_l -> ', C_l
        print 'C_d -> ', C_d

        plot_C_p_vs_x(C_pi_matrix, panels)
    
if __name__ == "__main__":
    main()
