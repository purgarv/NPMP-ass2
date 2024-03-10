import numpy as np
import os

from hill_functions import *  

# MASTER-SLAVE D FLIP-FLOP MODEL
def ff_ode_model(Y, T, params): 
    a, not_a, q, not_q, d, clk = Y
    alpha1, alpha2, alpha3, alpha4, delta1, delta2, Kd, n = params

    da_dt     = alpha1*(pow(d/Kd, n)/(1 + pow(d/Kd, n) + pow(clk/Kd, n) + pow(d/Kd, n)*pow(clk/Kd, n))) + alpha2*(1/(1 + pow(not_a/Kd, n))) - delta1 *a 
    dnot_a_dt = alpha1*(1/(1 + pow(d/Kd, n) + pow(clk/Kd, n) + pow(d/Kd, n)*pow(clk/Kd, n))) + alpha2*(1/(1 + pow(a/Kd, n))) - delta1*not_a   
    dq_dt     = alpha3*((pow(a/Kd, n)*pow(clk/Kd, n))/(1 + pow(a/Kd, n) + pow(clk/Kd, n) + pow(a/Kd, n)*pow(clk/Kd, n))) + alpha4*(1/(1 + pow(not_q/Kd, n))) - delta2*q  
    dnot_q_dt = alpha3*((pow(not_a/Kd, n)*pow(clk/Kd, n))/(1 + pow(not_a/Kd, n) + pow(clk/Kd, n) + pow(not_a/Kd, n)*pow(clk/Kd, n))) + alpha4*(1/(1 + pow(q/Kd, n))) - delta2*not_q   

    return np.array([da_dt, dnot_a_dt, dq_dt, dnot_q_dt]) 


def simplify_and(n, Kd, *args):
    result = 1
    for e in args:
        result *= activate_1(e, Kd, n)
    return result


def xor_gate(R1, R2, Kd1, n1):
    return activate_1((pow(R1/Kd1, n1)+pow(R2/Kd1, n1))/(1 + pow(R1/Kd1, n1)*pow(R2/Kd1, n1)), Kd1, n1)




def counter_3bit_asynchronous(Y, T, params):
    a1, not_a1, q1, not_q1, a2, not_a2, q2, not_q2, a3, not_a3, q3, not_q3 = Y

    clk = get_clock(T)

    d1 = not_q1
    d2 = not_q2
    d3 = not_q3

    Y_FF1 = [a1, not_a1, q1, not_q1, d1, clk]
    Y_FF2 = [a2, not_a2, q2, not_q2, d2, not_q1]
    Y_FF3 = [a3, not_a3, q3, not_q3, d3, not_q2]

    dY1 = ff_ode_model(Y_FF1, T, params)
    dY2 = ff_ode_model(Y_FF2, T, params)
    dY3 = ff_ode_model(Y_FF3, T, params)

    dY = np.append(np.append(dY1, dY2), dY3)

    return dY

def counter_3bit_synchronous(Y, T, params):
    a1, not_a1, q1, not_q1, a2, not_a2, q2, not_q2, a3, not_a3, q3, not_q3 = Y
    alpha1, alpha2, alpha3, alpha4, delta1, delta2, Kd, n = params
    clk = get_clock(T)

    d1 = not_q1
    d2 = xor_gate(q1, q2, Kd1=Kd, n1=n)*100
    d3 = xor_gate(simplify_and(n, Kd, q1, q2)*100, q3, Kd1=Kd, n1=n)*100
    # print(q1>50, q2>50, q3>50, d3>50)

    Y_FF1 = [a1, not_a1, q1, not_q1, d1, clk]
    Y_FF2 = [a2, not_a2, q2, not_q2, d2, clk]
    Y_FF3 = [a3, not_a3, q3, not_q3, d3, clk]
    dY1 = ff_ode_model(Y_FF1, T, params)
    dY2 = ff_ode_model(Y_FF2, T, params)
    dY3 = ff_ode_model(Y_FF3, T, params)

    dY = np.append(np.append(dY1, dY2), dY3)

    return dY


def counter_3bit_synchronous_ce(Y, T, params):
    a1, not_a1, q1, not_q1, a2, not_a2, q2, not_q2, a3, not_a3, q3, not_q3, clock_enable = Y
    alpha1, alpha2, alpha3, alpha4, delta1, delta2, Kd, n = params
    clk = simplify_and(n, Kd, get_clock(T), clock_enable)*100

    d1 = not_q1
    d2 = xor_gate(q1, q2, Kd1=Kd, n1=n)*100
    d3 = xor_gate(simplify_and(n, Kd, q1, q2)*100, q3, Kd1=Kd, n1=n)*100
    # print(clk ,q1>50, q2>50, q3>50, d3>50)
    
    Y_FF1 = [a1, not_a1, q1, not_q1, d1, clk]
    Y_FF2 = [a2, not_a2, q2, not_q2, d2, clk]
    Y_FF3 = [a3, not_a3, q3, not_q3, d3, clk]
    dY1 = ff_ode_model(Y_FF1, T, params)
    dY2 = ff_ode_model(Y_FF2, T, params)
    dY3 = ff_ode_model(Y_FF3, T, params)

    dY = np.append(np.append(np.append(dY1, dY2), dY3), clock_enable)

    return dY