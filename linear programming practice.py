# -*- coding: utf-8 -*-
"""
Created on Sat Jul 26 16:44:20 2025

@author: Lenovo
"""

from scipy.optimize import linprog
#objective function coefficients(min z)
c = [-1,-1,0,0,0,-1.65,-1.65,0,-2.3]
#Inequality constraints:A_ub*x = b_ub

A = [[5,0,0,0,0,10,0,0,0],# Machine A1 constraints     
     [0,7,0,0,0,0,9,0,12],# Machine A2 constraints 
     [0,0,6,0,0,0,0,8,0],# Machine B1 constraints 
     [0,0,0,4,0,0,0,0,11],# Machine B2 constraints 
     [0,0,0,0,7,0,0,0,0]]# Machine B3 constraints 
b = [6000,10000,4000,7000,4000,] #Maximal hour for each machine

#Equality constraints:A_eq*x = b_eq
A_eq = [[1,1,-1,-1,-1,0,0,0,0],# Coefficients for product I's process matching
        [0,0,0,0,0,1,1,-1,0]] # Coefficients for product II's process matching
b_eq = [0,0] # Right hand value of equality value

# Bounds: Each variable x1-x9 non_negative (no upper limit,so set to none)
bounds = [(0,None) for _ in range(9)]

#Solving the problem using the 'highs' solver(efficientfor linear programming)
result = linprog(
    c,
    A,
    b,
    A_eq,
    b_eq,
    bounds,
    method='highs')

if result.success:
    #Calculate maximum profit (recover the origial objective function)
    max_profit = -result.fun +1854 #Add back the consraint term
    # print results
    print('Solution successful!')
    print(f'The max profit is {max_profit:.2f}')
    print("Optimal production plan(variable values):")
    for i in range(9):
        var_name = f"x{i+1}"
        var_value = result.x[i]
        print(f'{var_name}={var_value:.2f}') # Optimal value for each variables
else:
    print("Solution failed.")
    print(f"Reason:{result.message}")


# Solution successful!
# The max profit is 4854.57
#  Optimal production plan(variable values):
# x1=1200.00
# x2=230.05
# x3=0.00
# x4=858.62
# x5=571.43
# x6=0.00
# x7=500.00
# x8=500.00
# x9=324.14













