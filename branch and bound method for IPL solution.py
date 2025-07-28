# -*- coding: utf-8 -*-
"""
Created on Mon Jul 28 10:35:12 2025

@author: Lenovo
"""
#%%%%%%%%
#Basic Frame
#relaxation->branching->boundary->recursion

import numpy as np
from scipy.optimize import linprog
#objective function min,similarly for max only need to change c into -c
def branch_bound_min(c,A_ub,b_ub,bounds,int_vars=None):
    """
    Parameters
    ----------
    c : list
        Coefficient of bjective function (c^T x).
    A_ub, b_ub : list
        Coefficient of inequlity constraint matrix
    bounds : list
        Variable bounds, list of (low,high).
    int_vars : 
        Indices of integer variables (default: all variables)

    Returns:
        best.x: optimal interger solution
        best.fun: optimal objective function value

    """
    # Set default integer variables if not specified
    if int_vars==None:
        int_vars=list(range(len(c)))
    #Initialize optimal values (infinity for minimization similarly negative infinity for maximization)    
    best_fun = np.inf
    best_x = None
    epsilon = 1e-6
    
    def solve_relaxation(A,b):
        """Solve the linear programming relaxation (ignore integrality)"""
        res = linprog(c,A,b,bounds=bounds,method='highs')
        if res.success:
            return res.x, res.fun
        else:
            return None, np.inf
    def branch(A,b):
        """Recursive branching function to explore subproblems"""
        nonlocal best_fun,best_x
        
        #step 1: Solve the relaxation of the current subproblem
        x,fun = solve_relaxation(A,b)
        
        #Prune: No feasible solution or worse than the current best
        if x is None or fun>=best_fun-epsilon: # Add a 'Tolerance Range',making the algorithm more robust under the floating-point precision limitations of the computer.
            return
        #step 2: Check if relaxation solution is integer
        is_integer = True
        for i in int_vars:
            #Check if value is close to an integer(within tolerance)
            if not np.isclose(x[i],round(x[i]),atol=epsilon):
                is_integer = False
                break
        if is_integer:
            #update best solution if the current is better
            if fun<best_fun-epsilon:
                best_fun = fun
                best_x = np.round(x).astype(int)# Ensuse the value and the type of the solution is integer
            return
        
        #Step 3: Select first non_integer variable for branching
        branch_var = None
        for i in int_vars:
            if not np.isclose(x[i],round(x[i]),atol=epsilon):
                branch_var = i
                break    
        
        #Calculating the branching point
        x_val = x[branch_var]
        floor_val = int(np.floor(x_val))
        ceiling_val = floor_val + 1
        
        #Branch 1: Add comstraint x[branch_var]<=floor_val
        new_A1 = np.vstack([A,[0]*len(c)])
        new_A1[-1,branch_var] = 1
        new_b1 = np.hstack([b,floor_val])
        branch(new_A1,new_b1)
        
        #Branch 2:Add comstraint x[branch_var]>=ceiling_val
        new_A2 = np.vstack([A,[0]*len(c)])
        new_A2[-1,branch_var] = -1
        new_b2 = np.hstack([b,-ceiling_val])
        branch(new_A2,new_b2)
        
    # Starting branching from the oriinal problem
    branch(A_ub,b_ub)
    return best_x,best_fun
        


#c = [20,10]
#A_ub = [[5,4],[2,5]]
#b_ub = [24,13]
#bounds = [(0,4),(0,None)]
#bound_branching_min(c,A_ub,b_ub,bounds)
#%%%%%%%%
#PuLP
      
import pulp
model = pulp.LpProblem('model' ,pulp.LpMinimize)    
x1 = pulp.LpVariable('x1',lowBound=0,cat='Integer')       
x2 = pulp.LpVariable('x2',lowBound=0,cat='Integer')
model += -20*x1 - 10*x2 #objective function
model += 5*x1 +4*x2 <=24 #constraint
model += 2*x1 +5*x2 <=13
model.solve(pulp.PULP_CBC_CMD(msg=0))  #CBC solver  
print('optimal solution: x1={}, x2={}'.format(x1.varValue,x2.varValue))
print('objective function solution:{}'.format(model.objective.value()))
        

       
        
        
    
    
    
    
    
    
    
    
    
    