def solve_homogeous_3dim(angle_1,angle_2,show=False) :
  import sympy as sp
  def angle_from_x_axis(angles_list): 
    new_angles = [0]
    for i in range(0, len(angles_list)-1):
      next_angle = new_angles[i] + (180 - angles_list[i])
      new_angles.append(next_angle)
    return new_angles
  def calculate_trig(angles_list):
    cos_values = [sp.cos(sp.rad(angle)) for angle in angle_from_x_axis(angles_list)]
    sin_values = [sp.sin(sp.rad(angle)) for angle in angle_from_x_axis(angles_list)]
    return cos_values, sin_values
  cos_results1, sin_results1 = calculate_trig(angle_1)
  cos_results2, sin_results2 = calculate_trig(angle_2)
  A = sp.Matrix([
    cos_results1,
    sin_results1,
    sin_results2,
    cos_results2])
  rref_matrix, pivot_columns = A.rref()
  if show == True :
    print("RREF Matrix:")
    sp.pretty_print(rref_matrix)
    print("Pivot Columns:", pivot_columns)
  else :
    pass
  dimension_null_space = 6-len(pivot_columns)
  for i in range(4) :
    row_i = rref_matrix[i,:]
    if str(sorted(row_i)) == '[0, 0, 0, 0, 0, 1]' :
      raise Exception ("Side can't be zero!!")
    else : 
      continue
  if dimension_null_space == 3 :
    from sympy.solvers.simplex import lpmax
    from sympy.abc import x, y
    free_var_col = {0,1,2,3,4,5} - set(pivot_columns)
    free_col1, free_col2, free_col3 = free_var_col
    col_1 = rref_matrix[:, free_col1]
    col_2 = rref_matrix[:, free_col2]
    col_3 = rref_matrix[:, free_col3]
    ineq_all = [(x * round(a,20) + y*round(b,20)+round(c,20) <= 0) for a, b, c in zip(col_1, col_2, col_3)] + [(x >= 0),(y>=0)]
    # print(ineq_all)
    constraint = [c.evalf() if isinstance(c, sp.Expr) else c for c in ineq_all]
    if show == True :
      print('constraint:',constraint)
    else :
      pass
    solve_sim_x = lpmax(x, constraint)
    solve_sim_y = lpmax(y, constraint)
    return solve_sim_x,solve_sim_y
  else :
    return sp.EmptySet