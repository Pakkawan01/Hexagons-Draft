from dimension_3 import solve_homogeous3
import pandas as pd
import sympy
import ast
sheet = ['case1','case2','case3','case6','case7','case8','case9','case10','case11']
for sheet_name_i in sheet :
    df = pd.read_excel('3-dim only.xlsx',sheet_name=sheet_name_i,header=None)
    print(f'stared {sheet_name_i}')
    solve_list = list()
    empty_list = list()
    Error_list = list()
    side_zeros = list()
    for i in range(len(df)) :
        angle = df.iloc[i,0]
        parsed_angle = ast.literal_eval(angle)
        angle1 = parsed_angle[0]
        angle2 = parsed_angle[1]
        try :
            solve = solve_homogeous3(angle1,angle2)
            solve_list.append((angle1,angle2,solve))
        except sympy.solvers.simplex.UnboundedLPError :
            Error_list.append((angle1,angle2))
        except sympy.solvers.simplex.InfeasibleLPError :
            empty_list.append((angle1,angle2))
        except Exception :
            side_zeros.append((angle1,angle2))
    pd.DataFrame(solve_list).to_excel('/Users/ramboair/My_Doc/Hexagons/new3dim/{}_solved.xlsx'.format(sheet_name_i),index=False)
    pd.DataFrame(empty_list).to_excel('/Users/ramboair/My_Doc/Hexagons/new3dim/{}_empty.xlsx'.format(sheet_name_i), index=False)
    pd.DataFrame(Error_list).to_excel('/Users/ramboair/My_Doc/Hexagons/new3dim/{}_unbounded.xlsx'.format(sheet_name_i), index=False)
    pd.DataFrame(side_zeros).to_excel('/Users/ramboair/My_Doc/Hexagons/new3dim/{}_sidezeros.xlsx'.format(sheet_name_i), index=False)
    print(f'finished {sheet_name_i}')