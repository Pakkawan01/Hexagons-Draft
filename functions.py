'Part1 Permutation Angle Set'
def my_permutation(list_of_angle):
  """
    Computes a set of unique permutations of 6 angles, 
    eliminating those that are equivalent under rotation or reversal.

    Args:
        list_of_angle (list or tuple of int): A sequence of exactly 6 angles 
            (e.g., [135, 105, 135, 135, 75, 135]). Each angle represents a position 
            in the original configuration.

    Returns:
        set of tuple: A set of unique permutations (as tuples), excluding permutations 
        that are rotations or reversed versions of each other.
  """
  from itertools import permutations
  set_of_list = set(permutations(list_of_angle,6))
  R = set()
  from collections import deque
  set_list = list(set_of_list)
  while set_list != list() :
    angle_i = set_list[0]
    ro_angle_i = deque(angle_i)
    re_angle_i = deque(tuple(reversed(angle_i)))
    try :
      set_list.remove(angle_i)
      set_list.remove(tuple(reversed(angle_i)))
    except ValueError :
      pass
    R.add(angle_i)
    for i in range(5) :
      ro_angle_i.rotate(-1)
      re_angle_i.rotate(-1)
      try :
        set_list.remove(tuple(ro_angle_i))
        set_list.remove(tuple(re_angle_i))
      except ValueError :
        pass
  return R

'Part2 Pair Permutation'
"แยก N,S class"
def n_class_or_s_class(angle_i) :
  """
    Classifies a 6-angle configuration based on its rotational and reflective symmetries.

    Given a sequence of 6 angles, this function determines its symmetry class 
    by examining all possible rotations and reflections. It returns the canonical 
    representative and a label indicating the symmetry group.

    Args:
        angle_i (list or tuple of int): A sequence of exactly 6 angles 
            (e.g., [135, 105, 135, 135, 75, 135]).

    Returns:
        tuple: A tuple containing:
            - representative_angle (tuple of int): A canonical form of the angle configuration.
            - class_label (str): One of the following class labels:
                * 'N'   - No symmetry (12 unique variants)
                * 'S3'  - 3 unique variants under rotation/reflection
                * 'S5'  - 2 unique variants under rotation/reflection
                * 'S2'  - 6 variants, but with a palindromic variant (reflected version equals original)
                * 'S4'  - 6 variants with 180-degree rotational symmetry
                * 'S1'  - 6 variants without specific symmetry patterns

    Note:
        The function uses rotation and reflection (reverse) to analyze equivalence.
        Only specific lengths (6 angles) are supported.

    Raises:
        ValueError: If input does not contain exactly 6 elements.
    """
  R = set()
  from collections import deque
  ro_angle_i = deque(angle_i)
  re_angle_i = deque(tuple(reversed(angle_i)))
  R.add(angle_i)
  R.add(tuple(reversed(angle_i)))
  for i in range(5) : 
    ro_angle_i.rotate(-1)
    re_angle_i.rotate(-1)
    R.add(tuple(ro_angle_i))
    R.add(tuple(re_angle_i))
  if angle_i == (120,120,120,120,120,120) :
    return angle_i, None
  if len(R) == 12 :
    return angle_i, 'N'
  elif len(R) == 3 :
    return angle_i, 'S3'
  elif len(R) == 2 :
    return angle_i, 'S5'
  elif len(R) == 6 :
    R_r = list()
    R_list = list(R)
    for i in R_list :
      if i == tuple(reversed(i)):
        R_r.append(i)
        R.remove(i)
    if len(R) == 4 :
      return R_r[0], 'S2'
    elif len(R) == 6 :
      list_R = list(R)
      ro_test = deque(list_R[0])
      for i in range(3):
        ro_test.rotate(-1)
      angle_test = tuple(ro_test)
      if list_R[0] == angle_test :
        return angle_i, 'S4'
      else :
        return angle_i, 'S1'

"การจับคู่ระหว่าง N,S"
def pair_i_class(angle_with_class1, angle_with_class2 = ((0), 'Dont Have')) :
  """
  Generate a list of angle pairs (angle_a, angle_b) for a given input with class.

  Args:
    angle_with_class1 (Tuple[int, str] or List[Tuple[int, str]]): First list of angle with class.
    angle_with_class2 (Tuple[int, str] or List[Tuple[int, str]]): Second list of angle with class, default: ((0), 'Dont Have'))

  Returns:
    List[Tuple[int, int]]: List of matching angle pairs.
  """
  from collections import deque
  list_angle = [angle_with_class1,angle_with_class2]
  order = {'N': 0, 'S1': 1, 'S2': 2, 'S3': 3, 'S4': 4, 'S5': 5, 'Dont Have':6}
  list_angle.sort(key=lambda x: order[x[1]])
  angle_1, class_1 = list_angle[0]
  angle_2, class_2 = list_angle[1]
  if class_1 == 'N' :
    ro_angle1 = deque(angle_1)
    re_angle1 = deque(tuple(reversed(angle_1)))
    dict_of_ro_re1 = dict()
    dict_of_ro_re1['0'] = angle_1
    dict_of_ro_re1["0'"] = tuple(reversed(angle_1))
    range_list_ro = ['1','2','3','4','5']
    range_list_re = ["1'","2'","3'","4'","5'"]
    for i in range_list_ro :
      ro_angle1.rotate(-1)
      dict_of_ro_re1[i] = tuple(ro_angle1)
    for i in range_list_re :
      re_angle1.rotate(-1)
      dict_of_ro_re1[i] = tuple(re_angle1)
    if class_2 == 'Dont Have' : #N
      pair_of_n = [
          (dict_of_ro_re1['0'],dict_of_ro_re1["0'"]),
          (dict_of_ro_re1['0'],dict_of_ro_re1["1'"]),
          (dict_of_ro_re1['0'],dict_of_ro_re1["2'"]),
          (dict_of_ro_re1['0'],dict_of_ro_re1["3'"]),
          (dict_of_ro_re1['0'],dict_of_ro_re1["4'"]),
          (dict_of_ro_re1['0'],dict_of_ro_re1["5'"]),
          (dict_of_ro_re1['0'],dict_of_ro_re1['1']),
          (dict_of_ro_re1['0'],dict_of_ro_re1['2']),
          (dict_of_ro_re1['0'],dict_of_ro_re1['3'])
      ]
      return pair_of_n
    elif class_2 == 'N' : #NN
      ro_angle2 = deque(angle_2)
      re_angle2 = deque(tuple(reversed(angle_2)))
      dict_of_ro_re2 = dict()
      dict_of_ro_re2['0'] = angle_2
      dict_of_ro_re2["0'"] = tuple(reversed(angle_2))
      range_list_ro = ['1','2','3','4','5']
      range_list_re = ["1'","2'","3'","4'","5'",]
      for i in range_list_ro :
        ro_angle2.rotate(-1)
        dict_of_ro_re2[i] = tuple(ro_angle2)
      for i in range_list_re :
        re_angle2.rotate(-1)
        dict_of_ro_re2[i] = tuple(re_angle2)
      pair_of_nn = [
          (dict_of_ro_re1['0'],dict_of_ro_re2['0']),
          (dict_of_ro_re1['0'],dict_of_ro_re2['1']),
          (dict_of_ro_re1['0'],dict_of_ro_re2['2']),
          (dict_of_ro_re1['0'],dict_of_ro_re2['3']),
          (dict_of_ro_re1['0'],dict_of_ro_re2['4']),
          (dict_of_ro_re1['0'],dict_of_ro_re2['5']),
          (dict_of_ro_re1["0'"],dict_of_ro_re2["0'"]),
          (dict_of_ro_re1["0'"],dict_of_ro_re2["1'"]),
          (dict_of_ro_re1["0'"],dict_of_ro_re2["2'"]),
          (dict_of_ro_re1["0'"],dict_of_ro_re2["3'"]),
          (dict_of_ro_re1["0'"],dict_of_ro_re2["4'"]),
          (dict_of_ro_re1["0'"],dict_of_ro_re2["5'"])
      ]
      return pair_of_nn
    elif class_2 == 'S1' : #NS1
      ro_angle2 = deque(angle_2)
      dict_of_ro2 = dict()
      dict_of_ro2['0'] = angle_2
      range_list_ro = ['1','2','3','4','5']
      for i in range_list_ro :
        ro_angle2.rotate(-1)
        dict_of_ro2[i] = tuple(ro_angle2)
      pair_of_ns1 = [
          (dict_of_ro_re1['0'],dict_of_ro2['0']),
          (dict_of_ro_re1['0'],dict_of_ro2['1']),
          (dict_of_ro_re1['0'],dict_of_ro2['2']),
          (dict_of_ro_re1['0'],dict_of_ro2['3']),
          (dict_of_ro_re1['0'],dict_of_ro2['4']),
          (dict_of_ro_re1['0'],dict_of_ro2['5'])
      ]
      return pair_of_ns1
    elif class_2 == 'S2' : #NS2
      ro_angle2 = deque(angle_2)
      dict_of_ro2 = dict()
      dict_of_ro2['0'] = angle_2
      range_list_ro = ['1','2','3','4','5']
      for i in range_list_ro :
        ro_angle2.rotate(-1)
        dict_of_ro2[i] = tuple(ro_angle2)
      pair_of_ns2 = [
          (dict_of_ro_re1['0'],dict_of_ro2['0']),
          (dict_of_ro_re1['0'],dict_of_ro2['1']),
          (dict_of_ro_re1['0'],dict_of_ro2['2']),
          (dict_of_ro_re1['0'],dict_of_ro2['3']),
          (dict_of_ro_re1['0'],dict_of_ro2['4']),
          (dict_of_ro_re1['0'],dict_of_ro2['5'])
      ]
      return pair_of_ns2
    elif class_2 == 'S4' : #NS4
      ro_angle2 = deque(angle_2)
      re_angle2 = deque(tuple(reversed(angle_2)))
      dict_of_ro_re2 = dict()
      dict_of_ro_re2['0'] = angle_2
      dict_of_ro_re2["0'"] = tuple(reversed(angle_2))
      range_list_ro = ['1','2']
      range_list_re = ["1'","2'"]
      for i in range_list_ro :
        ro_angle2.rotate(-1)
        dict_of_ro_re2[i] = tuple(ro_angle2)
      for i in range_list_re :
        re_angle2.rotate(-1)
        dict_of_ro_re2[i] = tuple(re_angle2)
      pair_of_ns4 = [
          (dict_of_ro_re1['0'],dict_of_ro_re2['0']),
          (dict_of_ro_re1['0'],dict_of_ro_re2['1']),
          (dict_of_ro_re1['0'],dict_of_ro_re2['2']),
          (dict_of_ro_re1["0'"],dict_of_ro_re2['0']),
          (dict_of_ro_re1["0'"],dict_of_ro_re2['1']),
          (dict_of_ro_re1["0'"],dict_of_ro_re2['2']),
      ]
      return pair_of_ns4
    elif class_2 == 'S5' : #NS5
      pair_of_ns5 = [(dict_of_ro_re1['0'],angle_2),(dict_of_ro_re1["0'"],angle_2)]
      return pair_of_ns5
  elif class_1 == 'S1' :
    ro_angle1 = deque(angle_1)
    dict_of_ro1 = dict()
    dict_of_ro1['0'] = angle_1
    range_list_ro = ['1','2','3','4','5']
    for i in range_list_ro :
      ro_angle1.rotate(-1)
      dict_of_ro1[i] = tuple(ro_angle1)
    if class_2 == 'Dont Have' : #S1
      pair_of_s1 = [
          (dict_of_ro1['0'],dict_of_ro1['1']),
          (dict_of_ro1['0'],dict_of_ro1['2']),
          (dict_of_ro1['0'],dict_of_ro1['3'])
      ]
      return pair_of_s1
    elif class_2 == 'S1' : #S1S1
      ro_angle2 = deque(angle_2)
      dict_of_ro2 = dict()
      dict_of_ro2['0'] = angle_2
      range_list_ro = ['1','2','3','4','5']
      for i in range_list_ro :
        ro_angle2.rotate(-1)
        dict_of_ro2[i] = tuple(ro_angle2)
      pair_of_s1s1 = [
          (dict_of_ro1['0'],dict_of_ro2['0']),
          (dict_of_ro1['0'],dict_of_ro2['1']),
          (dict_of_ro1['0'],dict_of_ro2['2']),
          (dict_of_ro1['0'],dict_of_ro2['3']),
      ]
      return pair_of_s1s1
    elif class_2 == 'S2' : #S1S2
      ro_angle2 = deque(angle_2)
      dict_of_ro2 = dict()
      dict_of_ro2['0'] = angle_2
      range_list_ro = ['1','2','3','4','5']
      for i in range_list_ro :
        ro_angle2.rotate(-1)
        dict_of_ro2[i] = tuple(ro_angle2)
      pair_of_s1s2 = [(dict_of_ro1['0'],dict_of_ro2['0']),(dict_of_ro1['0'],dict_of_ro2['1']),(dict_of_ro1['0'],dict_of_ro2['2'])]
      return pair_of_s1s2
    elif class_2 == 'S3' : #S1S3
      ro_angle2 = deque(angle_2)
      dict_of_ro2 = dict()
      dict_of_ro2['0'] = angle_2
      range_list_ro = ['1','2']
      for i in range_list_ro :
        ro_angle2.rotate(-1)
        dict_of_ro2[i] = tuple(ro_angle2)
      pair_of_s1s3 = [(dict_of_ro1['0'],dict_of_ro2['0']),(dict_of_ro1['1'],dict_of_ro2['0'])]
      return pair_of_s1s3
    elif class_2 == 'S4' : #S1S4
      ro_angle2 = deque(angle_2)
      re_angle2 = deque(tuple(reversed(angle_2)))
      dict_of_ro_re2 = dict()
      dict_of_ro_re2['0'] = angle_2
      dict_of_ro_re2["0'"] = tuple(reversed(angle_2))
      range_list_ro = ['1','2']
      range_list_re = ["1'","2'"]
      for i in range_list_ro :
        ro_angle2.rotate(-1)
        dict_of_ro_re2[i] = tuple(ro_angle2)
      for i in range_list_re :
        re_angle2.rotate(-1)
        dict_of_ro_re2[i] = tuple(re_angle2)
      pair_of_s1s4 = [
          (dict_of_ro1['0'],dict_of_ro_re2['0']),
          (dict_of_ro1['1'],dict_of_ro_re2['0']),
          (dict_of_ro1['2'],dict_of_ro_re2['0'])
      ]
      return pair_of_s1s4
    elif class_2 == 'S5' : #S1S5
      pair_of_s1s5 = [(dict_of_ro1['0'],angle_2),(dict_of_ro1['1'],angle_2)]
      return pair_of_s1s5
  elif class_1 == 'S2' :
    ro_angle1 = deque(angle_1)
    dict_of_ro1 = dict()
    dict_of_ro1['0'] = angle_1
    range_list_ro = ['1','2','3','4','5']
    for i in range_list_ro :
      ro_angle1.rotate(-1)
      dict_of_ro1[i] = tuple(ro_angle1)
    if class_2 == 'Dont Have' : #S2
      pair_of_s2 = [
          (dict_of_ro1['0'],dict_of_ro1['1']),
          (dict_of_ro1['0'],dict_of_ro1['2']),
          (dict_of_ro1['0'],dict_of_ro1['3'])
      ]
      return pair_of_s2
    elif class_2 == 'S2' : #S2S2
      ro_angle2 = deque(angle_2)
      dict_of_ro2 = dict()
      dict_of_ro2['0'] = angle_2
      range_list_ro = ['1','2','3','4','5']
      for i in range_list_ro :
        ro_angle2.rotate(-1)
        dict_of_ro2[i] = tuple(ro_angle2)
      pair_of_s2s2 = [
          (dict_of_ro1['0'],dict_of_ro2['0']),
          (dict_of_ro1['0'],dict_of_ro2['1']),
          (dict_of_ro1['0'],dict_of_ro2['2']),
          (dict_of_ro1['0'],dict_of_ro2['3'])
      ]
      return pair_of_s2s2
    elif class_2 == 'S3' : #S2S3
      ro_angle2 = deque(angle_2)
      dict_of_ro2 = dict()
      dict_of_ro2['0'] = angle_2
      range_list_ro = ['1','2']
      for i in range_list_ro :
        ro_angle2.rotate(-1)
        dict_of_ro2[i] = tuple(ro_angle2)
      pair_of_s2s3 = [
          (dict_of_ro1['0'],dict_of_ro2['0']),
          (dict_of_ro1['0'],dict_of_ro2['1'])
      ]
      return pair_of_s2s3
    elif class_2 == 'S4' : #S2S4
      ro_angle2 = deque(angle_2)
      re_angle2 = deque(tuple(reversed(angle_2)))
      dict_of_ro_re2 = dict()
      dict_of_ro_re2['0'] = angle_2
      dict_of_ro_re2["0'"] = tuple(reversed(angle_2))
      range_list_ro = ['1','2']
      range_list_re = ["1'","2'"]
      for i in range_list_ro :
        ro_angle2.rotate(-1)
        dict_of_ro_re2[i] = tuple(ro_angle2)
      for i in range_list_re :
        re_angle2.rotate(-1)
        dict_of_ro_re2[i] = tuple(re_angle2)
      pair_of_s2s4 = [
          (dict_of_ro1['0'],dict_of_ro_re2['0']),
          (dict_of_ro1['1'],dict_of_ro_re2['0']),
          (dict_of_ro1['2'],dict_of_ro_re2['0'])
      ]
      return pair_of_s2s4
  elif class_1 == 'S3' :
    ro_angle1 = deque(angle_1)
    dict_of_ro1 = dict()
    dict_of_ro1['0'] = angle_1
    range_list_ro = ['1','2']
    for i in range_list_ro :
      ro_angle1.rotate(-1)
      dict_of_ro1[i] = tuple(ro_angle1)
    if class_2 == 'Dont Have' : #S3
      pair_of_s3 = [(dict_of_ro1['0'],dict_of_ro1['1'])]
      return pair_of_s3
  elif class_1 == 'S4' :
    ro_angle1 = deque(angle_1)
    re_angle1 = deque(tuple(reversed(angle_1)))
    dict_of_ro_re1 = dict()
    dict_of_ro_re1['0'] = angle_1
    dict_of_ro_re1["0'"] = tuple(reversed(angle_1))
    range_list_ro = ['1','2']
    range_list_re = ["1'","2'"]
    for i in range_list_ro :
      ro_angle1.rotate(-1)
      dict_of_ro_re1[i] = tuple(ro_angle1)
    for i in range_list_re :
      re_angle1.rotate(-1)
      dict_of_ro_re1[i] = tuple(re_angle1)
    if class_2 == 'Dont Have' : #S4
      pair_of_s4 = [
          (dict_of_ro_re1['0'],dict_of_ro_re1['1']),
          (dict_of_ro_re1['0'],dict_of_ro_re1["0'"]),
          (dict_of_ro_re1['0'],dict_of_ro_re1["1'"]),
          (dict_of_ro_re1['0'],dict_of_ro_re1["2'"]),
      ]
      return pair_of_s4
  elif class_1 == 'S5' : #S5
    pair_of_s5 = [(angle_1,tuple(reversed(angle_1)))]
    return pair_of_s5

"การจับคู่ทั้งหมด"
def pair_permutation(set_angle_with_class) :
  """
Generate all angle sequence pairs based on transformation rules from a set of sequences with class labels.

This function takes a set of angle sequences, where each sequence is associated with a class ID defining
its transformation behavior. It generates pairs in two ways:

  1. For each single sequence, generate all transformed versions (e.g., via rotation or reflection)
   according to its class using `pair_i_class(angle_with_class)`.
  2. For all unique combinations of two different sequences, generate all possible pairings
   based on their class rules using `pair_i_class(seq1, seq2)`.

  Args:
    set_angle_with_class (Set[Tuple[List[int], int]]): 
        A set of tuples, where each tuple consists of:
        - A list of angles representing a shape or configuration.
        - An integer class ID indicating transformation rules.

  Returns:
    List[Tuple[List[int], List[int]]]: 
        A list of all angle sequence pairs produced from the input,
        including both within-sequence and between-sequence combinations,
        as determined by the transformation rules for each class.

  """

  list_angle_with_class = list(set_angle_with_class)
  pair1 = pair_i_class(list_angle_with_class[0])
  for i in list_angle_with_class[1:] :
    pair1 = pair1 + pair_i_class(i)
  from itertools import combinations
  combi = list(combinations(list_angle_with_class,2))
  pair2 = pair_i_class(combi[0][0],combi[0][1])
  for i in combi[1:] :
    pair2 = pair2 + pair_i_class(i[0],i[1])
  all_pair = pair1 + pair2
  return all_pair

"Part3 solve Ax=0"
def solve_homogeous(angle_1,angle_2,show=False) :
  """
Solve a homogeneous system of linear equations derived from two sequences of angles,
and return the range of values for the free variable(s) that satisfy a system of inequalities.

This function transforms each input list of angles into directional vectors using trigonometric
calculations, constructs a matrix from these vectors, and solves the corresponding homogeneous 
system using sympy's RREF (Reduced Row Echelon Form). If the null space of the matrix has 
dimension 2, it finds a symbolic solution set (interval) where a linear combination of the two
free columns satisfies a set of inequalities, including positivity.

  Args:
    angle_1 (List[int]): A list of angles (in degrees), representing the first path or structure.
    angle_2 (List[int]): A list of angles (in degrees), representing the second path or structure.
    show (bool, optional): If True, prints the RREF matrix and pivot columns for debugging purposes. Defaults to False.

  Returns:
    sympy.Set: A symbolic interval representing the valid values of the free variable `s`
               that satisfy all derived inequalities. Returns an empty set if the system has
               no solution or the null space does not have dimension 2.

  Note:
    - Each angle list is assumed to describe a chain of directions where the next direction
      is determined by reflecting the previous angle.
    - The solution is based on solving inequalities derived from linear combinations
      of the free vectors in the null space.
  """

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
  if dimension_null_space == 2 :
    free_var_col = {0,1,2,3,4,5} - set(pivot_columns)
    free_col1, free_col2 = free_var_col
    col_1 = rref_matrix[:, free_col1]
    col_2 = rref_matrix[:, free_col2]
    s = sp.symbols('s') 
    ineq_all = [(s * a + b < 0) for a, b in zip(col_1, col_2)] + [(s > 0)]
    solution = sp.S.Reals
    for ineq in ineq_all :
      sol = sp.solveset(ineq, s, domain=sp.S.Reals)
      try :
        solution = solution.intersect(sol)
        continue
      except Exception :
        return sp.EmptySet
    return solution
  else :
    return sp.EmptySet

def solve_homogeous_dim(angle_1,angle_2) :
  """
Compute the dimension of the null space of a homogeneous system constructed from two sequences of angles.

This function transforms two sequences of angles into directional vectors using trigonometric projections.
It then builds a matrix from these vectors, reduces it to row-reduced echelon form (RREF), and calculates 
the dimension of its null space, which indicates the number of free variables (degrees of freedom) 
in the system.

Args:
    angle_1 (List[int]): A list of angles (in degrees), representing the first angle sequence.
    angle_2 (List[int]): A list of angles (in degrees), representing the second angle sequence.

Returns:
    int: The dimension of the null space (i.e., number of free variables) of the homogeneous system,
         which ranges from 0 to 6.
        
  """

  import sympy as sp
  def angle_from_x_axis(angles_list):
    new_angles = [0]
    for i in range(0, len(angles_list) - 1):
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
  dimension_null_space = 6 - len(pivot_columns)
  return dimension_null_space

def solve_homogeous_3dim(angle_1,angle_2,show=False) :
  """
Solve a homogeneous system of linear equations in 3D based on two sequences of angles,
and compute the maximum feasible values of two free variables under linear inequality constraints.

This function constructs a 3D system by transforming angle sequences into directional vectors,
forms a matrix of these vectors, and reduces it to RREF. If the null space has dimension 3, 
it formulates a set of linear inequalities and uses linear programming to find the maximum 
feasible values of `x` and `y` under the derived constraints.

  Args:
    angle_1 (List[int]): A list of angles (in degrees), representing the first 3D path or structure.
    angle_2 (List[int]): A list of angles (in degrees), representing the second 3D path or structure.
    show (bool, optional): If True, displays the RREF matrix and constraint set for debugging. Defaults to False.

  Returns:
    Tuple[sympy.Expr, sympy.Expr] or sympy.EmptySet:
        - A tuple containing the maximum feasible values for variables `x` and `y` under the constraint system,
          as calculated by linear programming (lpmax from sympy).
        - Returns `EmptySet` if the system has no solution or the null space is not 3-dimensional.

  Raises:
    Exception: If the system results in any zero side length, which violates the problem's physical constraints.

  Note:
    - This function assumes a system of 6 equations in 6 variables, derived from cosine and sine values
      of the direction chains formed from input angles.
    - It uses symbolic linear programming (`lpmax`) to solve for optimal values of `x` and `y`
      constrained by geometric relationships.
  """

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