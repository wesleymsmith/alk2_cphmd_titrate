
# Necessary function for VMD_pH.py

# Script written by Jason Swails ca. 2009
# You are free to choose an open source license (e.g., GPL, LGPL, BSD, etc.)
# This is not sellable, anyway.

def GetAtoms(resname, startatom, state):

# NOTE, the atom numbers returned are numbered for VMD, which means that the numbering
# starts from ZERO.  This is going to be one less than the numbering you'd get from
# AMBER

   if resname == 'AS4':
      if state == 0:
         return [startatom + 9, startatom + 12, startatom + 13, startatom + 14]
      elif state == 1:
         return [startatom + 12, startatom + 13, startatom + 14]
      elif state == 2:
         return [startatom + 9, startatom + 13, startatom + 14] 
      elif state == 3:
         return [startatom + 9, startatom + 12, startatom + 14] 
      elif state == 4:
         return [startatom + 9, startatom + 12, startatom + 13] 
   elif resname == 'GL4':
      if state == 0:
         return [startatom + 12, startatom + 15, startatom + 16, startatom + 17] 
      elif state == 1:
         return [startatom + 15, startatom + 16, startatom + 17] 
      elif state == 2:
         return [startatom + 12, startatom + 16, startatom + 17] 
      elif state == 3:
         return [startatom + 12, startatom + 15, startatom + 17] 
      elif state == 4:
         return [startatom + 12, startatom + 15, startatom + 16] 
   elif resname == 'HIP':
      if state == 0:
         return []
      elif state == 1:
         return [startatom + 12]
      elif state == 2:
         return [startatom + 8]
   elif resname == 'CYS':
      if state == 0:
         return []
      elif state == 1:
         return [startatom + 7]
   elif resname == 'LYS':
      if state == 0:
         return []
      elif state == 1:
         return [startatom + 16]
   elif resname == 'TYR':
      if state == 0:
         return []
      else:
         return [startatom + 13]
