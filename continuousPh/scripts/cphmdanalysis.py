from MDAnalysis.analysis.distances import dist
import math
import copy
from scipy.optimize import curve_fit

def henderson_hasselbalch_for_fitting(ph, hill, pka):
    return 1 / (1 + 10**(hill*(pka - ph)))

def compute_pkas(phs, lambda_files):
    data = []
    for i in range(len(lambda_files[0].s)):
        s_values = [lambda_file.s[i] for lambda_file in lambda_files]
        diffs = [abs(s_value - 0.5) for s_value in s_values]
        min_diff = min(diffs)
        if min_diff < 0.4:
            initial_params = [1, phs[diffs.index(min_diff)]]
            fit = curve_fit(henderson_hasselbalch_for_fitting, phs, s_values, p0 = initial_params)
            data.append([fit[0][0], math.sqrt(fit[1][0][0]), fit[0][1], math.sqrt(fit[1][1][1])])
        else:
            data.append([math.nan, math.nan, math.nan])
    return data

class lambda_file:
    
    def __init__(self, file_path = ''):
        with open(file_path, 'r') as lf:
            line = lf.readline()
            tokens = line.split()
            self.ntitr = int(tokens[-1])
            line = lf.readline()
            tokens = line.split()
            self.titrreses = [int(token) for j, token in enumerate(tokens)
                              if (j > 1 and token != tokens[j - 1])]
            self.titrcols = [j - 1 for j, token in enumerate(tokens)
                             if (j > 1 and token != tokens[j - 1])]
            line = lf.readline()
            tokens = line.split()
            self.spgrps = [int(token) for j, token in enumerate(tokens)
                           if (j > 1 and (int(token) == 0 or int(token) == 1
                                          or int(token) == 3))]
            self.steps = []    
            self.ititrs = [titrcol - 1 for titrcol in self.titrcols]
            line = lf.readline()
            self.lambdavals = [[] for j in range(self.ntitr)]
            for line in lf:
                tokens = line.split()
                self.steps.append(int(tokens[0]))
                for j in range(self.ntitr):
                    self.lambdavals[j].append(float(tokens[1 + j]))
                    
    def compute_s_values(self, smalllambda = 0.2, biglambda = 0.8, minstep = -1, maxstep = -1):
        self.running_s = [[] for titrres in self.titrreses]
        self.s = []
        self.mixed = []
        self.micros = [[] for titrres in self.titrreses]
        total_frames = len(self.lambdavals[0])
        index_min = 0
        index_max = len(self.lambdavals[0])
        if minstep > 0:
            for i, step in enumerate(self.steps):
                if step > minstep:
                    index_min = i
                    break
        if maxstep > 0:
            for i, step in enumerate(self.steps):
                if step > maxstep:
                    index_max = i
                    break
        for i, titrres in enumerate(self.titrreses):
            spgrp = self.spgrps[i]
            titrcol = self.titrcols[i] - 1
            number_protonated = 0
            number_unprotonated = 0
            if spgrp != 0:
                number_protonated_taut1 = 0
                number_unprotonated_taut1 = 0
                number_protonated_taut2 = 0
                number_unprotonated_taut2 = 0
            for j, mylambda in enumerate(self.lambdavals[titrcol]):
                if j >= index_min and j <= index_max:
                    if spgrp is 0:
                        if mylambda < smalllambda:
                            number_protonated += 1
                        elif mylambda > biglambda:
                            number_unprotonated += 1
                    else:
                        myx = self.lambdavals[titrcol + 1][j]
                        if mylambda < smalllambda and myx < smalllambda:
                            number_protonated_taut1 += 1
                            number_protonated += 1
                        elif mylambda < smalllambda and myx > biglambda:
                            number_protonated_taut2 += 1
                            number_protonated += 1
                        elif mylambda > biglambda and myx < smalllambda:
                            number_unprotonated_taut1 += 1
                            number_unprotonated += 1
                        elif mylambda > biglambda and myx > biglambda:
                            number_unprotonated_taut2 += 1
                            number_unprotonated += 1
                    if number_protonated + number_unprotonated > 0:
                        self.running_s[i].append(number_unprotonated / (number_protonated + number_unprotonated))
                    else:
                        self.running_s[i].append(0)
            self.mixed.append((total_frames - number_protonated - number_unprotonated) / total_frames)
            if number_protonated + number_unprotonated > 0:
                self.s.append(number_unprotonated / (number_protonated + number_unprotonated))
            else:
                self.s.append(0)
            if spgrp is 1:
                if (number_unprotonated_taut1 + number_protonated) > 0:
                    self.micros[i].append(number_unprotonated_taut1 / (number_unprotonated_taut1
                                                                       + number_protonated))
                else:
                    self.micros[i].append(1.0)
                if number_unprotonated_taut2 + number_protonated > 0:
                    self.micros[i].append(number_unprotonated_taut2 / (number_unprotonated_taut2 +
                                                                      number_protonated))
                else:
                    self.micros[i].append(1.0)
            elif spgrp == 3:
                if number_unprotonated + number_protonated_taut1 > 0:
                    self.micros[i].append(number_unprotonated / (number_unprotonated +
                                                                number_protonated_taut1))
                else:
                    self.micros[i].append(1.0)
                if number_unprotonated + number_protonated_taut2 > 0:
                    self.micros[i].append(number_unprotonated / (number_unprotonated +
                                                                number_protonated_taut2))
                else:
                    self.micros[i].append(1.0)

# Here, I am just reading in in the pH values. May want to read
#other information in the future.
class rem_file:
    def __init__(self, file_path, swap_freq = 500):
        self.swap_freq = swap_freq
        self.phs = []
        with open(file_path, 'r') as rem:
            while True:
                line = rem.readline()
                if 'exchange' in line:
                    break
            temp_phs = []
            while True:
                line = rem.readline()
                if 'exchange' in line:
                    break
                temp_phs.append(float(line.split()[3]))
            self.phs.append(temp_phs)
            while True:
                temp_phs = []
                for i in range(len(self.phs[0])):
                    line = rem.readline()
                    if line == '':
                        break
                    temp_phs.append(float(line.split()[3]))
                line = rem.readline()
                if line == '':
                    break
                self.phs.append(temp_phs)

#This function assumes that the lambda files and remfile came from a normal replica-exchange simulation
def wrap_lambda_files(remfile, lambda_file_array, ph_array, smalllambda = 0.2, biglambda = 0.8, minstep = -1,
                      maxstep = -1):
    wrapped_lambda_files = copy.deepcopy(lambda_file_array)
    for wrapped_lambda_file in wrapped_lambda_files:
        wrapped_lambda_file.lambdavals = [[] for i in range(len(lambda_file_array[0].lambdavals))]
        wrapped_lambda_file.steps = []
    for i, step in enumerate(lambda_file_array[0].steps):
        if step < len(remfile.phs) * remfile.swap_freq:
            rem_index = math.floor(step / remfile.swap_freq)
            for j, wrapped_lambda_file in enumerate(wrapped_lambda_files):
                wrapped_lambda_file.steps.append(step)
                rep_index = remfile.phs[rem_index].index(ph_array[j])
                for k in range(len(wrapped_lambda_file.lambdavals)):
                    wrapped_lambda_file.lambdavals[k].append(lambda_file_array[rep_index].lambdavals[k][i])
    for wrapped_lambda_file in wrapped_lambda_files:
        wrapped_lambda_file.compute_s_values(smalllambda, biglambda, minstep, maxstep)
    return wrapped_lambda_files


#Dangerous - a work in progress
def hbond_calculations(universe, lambda_data, resid, nc_freq, lambda_freq, min_step = -1, max_step = -1, hbond_cutoff = 2.4, salt_bridge_cutoff = 4.0,
                       biglambda = 0.8, smalllambda = 0.2):
    if salt_bridge_cutoff > hbond_cutoff:
        cutoff = salt_bridge_cutoff
    else:
        cutoff = hbond_cutoff
    main_resname = universe.select_atoms('resid {}'.format(resid))[0].resname
    if main_resname == 'AS2' or main_resname == 'GL2':
        titrres = lambda_data.titrreses.index(int(resid))
        lambda_index = lambda_data.ititrs[titrres]
        x_index = lambda_index + 1
    donated_hbonds = 0
    accepted_hbonds = 0
    salt_bridges = 0
    mixed_frames = 0
    step_ratio = nc_freq / lambda_freq
    if min_step > 0:
        min_iter = math.floor(min_step / nc_freq)
    else:
        min_iter = 0
    if max_step > 0:
        max_iter = math.floor(max_step / nc_freq)
    else:
        max_iter = len(universe.trajectory)
    frames = max_iter - min_iter
    for i, ts in enumerate(universe.trajectory[min_iter:max_iter]):
        lambda_step_index = math.floor(i * step_ratio)
        if main_resname == 'AS2' or main_resname == 'GL2':
            if lambda_step_index >= len(lambda_data.lambdavals[lambda_index]):
                print('WTH {} {} {} {} {} {}'.format(step_ratio, max_iter, i, len(lambda_data.lambdavals[lambda_index]), lambda_step_index, len(universe.trajectory)))
            mylambda = lambda_data.lambdavals[lambda_index][lambda_step_index]
            myx = lambda_data.lambdavals[x_index][lambda_step_index]
            mixed = False
            if mylambda > biglambda and myx < smalllambda:
                protonated = False
                taut_state = 1
            elif mylambda > biglambda and myx > biglambda:
                protonated = False
                taut_state = 2
            elif mylambda < smalllambda and myx < smalllambda:
                protonated = True
                taut_state = 1
            elif mylambda < smalllambda and myx > biglambda:
                protonated = True
                taut_state = 2
            else:
                mixed = True
                mixed_frames += 1
        if not mixed:
            nearby_atoms = universe.select_atoms('around {} resid {}'.format(cutoff, resid))
            nearby_residues = []
            nearby_resnames = []
            for nearby_atom in nearby_atoms:
                nearby_resid = nearby_atom.resid
                resname = nearby_atom.resname
                if nearby_resid not in nearby_residues and nearby_resid != resid:
                    nearby_residues.append(nearby_resid)
                    nearby_resnames.append(resname)
            for j, nearby_resid in enumerate(nearby_residues):
                resname = nearby_resnames[j]
                if main_resname == 'AS2' or main_resname == 'GL2':
                    if main_resname == 'AS2':
                        od1_selection = universe.select_atoms('resid {} and name OD1'.format(resid))
                        od2_selection = universe.select_atoms('resid {} and name OD2'.format(resid))
                        hd1_selection = universe.select_atoms('resid {} and name HD1'.format(resid))
                        hd2_selection = universe.select_atoms('resid {} and name HD2'.format(resid))
                    else:
                        od1_selection = universe.select_atoms('resid {} and name OE1'.format(resid))
                        od2_selection = universe.select_atoms('resid {} and name OE2'.format(resid))
                        hd1_selection = universe.select_atoms('resid {} and name HE1'.format(resid))
                        hd2_selection = universe.select_atoms('resid {} and name HE2'.format(resid))
                    if not resname == 'PRO':
                        donor_selection = universe.select_atoms('resid {} and name H'.format(nearby_resid))
                        distances = [dist(od1_selection, donor_selection)[2][0],
                                     dist(od2_selection, donor_selection)[2][0]]
                        if not protonated:
                            if min(distances) < hbond_cutoff:
                                accepted_hbonds += 1
                        else:
                            if distances[0] < hbond_cutoff and taut_state == 2:
                                accepted_hbonds += 1
                            elif distances[1] < hbond_cutoff and taut_state == 1:
                                accepted_hbonds += 1
                    if protonated:
                        acceptor_selection = universe.select_atoms('resid {} and name O'.format(nearby_resid))
                        if taut_state == 1:
                            donor_selection = hd1_selection
                        else:
                            donor_selection = hd2_selection
                        distance = dist(donor_selection, acceptor_selection)[2][0]
                        if distance < hbond_cutoff:
                            donated_hbonds += 1
                    if not protonated and resname == 'ARG':
                        hh11_selection = universe.select_atoms('resid {} and name HH11'.format(nearby_resid))
                        hh12_selection = universe.select_atoms('resid {} and name HH12'.format(nearby_resid))
                        hh21_selection = universe.select_atoms('resid {} and name HH21'.format(nearby_resid))
                        hh22_selection = universe.select_atoms('resid {} and name HH22'.format(nearby_resid))
                        distance = min([dist(hh11_selection, od1_selection)[2][0],
                                        dist(hh12_selection, od1_selection)[2][0],
                                        dist(hh21_selection, od1_selection)[2][0],
                                        dist(hh22_selection, od1_selection)[2][0],
                                        dist(hh12_selection, od2_selection)[2][0],
                                        dist(hh12_selection, od2_selection)[2][0],
                                        dist(hh21_selection, od2_selection)[2][0],
                                        dist(hh22_selection, od2_selection)[2][0]])
                        if distance < salt_bridge_cutoff:
                            salt_bridges += 1
                    if protonated and resname == 'ARG' and taut_state == 1:
                        hh11_selection = universe.select_atoms('resid {} and name HH11'.format(nearby_resid))
                        hh12_selection = universe.select_atoms('resid {} and name HH12'.format(nearby_resid))
                        hh21_selection = universe.select_atoms('resid {} and name HH21'.format(nearby_resid))
                        hh22_selection = universe.select_atoms('resid {} and name HH22'.format(nearby_resid))
                        distance = min([dist(hh12_selection, od2_selection)[2][0],
                                        dist(hh12_selection, od2_selection)[2][0],
                                        dist(hh21_selection, od2_selection)[2][0],
                                        dist(hh22_selection, od2_selection)[2][0]])
                        if distance < hbond_cutoff:
                            accepted_hbonds += 1
                    if protonated and resname == 'ARG' and taut_state == 2:
                        hh11_selection = universe.select_atoms('resid {} and name HH11'.format(nearby_resid))
                        hh12_selection = universe.select_atoms('resid {} and name HH12'.format(nearby_resid))
                        hh21_selection = universe.select_atoms('resid {} and name HH21'.format(nearby_resid))
                        hh22_selection = universe.select_atoms('resid {} and name HH22'.format(nearby_resid))
                        distance = min([dist(hh12_selection, od1_selection)[2][0],
                                        dist(hh12_selection, od1_selection)[2][0],
                                        dist(hh21_selection, od1_selection)[2][0],
                                        dist(hh22_selection, od1_selection)[2][0]])
                        if distance < hbond_cutoff:
                            accepted_hbonds += 1
                    #Add non-titrating HIS?
                    if resname == 'HIP':
                        hip_titrres = lambda_data.titrreses.index(int(nearby_resid))
                        hip_lambda_index = lambda_data.ititrs[hip_titrres]
                        hip_lambda = lambda_data.lambdavals[hip_lambda_index][lambda_step_index]
                        hip_x = lambda_data.lambdavals[hip_lambda_index + 1][lambda_step_index]
                        other_hd1_selection = universe.select_atoms('resid {} and name HD1'.format(nearby_resid))
                        other_he2_selection = universe.select_atoms('resid {} and name HE2'.format(nearby_resid))
                        if hip_lambda > biglambda:
                            if hip_x < smalllambda:
                                if protonated and taut_state == 1:
                                    distance = dist(other_hd1_selection, od2_selection)[2][0]
                                    if distance < hbond_cutoff:
                                        accepted_hbonds += 1
                                    distance = dist(other_he2_selection, od1_selection)[2][0]
                                    if distance < hbond_cutoff:
                                        donated_hbonds += 1
                                elif protonated and taut_state == 2:
                                    distance = dist(other_hd1_selection, od1_selection)[2][0]
                                    if distance < hbond_cutoff:
                                        accepted_hbonds += 1
                                    distance = dist(other_he2_selection, od2_selection)[2][0]
                                    if distance < hbond_cutoff:
                                        donated_hbonds += 1
                                else:
                                    distance = min([dist(other_hd1_selection, od1_selection)[2][0],
                                                    dist(other_hd1_selection, od2_selection)[2][0]])
                                    if distance < hbond_cutoff:
                                        accepted_hbonds += 1
                            # Here I am assuming that if the H is not on ND1 it is on NE2. Change?
                            else:                                
                                if protonated and taut_state == 1:
                                    distance = dist(other_he2_selection, od2_selection)[2][0]
                                    if distance < hbond_cutoff:
                                        accepted_hbonds += 1
                                elif protonated and taut_state == 2:
                                    distance = dist(other_he2_selection, od1_selection)[2][0]
                                    if distance < hbond_cutoff:
                                        accepted_hbonds += 1
                                else:
                                    distance = min([dist(other_he2_selection, od1_selection)[2][0],
                                                    dist(other_he2_selection, od2_selection)[2][0]])
                                    if distance < hbond_cutoff:
                                        accepted_hbonds += 1
                        #Here I am assuming that if hip is not unprotonated, it can form a salt bridge. Change?
                        elif not protonated:
                            distance = min([dist(other_hd1_selection, od1_selection)[2][0],
                                            dist(other_he2_selection, od1_selection)[2][0],
                                            dist(other_hd1_selection, od2_selection)[2][0],
                                            dist(other_he2_selection, od2_selection)[2][0]])
                            if distance < salt_bridge_cutoff:
                                salt_bridges += 1
                        elif taut_state == 1:
                            distance = min([dist(other_hd1_selection, od2_selection)[2][0],
                                            dist(other_he2_selection, od2_selection)[2][0]])
                            if distance < hbond_cutoff:
                                accepted_hbonds += 1
                        else:
                            distance = min([dist(other_hd1_selection, od2_selection)[2][0],
                                            dist(other_he2_selection, od2_selection)[2][0]])
                            if distance < hbond_cutoff:
                                accepted_hbonds += 1
                    if not protonated and resname == 'LYS':
                        hz1_selection = universe.select_atoms('resid {} and name HZ1'.format(nearby_resid))
                        hz2_selection = universe.select_atoms('resid {} and name HZ2'.format(nearby_resid))
                        hz3_selection = universe.select_atoms('resid {} and name HZ3'.format(nearby_resid))
                        distance = min([dist(hz1_selection, od1_selection)[2][0],
                                        dist(hz2_selection, od1_selection)[2][0],
                                        dist(hz3_selection, od1_selection)[2][0],
                                        dist(hz1_selection, od2_selection)[2][0],
                                        dist(hz2_selection, od2_selection)[2][0],
                                        dist(hz3_selection, od2_selection)[2][0]])
                        if distance < salt_bridge_cutoff:
                            salt_bridges += 1
                    if protonated and resname == 'LYS' and taut_state == 1:
                        hz1_selection = universe.select_atoms('resid {} and name HZ1'.format(nearby_resid))
                        hz2_selection = universe.select_atoms('resid {} and name HZ2'.format(nearby_resid))
                        hz3_selection = universe.select_atoms('resid {} and name HZ3'.format(nearby_resid))
                        distance = min([dist(hz1_selection, od2_selection)[2][0],
                                        dist(hz2_selection, od2_selection)[2][0],
                                        dist(hz3_selection, od2_selection)[2][0]])
                        if distance < hbond_cutoff:
                            accepted_hbonds += 1
                    if protonated and resname == 'LYS' and taut_state == 2:
                        hz1_selection = universe.select_atoms('resid {} and name HZ1'.format(nearby_resid))
                        hz2_selection = universe.select_atoms('resid {} and name HZ2'.format(nearby_resid))
                        hz3_selection = universe.select_atoms('resid {} and name HZ3'.format(nearby_resid))
                        distance = min([dist(hz1_selection, od1_selection)[2][0],
                                        dist(hz2_selection, od1_selection)[2][0],
                                        dist(hz3_selection, od1_selection)[2][0]])
                        if distance < hbond_cutoff:
                            accepted_hbonds += 1
                    if resname == 'AS2' or resname == 'GL2':
                        other_titrres = lambda_data.titrreses.index(int(nearby_resid))
                        other_lambda_index = lambda_data.ititrs[other_titrres]
                        other_lambda = lambda_data.lambdavals[other_lambda_index][lambda_step_index]
                        other_x = lambda_data.lambdavals[other_lambda_index + 1][lambda_step_index]
                        if resname == 'AS2':
                            other_od1_selection = universe.select_atoms('resid {} and name OD1'.format(nearby_resid))
                            other_od2_selection = universe.select_atoms('resid {} and name OD2'.format(nearby_resid))
                            other_hd1_selection = universe.select_atoms('resid {} and name HD1'.format(nearby_resid))
                            other_hd2_selection = universe.select_atoms('resid {} and name HD2'.format(nearby_resid))
                        else:
                            other_od1_selection = universe.select_atoms('resid {} and name OE1'.format(nearby_resid))
                            other_od2_selection = universe.select_atoms('resid {} and name OE2'.format(nearby_resid))
                            other_hd1_selection = universe.select_atoms('resid {} and name HE1'.format(nearby_resid))
                            other_hd2_selection = universe.select_atoms('resid {} and name HE2'.format(nearby_resid))
                        if other_lambda > biglambda and protonated:
                            if taut_state == 1:
                                distance = min([dist(hd1_selection, other_od1_selection)[2][0],
                                                dist(hd1_selection, other_od2_selection)[2][0]])
                                if distance < hbond_cutoff:
                                    donated_hbonds += 1
                            else:
                                distance = min([dist(hd2_selection, other_od1_selection)[2][0],
                                                dist(hd2_selection, other_od2_selection)[2][0]])
                                if distance < hbond_cutoff:
                                    donated_hbonds += 1
                        #Again, assuming that if lambda < biglambda, then is protonated.
                        elif other_lambda < biglambda and not protonated:
                            if other_x < smalllambda:
                                distance = min([dist(other_hd1_selection, od1_selection)[2][0],
                                               dist(other_hd1_selection, od2_selection)[2][0]])
                                if distance < hbond_cutoff:
                                    accepted_hbonds += 1
                            else:
                                distance = min([dist(other_hd2_selection, od1_selection)[2][0],
                                               dist(other_hd2_selection, od2_selection)[2][0]])
                                if distance < hbond_cutoff:
                                    accepted_hbonds += 1
                        elif other_lambda < biglambda and protonated:
                            if other_x < smalllambda and taut_state == 1:
                                distance = dist(other_hd1_selection, od2_selection)[2][0]
                                if distance < hbond_cutoff:
                                    accepted_hbonds += 1
                                distance = dist(hd1_selection, other_od2_selection)[2][0]
                                if distance < hbond_cutoff:
                                    donated_hbonds += 1
                            elif other_x < smalllambda and taut_state == 2:
                                distance = dist(other_hd1_selection, od1_selection)[2][0]
                                if distance < hbond_cutoff:
                                    accepted_hbonds += 1
                                distance = dist(hd2_selection, other_od2_selection)[2][0]
                                if distance < hbond_cutoff:
                                    donated_hbonds += 1
                            elif other_x > smalllambda and taut_state == 1:
                                distance = dist(other_hd2_selection, od2_selection)[2][0]
                                if distance < hbond_cutoff:
                                    accepted_hbonds += 1
                                distance = dist(hd1_selection, other_od1_selection)[2][0]
                                if distance < hbond_cutoff:
                                    donated_hbonds += 1
                            elif other_x > smalllambda and taut_state == 2:
                                distance = dist(other_hd2_selection, od1_selection)[2][0]
                                if distance < hbond_cutoff:
                                    accepted_hbonds += 1
                                distance = dist(hd2_selection, other_od1_selection)[2][0]
                                if distance < hbond_cutoff:
                                    donated_hbonds += 1
                    if resname == 'SER':
                        hg_selection = universe.select_atoms('resid {} and name HG'.format(nearby_resid))
                        if not protonated:
                            distance = min([dist(hg_selection, od1_selection)[2][0],
                                           dist(hg_selection, od2_selection)[2][0]])
                        elif taut_state == 1:
                            distance = dist(hg_selection, od2_selection)[2][0]
                        elif taut_state == 2:
                            distance = dist(hg_selection, od1_selection)[2][0]
                        if distance < hbond_cutoff:
                            accepted_hbonds += 1
                    if resname == 'THR':
                        hg_selection = universe.select_atoms('resid {} and name HG1'.format(nearby_resid))
                        if not protonated:
                            distance = min([dist(hg_selection, od1_selection)[2][0],
                                           dist(hg_selection, od2_selection)[2][0]])
                        elif taut_state == 1:
                            distance = dist(hg_selection, od2_selection)[2][0]
                        elif taut_state == 2:
                            distance = dist(hg_selection, od1_selection)[2][0]
                        if distance < hbond_cutoff:
                            accepted_hbonds += 1
                    if resname == 'ASN':
                        other_od1_selection = universe.select_atoms('resid {} and name OD1'.format(nearby_resid))
                        other_hd21_selection = universe.select_atoms('resid {} and name HD21'.format(nearby_resid))
                        other_hd22_selection = universe.select_atoms('resid {} and name HD22'.format(nearby_resid))
                        if not protonated:
                            distance = min([dist(other_hd21_selection, od1_selection)[2][0],
                                           dist(other_hd22_selection, od1_selection)[2][0],
                                           dist(other_hd21_selection, od2_selection)[2][0],
                                           dist(other_hd22_selection, od2_selection)[2][0]])
                            if distance < hbond_cutoff:
                                accepted_hbonds += 1
                        elif taut_state == 1:
                            distance = min([dist(other_hd21_selection, od2_selection)[2][0],
                                           dist(other_hd22_selection, od2_selection)[2][0]])
                            if distance < hbond_cutoff:
                                accepted_hbonds += 1
                            distance = dist(other_od1_selection, hd1_selection)[2][0]
                            if distance < hbond_cutoff:
                                donated_hbonds += 1
                        elif taut_state == 2:
                            distance = min([dist(other_hd21_selection, od1_selection)[2][0],
                                           dist(other_hd22_selection, od1_selection)[2][0]])
                            if distance < hbond_cutoff:
                                accepted_hbonds += 1
                            distance = dist(other_od1_selection, hd2_selection)[2][0]
                            if distance < hbond_cutoff:
                                donated_hbonds += 1
                    if resname == 'GLN':
                        other_od1_selection = universe.select_atoms('resid {} and name OE1'.format(nearby_resid))
                        other_hd21_selection = universe.select_atoms('resid {} and name HE21'.format(nearby_resid))
                        other_hd22_selection = universe.select_atoms('resid {} and name HE22'.format(nearby_resid))
                        if not protonated:
                            distance = min([dist(other_hd21_selection, od1_selection)[2][0],
                                           dist(other_hd22_selection, od1_selection)[2][0],
                                           dist(other_hd21_selection, od2_selection)[2][0],
                                           dist(other_hd22_selection, od2_selection)[2][0]])
                            if distance < hbond_cutoff:
                                accepted_hbonds += 1
                        elif taut_state == 1:
                            distance = min([dist(other_hd21_selection, od2_selection)[2][0],
                                           dist(other_hd22_selection, od2_selection)[2][0]])
                            if distance < hbond_cutoff:
                                accepted_hbonds += 1
                            distance = dist(other_od1_selection, hd1_selection)[2][0]
                            if distance < hbond_cutoff:
                                donated_hbonds += 1
                        elif taut_state == 2:
                            distance = min([dist(other_hd21_selection, od1_selection)[2][0],
                                           dist(other_hd22_selection, od1_selection)[2][0]])
                            if distance < hbond_cutoff:
                                accepted_hbonds += 1
                            distance = dist(other_od1_selection, hd2_selection)[2][0]
                            if distance < hbond_cutoff:
                                donated_hbonds += 1
                    #Here assuming CYS is not titrating. Fix this later.
                    if resname == 'CYS':
                        hg_selection = universe.select_atoms('resid {} and name HG'.format(nearby_resid))
                        if not protonated:
                            distance = min([dist(hg_selection, od1_selection)[2][0],
                                           dist(hg_selection, od2_selection)[2][0]])
                        elif taut_state == 1:
                            distance = dist(hg_selection, od2_selection)[2][0]
                        elif taut_state == 2:
                            distance = dist(hg_selection, od1_selection)[2][0]
                        if distance < hbond_cutoff:
                            accepted_hbonds += 1
                    if resname == 'PRO':
                        n_selection = universe.select_atoms('resid {} and name N'.format(nearby_resid))
                        if protonated:
                            if taut_state == 1:
                                h_selection = hd1_selection
                            else:
                                h_selection = hd2_selection
                            distance = dist(n_selection, h_selection)[2][0]
                            if distance < hbond_cutoff:
                                donated_hbonds += 1
                    #Here ignoring that TYR can also accept
                    if resname == 'TYR':
                        h_selection = universe.select_atoms('resid {} and name HH'.format(nearby_resid))
                        if not protonated:
                            distance = min([dist(h_selection, od1_selection)[2][0],
                                           dist(h_selection, od2_selection)[2][0]])
                        elif taut_state == 1:
                            distance = dist(h_selection, od2_selection)[2][0]
                        else:
                            distance = dist(h_selection, od1_selection)[2][0]
                        if distance < hbond_cutoff:
                            accepted_hbonds += 1
                    if resname == 'TRP':
                        h_selection = universe.select_atoms('resid {} and name HE1'.format(nearby_resid))
                        if not protonated:
                            distance = min([dist(h_selection, od1_selection)[2][0],
                                           dist(h_selection, od2_selection)[2][0]])
                        elif taut_state == 1:
                            distance = dist(h_selection, od2_selection)[2][0]
                        else:
                            distance = dist(h_selection, od1_selection)[2][0]
                        if distance < hbond_cutoff:
                            accepted_hbonds += 1
    return accepted_hbonds, donated_hbonds, salt_bridges, frames, mixed_frames
