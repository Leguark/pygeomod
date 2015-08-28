'''Create and analyse gravity forward models from geogrids using Noddy
Created on Feb 8, 2015

@author: flow
'''

import os, sys
import pickle
import numpy as np
import matplotlib.pyplot as plt


sys.path.append("/Users/flow/git/pygeomod/pygeomod")
sys.path.append("/Users/flow/git/pynoddy/pynoddy")

import geogrid
import subprocess

import pynoddy.history
import pynoddy.output



def determine_probabilities(model_dir):
    all_probs = {}
    i = 0
    for f in os.listdir(model_dir):
        if os.path.splitext(f)[1] == ".pkl" and "Sandstone" in f:
            #===================================================================
            # Load grid
            #===================================================================
            print("Work on file %s" % f)
            grid_file = open(os.path.join(model_dir, f), "r")
            grid_ori = pickle.load(grid_file)
            grid_file.close()
            #===================================================================
            # Extract subgrid
            #===================================================================
            subrange = (40,200,30,250,0,80)
            grid = grid_ori.extract_subgrid(subrange)
            #===================================================================
            # Initiate probability grids
            #===================================================================
            if i == 0:
                vals = grid.unit_ids
                for val in vals:
                    if not all_probs.has_key(val):
                        all_probs[val] = np.zeros_like(grid.grid, dtype = "float")
            i += 1
            #===================================================================
            # Add to probability grid
            #===================================================================
            for val in vals:
                all_probs[val] += (grid.grid == val)

    return all_probs, i


def analyse_geophysics(model, **kwds):
    """Simulate potential-fields and use for model analysis

    It is possible to directly define filter for processing of gravity

    **Arguments**:
        - *model_dir*: directory containing sub-directories with uncertainty runs
                    (uncertainty_run_01, etc.);

    **Optional keywords**:
        - *grav_min* = float : reject models with a grav value lower than this
        - *grav_max* = float : reject models with a grav value larger than this
    """
    # os.chdir(model_dir)
    all_gravs = []
    all_gravs_filtered = []
    all_mags = []
    all_mags_filtered = []
    all_probs = {}
    all_probs_filtered = {}
    i_all = 0
    i_filtered = 0
    used_grids = []
    used_grids_filtered = []
    f = model
    #for f in os.listdir(model_dir):
    #    if os.path.splitext(f)[1] == ".pkl" and "Sandstone" in f:
    #===================================================================
    # Load grid
    #===================================================================
#    grid_file = open(os.path.join(model_dir, f), "r")
#    grid_ori = pickle.load(grid_file)
#    grid_file.close()
    #===================================================================
    # Extract subgrid
    #===================================================================
#    subrange = (40,200,30,250,0,80)
#    grid = grid_ori.extract_subgrid(subrange)
    # grid = grid_ori
    # substitute 0 with something else in grid ids
    tmp_grid = np.zeros_like(grid.grid)
    tmp_grid[grid.grid == 0] += 1
    print tmp_grid.shape
    grid.set_basename(f.split(".")[0])
    print "Basename"
    print grid.basename
    grid.grid += tmp_grid
#             n_cells = np.prod(grid.grid.shape)
    grid.determine_geology_ids()
    #===================================================================
    # # set densities and magnetic susceptibilities
    #===================================================================
    densities = {0: 0.1,
                 1: 2610.,
                 2: 2920.,
                 3: 3100.,
                 4: 2920.,
                 5: 2610.}
    sus = {0: 0.001,
           1: 0.001,
           2: 0.001,
           3: 0.1,
           4: 0.001,
           5: 0.001}
    grid.set_densities(densities)
    grid.set_sus(sus)
    grid.write_noddy_files(gps_range = 0.0)
    print grid.unit_ids
    sim_type = "ANOM_FROM_BLOCK"
    history = grid.basename + "_base.his"
    output_name = grid.basename

    # save grid as vtk for testing:
    # grid_ori.export_to_vtk(vtk_filename = grid.basename)
    #===================================================================
    # Run gravity forward modeling
    #===================================================================
    out =  subprocess.Popen(['noddy.exe', history, output_name, sim_type],
                shell=False, stderr=subprocess.PIPE,
                stdout=subprocess.PIPE).stdout.read()

    #===================================================================
    # Initiate probability grids
    #===================================================================
    if i_all == 0:
        vals = grid.unit_ids
        for val in vals:
            if not all_probs.has_key(val):
                all_probs[val] = np.zeros_like(grid.grid, dtype = "float")
                all_probs_filtered[val] = np.zeros_like(grid.grid, dtype = "float")

    #===================================================================
    # Create plot and store data
    #===================================================================
    geophys = pynoddy.output.NoddyGeophysics(grid.basename)

#=================================================
#             Check gravity constraints
#=================================================
    filter_val = True
    if kwds.has_key("grav_max"):
        if np.max(geophys.grv_data) > kwds['grav_max']:
            filter_val = False
    if kwds.has_key("grav_min"):
        if np.min(geophys.grv_data) < kwds['grav_min']:
            filter_val = False
    if filter_val:
        all_gravs_filtered.append(geophys.grv_data)
        all_mags_filtered.append(geophys.mag_data)
        used_grids_filtered.append("%s/%s" % (model_dir, grid.basename))
#                 test_grid = np.zeros_like(grid.grid)
        for val in vals:
            all_probs_filtered[val] += (grid.grid == val)
#                     test_grid += grid.grid == val
        # check probability
#                 assert(np.sum(test_grid) == n_cells)
        i_filtered += 1
    all_gravs.append(geophys.grv_data)
    all_mags.append(geophys.mag_data)
    used_grids.append("%s/%s" % (model_dir, grid.basename))
#             test_grid = np.zeros_like(grid.grid)
    for val in vals:
        all_probs[val] += (grid.grid == val)
#                 test_grid += grid.grid == val
#             assert(np.sum(test_grid) == n_cells)
    i_all += 1

    #===================================================================
    # Export to vtk for test
    #===================================================================
#             grid_out = pynoddy.output.NoddyOutput(grid.basename)
#             grid_out.export_to_vtk(vtk_filename = grid.basename)


    #=======================================================================
    # Analyse statistics for all simulated grids
    #=======================================================================
    # all_gravs = np.array(all_gravs)
    return all_gravs, all_mags, used_grids, all_probs, i_all,\
           all_gravs_filtered, all_mags_filtered, used_grids_filtered, all_probs_filtered, i_filtered
#     f_all = open("all_gravs.pkl", 'w')
#     pickle.dump(all_gravs, f_all)
#     f_all.close()
#     return all_gravs


if __name__ == '__main__':
    analyse_geophys = True
    analyse_entropy = False
    all_gravs_exp = []
    all_mags_exp = []
    all_gravs_filtered_exp = []
    all_mags_filtered_exp = []
    used_grids_exp = []
    used_grids_filtered_exp = []
    i_all_exp = 0
    i_filtered_exp = 0
    os.chdir(r'/Users/flow/git/paper_sandstone/workdir/hres/')
    if analyse_geophys:
        i = 0
        print("Analyse gravity and magnetics")
        for f1 in os.listdir('.'):
            if "uncertainty_run" in f1:
                print("Analyse file %s" % f1)
                all_gravs, all_mags, used_grids, all_probs, i_all,\
                all_gravs_filtered, all_mags_filtered, used_grids_filtered,\
                all_probs_filtered, i_filtered = analyse_geophysics(f1, grav_max = 15.298)
#                 grav, mag, used_grids = analyse_geophysics(f1, grav_max = 17.)
                all_gravs_exp += all_gravs
                all_mags_exp += all_mags
                # now: add filtered results
                all_gravs_filtered_exp += all_gravs_filtered
                all_mags_filtered_exp += all_mags_filtered
                # add to used grids:
                used_grids_exp += used_grids
                used_grids_filtered_exp += used_grids_filtered
                # Finally: add probability grids
                if i == 0: # Initialise output dictionary
                    all_probs_exp = all_probs
                    all_probs_filtered_exp = all_probs_filtered
                else:
                    for key in all_probs.keys():
                        all_probs_exp[key] += all_probs[key]
                    for key in all_probs_filtered.keys():
                        all_probs_filtered_exp[key] += all_probs_filtered[key]
                # update counters
                i_all_exp += i_all
                i_filtered_exp += i_filtered
                i += 1


        #=======================================================================
        # Postprocess probabilities and calculate entropy
        #=======================================================================

        print("\n\tPostprocessing results...\n\n")

        # scale to probability
        for key in all_probs_exp.keys():
            print("n experiments: %d\n" % i_all_exp)
            print("n filtered experiments: %d\n" % i_filtered_exp)
            all_probs_exp[key] = all_probs_exp[key] / float(i_all_exp)
            all_probs_filtered_exp[key] = all_probs_filtered_exp[key] / float(i_filtered_exp)

        # calculate entropy
        h_exp = np.zeros_like(all_probs_exp[all_probs_exp.keys()[0]])
        nx, ny, nz = h_exp.shape
        for pg in all_probs_exp.values():
            for i in range(nx):
                for j in range(ny):
                    for k in range(nz):
                        if pg[i,j,k] != 0:
                            h_exp[i,j,k] -= pg[i,j,k] * np.log2(pg[i,j,k])

        h_filtered_exp = np.zeros_like(all_probs_filtered_exp[all_probs_filtered_exp.keys()[0]])
        nx, ny, nz = h_filtered_exp.shape
        for pg in all_probs_filtered_exp.values():
            for i in range(nx):
                for j in range(ny):
                    for k in range(nz):
                        if pg[i,j,k] != 0:
                            h_filtered_exp[i,j,k] -= pg[i,j,k] * np.log2(pg[i,j,k])

        print("\n\tPickle and save files...\n\n")

        #=======================================================================
        # Save all that to files
        #=======================================================================
        pickle.dump(all_probs_exp, open("all_probs_exp.pkl", 'w'))
        pickle.dump(all_probs_filtered_exp, open("all_probs_exp_filtered.pkl", 'w'))
        pickle.dump(all_gravs_exp, open("all_gravs_exp.pkl", 'w'))
        pickle.dump(all_gravs_filtered_exp, open("all_gravs_exp_filtered.pkl", 'w'))
        pickle.dump(all_mags_exp, open("all_mags_exp.pkl", 'w'))
        pickle.dump(all_mags_filtered_exp, open("all_mags_exp_filtered.pkl", 'w'))
        pickle.dump(h_exp, open("h_exp.pkl", 'w'))
        pickle.dump(h_filtered_exp, open("h_filtered_exp.pkl", 'w'))
        pickle.dump(used_grids_exp, open("used_grids_exp.pkl", 'w'))
        pickle.dump(used_grids_filtered_exp, open("used_grids_filtered_exp.pkl", 'w'))

#     #             break
#         f_all = open("all_gravs.pkl", 'w')
#         pickle.dump(all_gravs, f_all)
#         f_all.close()
#         f_mag = open("all_mags.pkl", 'w')
#         pickle.dump(all_mags, f_mag)
#         f_mag.close()
    if analyse_entropy:
        print("Analyse information entropy of grids")
        i = 0
        for f1 in os.listdir('.'):
            if "uncertainty_run" in f1:
                print("Analyse file %s" % f1)
                if i == 0:
                    all_probs, j = determine_probabilities(f1)
                else: # add to previous dict:
                    tmp_probs, j1 = determine_probabilities(f1)
                    j += j1 # for total count
                    for key in all_probs.keys():
                        all_probs[key] += tmp_probs[key]
                i += 1
#                 if i == 2: break
        # scale to probability
        for key in all_probs.keys():
            all_probs[key] = all_probs[key] / float(j)
        f_g = open("all_probs.pkl", "w")
        pickle.dump(all_probs, f_g)
        f_g.close()
