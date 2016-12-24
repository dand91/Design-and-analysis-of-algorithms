from __future__ import  division
from heapq import heappush as h_push
from heapq import heappop as h_pop
import math
import untangle
import scipy.signal
import copy
import matplotlib.pyplot as plt

#plt.ion()

from data_structures import *

count_correct = True

flatfrog_data = True

plot_best = True
plot_convergence = False
plot_averaged_dendo = True
plot_averaged_detrend_dendo = False
plot_averaged_internal = True
plot_averaged_detrend_internal = False
plot_averaged_external = True
plot_averaged_detrend_external = False

complete_list = True
single_list = False

def HAC(data,cluster_list,colors):

    # Initiate variables

    all_guess = 0
    correct = 0
    corr_id = 0
    tot_id = 0

    c_h = []

    N = len(data)

    id_counter = 0

    c_l = [None] * N

    iter = 0

    d_l = []

    d_l.append(0)

    a_i_l = []

    e_i_l = []

    h = []

    # Create new cluster for each data point

    for n in xrange(N):

        n_c = cl(id_counter)

        n_c.add(data[id_counter])

        c_l[n_c.id] = n_c

        id_counter += 1

    # Calculate initial distances between clusters

    for i in xrange(0,N):

        for j in xrange(i+1,N):

                d = c_l[i].c_dist(c_l[j])

                n_c_entry = c_e(d, c_l[i], c_l[j])

                h_push(h,n_c_entry)

    while iter < N-2:

        print "Iteration ", iter

        iter += 1

        # Find smallest value

        max_entry = h_pop(h)

        d_l.append(max_entry.dist)

        # Create new cluster

        n_c = cl(max_entry.ci.id)

        # Add points from previous clusters

        n_c.add(max_entry.ci.get_p_l())
        n_c.add(max_entry.cj.get_p_l())

        # Remove clusters who is to be merged

        c_l[max_entry.cj.id] = None
        c_l[max_entry.ci.id] = None

        # Remove old distances

        t_h = copy.deepcopy(h)

        h = [h_e for h_e in h if \
                 h_e.ci.id != max_entry.ci.id \
             and h_e.ci.id != max_entry.cj.id \
             and h_e.cj.id != max_entry.cj.id \
             and h_e.cj.id != max_entry.ci.id]

        # Calculate new distances

        for c in [x for x in c_l if x is not None]:

            k_d = []

            for h_e in t_h:

                h_e_l = ( h_e.ci.id, h_e.cj.id )

                if max_entry.ci.id in h_e_l and c.id in h_e_l:

                    k_d.append(h_e.dist)

                if max_entry.cj.id in h_e_l and c.id in h_e_l:

                    k_d.append(h_e.dist)

            assert len(k_d) == 2

            assert complete_list is not single_list

            if complete_list:

                d = 0.5 * (k_d[0] + k_d[1]) + 0.5 * math.fabs(k_d[0] - k_d[1])

            if single_list:

                d = 0.5 * (k_d[0] + k_d[1]) - 0.5 * math.fabs(k_d[0] - k_d[1])

            n_entry = c_e(d, n_c, c)

            h_push(h,n_entry)

        # Add new cluster

        c_l[n_c.id] = n_c

        # Save cluster

        history = [x for x in c_l if x != None]

        c_h.insert(iter-1, copy.deepcopy(history))

        # Calculate internal distances

        a_i_l_t = []

        for c in [x for x in c_l if x is not None]:

            a_i_l_t.append(c.i_dist())

        a_i_l.append(np.mean(a_i_l_t))

        # Calculate external distances

        e_i_l_t = []

        for ci in [x for x in c_l if x is not None]:
            for cj in [x for x in c_l if x is not None]:

                e_i_l_t.append(ci.c_dist(cj))

        e_i_l.append(np.mean(e_i_l_t))

    d_l_prev = []

    for index in xrange(1,len(d_l)):

        d_l_prev.append( d_l[index]-d_l[index-1] )

    if plot_convergence:

        plt.figure(110)
        plt.plot(d_l_prev)

    print "GUESSING" ,d_l_prev.index(max(d_l_prev))

    if d_l_prev.index(max(d_l_prev)) == 16:

        correct += 1

    all_guess += 1

    print "LENDATA" , len(data)

    for ch in c_h:

        if len(ch) == len(data) - d_l_prev.index(max(d_l_prev)):

            # Plot clusters

            cluster_list, colors, corr_id, tot_id = plot_and_count(len(data) - d_l_prev.index(max(d_l_prev)), ch, colors, cluster_list, corr_id, tot_id, plot_best)

    for i in range(18-len(d_l_prev)):

        d_l_prev.append(0)

    for i in range(18-len(a_i_l)):

        a_i_l.append(0)

    for i in range(18 - len(e_i_l)):

        e_i_l.append(0)

    return cluster_list,all_guess,correct,colors,corr_id, tot_id,d_l_prev,e_i_l,a_i_l

if __name__ == '__main__':

    # Initiate variables

    colors = []

    tot = 0
    tot_corr = 0
    corr_id = 0
    tot_id = 0

    cluster_list = []

    frame_list = []

    counter = 0

    a_d_l = []
    a_e_l = []
    a_i_l = []

    # Load data set

    if not flatfrog_data:

        pass

    if flatfrog_data:

        obj = untangle.parse('data/data1.xml')

        for frame in obj.root.frame:

            point_list = []

            try:

                for coord in frame.coord:

                    point_list.append(ff_p(int(coord["trace_id"])
                                           , float(coord["x_mm"])
                                           , float(coord["y_mm"])
                                           , float(coord["vx_mm_s"])
                                           , float(coord["vy_mm_s"])
                                           , float(coord["pressure"])
                                           , float(coord["A"])))

                frame_list.append(point_list)

            except IndexError as e:

                pass

    # Run algorithm

        # data 2 => 150:len(frame_list)-100
        # data 3 => 25:len(frame_list)-130

    for point_list in frame_list:

        if len(point_list) > 2:

            cluster_list, all_guess, correct, colors, k_corr_id, k_tot_id,d_l_t,e_l_t,i_l_t = HAC(point_list, cluster_list, colors)

            a_d_l.append(d_l_t)
            a_i_l.append(i_l_t)
            a_e_l.append(e_l_t)

            tot += all_guess
            tot_corr += correct

            corr_id += k_corr_id
            tot_id += k_tot_id

            counter += 1

            print "COUNTER", counter

    def mean(a):

        return sum(a) / len(a)

    # Plot selected graphs

    if plot_averaged_dendo:

        plt.figure()
        plt.scatter(xrange(len(map(mean, zip(*a_d_l)))), map(mean, zip(*a_d_l)))
        plt.title("Averaged Dendogram Distance")

    if plot_averaged_detrend_dendo:

        plt.figure()
        ffty = scipy.signal.detrend(list(map(mean, zip(*a_d_l))))
        plt.scatter(xrange(len(map(mean, zip(*a_d_l)))), ffty)
        plt.title("Averaged Detrand Dendogram Distance")

    if plot_averaged_external:

        plt.figure()
        plt.scatter(xrange(len(map(mean, zip(*a_e_l)))), map(mean, zip(*a_e_l)))
        plt.title("Averaged External Distance")

    if plot_averaged_detrend_external:

        plt.figure()
        ffty = scipy.signal.detrend(list(map(mean, zip(*a_e_l))))
        plt.scatter(xrange(len(map(mean, zip(*a_e_l)))), ffty)
        plt.title("Averaged External Detrend Distance")

    if plot_averaged_internal:

        plt.figure()
        plt.scatter(xrange(len(map(mean, zip(*a_i_l)))), map(mean, zip(*a_i_l)))
        plt.title("Averaged Internal Distance")

    if plot_averaged_detrend_internal:

        plt.figure()
        ffty = scipy.signal.detrend(list(map(mean, zip(*a_i_l))))
        plt.scatter(xrange(len(map(mean, zip(*a_i_l)))), ffty)
        plt.title("Averaged Internal Detrend Distance")

    if count_correct:

        print "CORRECT GUESS NUMBER CLUSTER RATIO: ", tot_corr / tot
        print "CORRECT GUESS FINGER RATIO: ", corr_id / tot_id

    plt.show()
