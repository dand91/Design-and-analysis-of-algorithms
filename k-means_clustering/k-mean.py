from __future__ import division
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import random
import copy
import scipy.signal
import untangle
import sys
import scipy.signal

from data_structures import *

count_correct = True

flatfrog_data = True

plot_center = False

zero_guard = 0.000001
convergence_tolerance = 0.0001

plot_each_iteration = False
plot_best = False
plot_convergence = False
plot_silhouette = False
plot_BIC = False
plot_average_BIC = False
plot_average_detrend_BIC = False

plt.ion()

def calculate_BIC(point_list, cluster_list):

    #number of clusters

    k = len(cluster_list)

    #size of data set

    R = len(point_list)

    if flatfrog_data:

        M = 5

    else:

        M = 2

    logl = []

    sigma_2 = 0

    for cluster in cluster_list:

        cluster_points = cluster.get_p_l()
        cluster_center = cluster.get_center()

        for cluster_point in cluster_points:

            sigma_2 += cluster_point.dist(cluster_center) ** 2

    sigma_2_norm = (1.0 / (zero_guard + R - k) / M) * sigma_2 + zero_guard

    for cluster in cluster_list:

        cluster_center = cluster.get_center()

        R_n = len(cluster.get_p_l())

        for cluster_point in cluster.get_p_l():

            logl.append(np.log(1.0/np.sqrt(2*np.pi*sigma_2_norm)) - (1.0/(2*sigma_2_norm)) * ( cluster_point.dist(cluster_center) ** 2 ) + np.log(R_n/R))

    BIC = np.sum(logl) - 0.5 * np.log(R) * k * (M+1)

    return BIC

def k_mean_initiation(point_list,cluster_list):

    center_dist_list = []

    for point in point_list:

        cluster_dist_list = []

        for cluster in cluster_list:

            cluster_dist_list.append(math.pow(point.dist(cluster.get_center()),2))

        center_dist_list.append(min(cluster_dist_list))

    dist_sum = sum(center_dist_list)

    point_prob = []

    for dist in center_dist_list:

        point_prob.append(dist / (dist_sum + zero_guard))

    new_cluster = cl_ff()

    index = np.where(np.random.multinomial(1, point_prob))[0][0]

    new_cluster.set_center(point_list[index])

    cluster_list.append(new_cluster)

    return cluster_list

def k_mean(max_iter,point_list,cluster_list,color_iter):

    iteration = 0

    center_delta_dist_list = []

    old_cluster_center_list = []

    while iteration <= max_iter:

        iteration += 1

        print "Iteration", iteration

        # Clear cluster list to calculate new center

        for cluster in cluster_list:

            cluster.clear_list()

        # Iterate over all points to find closes cluster.

        for point in point_list:

            cluster_dist_list = []

            for cluster in cluster_list:

                cluster_dist_list.append(point.dist(cluster.get_center()))

            cluster_list[cluster_dist_list.index(min(cluster_dist_list))].add_point(point)

        # Compute new center

        new_cluster_center_list = []

        for cluster in cluster_list:

            new_center = cluster.compute_center()

            new_cluster_center_list.append(new_center)

        if len(old_cluster_center_list) == 0:

            old_cluster_center_list = new_cluster_center_list

        else:

            dist_sum = 0

            for old_center, new_center in zip(old_cluster_center_list, new_cluster_center_list):

                dist_sum += old_center.dist(new_center)

            center_delta_dist_list.append(dist_sum / k)

            old_cluster_center_list = new_cluster_center_list

        list_len = len(center_delta_dist_list)

        if list_len > 2 and math.fabs(

                        center_delta_dist_list[list_len - 2] - center_delta_dist_list[list_len - 1]) < convergence_tolerance:

            print "Found in iteration = " + str(iteration)

            iteration = max_iter + 1

    if iteration != max_iter + 1:

        print "Not found, using max iteration = " + str(iteration) + " " + str(max_iter + 1)

    if plot_convergence:

        # Plot center convergence

        plot_center_convergence(center_delta_dist_list, next(color_iter))

    return cluster_list

def plot_center_convergence(center_delta_dist_list, color):

    plt.figure(max_k + 1)

    plt.title("Convergence plot")
    plt.xlabel("Iterations")
    plt.ylabel("Center distance change")

    plt.plot(center_delta_dist_list, color=color, alpha=0.75, label='k = ' + str(k))

    plt.legend(loc='lower right')

def plot_scatter(k, cluster_list, colors, pcl, corr_id, tot_id):

    color_cluster_history = []

    v_l = []

    for p in pcl:

        d_l = []

        for c in cluster_list:

            if c not in v_l:

                d_l.append((p.get_center().dist(c.get_center()), c,p))

        if len(d_l) > 0:

            if min(d_l)[0] < 10:

                color_cluster_history.append(min(d_l))
                color_cluster = min(d_l)[1]

                v_l.append(color_cluster)

                if p.color is None:

                    if len(colors) == 0:

                        colors = ["red", "green", "blue", "yellow", "purple", "orange", "magenta", "darkkhaki",
                                  "orchid", "skyblue"]

                    color_cluster.color = colors.pop()

                else:
                    color_cluster.color = p.color

    for d, c, p in color_cluster_history:

        c_id = c.get_trace_id()
        c_id.append(sys.maxint)
        p_id = p.get_trace_id()
        p_id.append(sys.maxint)

        corr_id += len(list(set(c_id).intersection(p_id))) -1

    for c in cluster_list:

        tot_id += len(c.get_trace_id())

    for cluster in cluster_list:

        if plot_best:

            plt.figure(100)

            X = cluster.get_x_mm()
            Y = cluster.get_y_mm()

            plt.title("Scatter plot, k = " + str(k))
            plt.xlabel("X")
            plt.ylabel("Y")

            plt.scatter(X, Y, color=cluster.color, alpha=0.75)

        if plot_center:

            plt.plot(cluster.get_center().x, cluster.get_center().y, marker='x', color='black')

        plt.pause(0.05)

    return cluster_list, corr_id, tot_id


def calculate_silhouette(point_list,cluster_list):

    a = []
    b = []

    counter = 0

    for point in point_list:

        for cluster in cluster_list:

            counter += 1

            cluster_points = cluster.get_p_l()

            t = []
            t_b = []

            cluster_found = False

            for cluster_point in cluster_points:

                dist = point.dist(cluster_point)

                t.append(dist)

                if dist == 0:

                    cluster_found = True

            if cluster_found:

                a.extend(t)

            else:

                try:
                    t_b.append(sum(t) / len(t))

                except:

                    t_b.append(0)

        if len(t_b) == 0:

            b.append(0)

        else:

            b.append(min(t_b))

    a_mean = sum(a) / len(a)

    if len(b) == 0:

        b_mean = 0

    else:

        b_mean = sum(b) / len(b)

    return (b_mean - a_mean + zero_guard) / (max(b_mean, a_mean) + zero_guard)


orig_x = []
orig_y = []

frame_list = []

if not flatfrog_data:

    point_list = []

    with open("data/R15.txt") as f:

        for line in f:

            x_y = line.split()

            orig_x.append(float(x_y[0]))
            orig_y.append(float(x_y[1]))

            point_list.append(cl_point(float(x_y[0]), float(x_y[1])))

    frame_list.append(point_list)

if flatfrog_data:

    obj = untangle.parse('data/data1.xml')


    for frame in obj.root.frame:

        print "iter"
        temp_list = []
        try:
            for coord in frame.coord:
                temp_list.append(

                    ff_point(int(coord["trace_id"]), float(coord["x_mm"]), float(coord["y_mm"]), float(coord["vx_mm_s"]),
                             float(coord["vy_mm_s"]), float(coord["pressure"]), float(coord["A"]))

                )

            frame_list.append(temp_list)

        except:
            pass


prev_cluster_list = []
bic_average = []

corr_id = 0
tot_id = 0

colors = []
frame_counter = 0

tot_clust_nbr = 0
corr_clust_nbr = 0

# data 2 => 150:len(frame_list)-100
# data 3 => 25:len(frame_list)-130

for point_list in frame_list:

    frame_counter += 1

    print "Frame number: ", frame_counter

    # Constants and initiations

    min_k = 1
    max_k = 7

    color_iter = iter(cm.rainbow(np.linspace(0, 1, max_k)))

    silhouette_list = []
    bic_list = []

    # k-mean++

    max_iter = 50

    bic_pos = 1

    cluster_matrix = []

    for k in xrange(min_k,max_k+1):

        # Add first cluster uniformly

        cluster_list = []

        new_cluster = cl_ff()

        new_cluster.set_center(point_list[random.randint(0, len(point_list) - 1)])

        cluster_list.append(new_cluster)

        print "Running k = " + str(k)

        for kk in xrange(1,k):

            cluster_list = k_mean_initiation(point_list, cluster_list)  # Initiation

        cluster_list = k_mean(max_iter,point_list,cluster_list,color_iter)

        cluster_matrix.insert(k, copy.deepcopy(cluster_list))

        if plot_each_iteration:

            plot_scatter(k, cluster_list, colors)

        silhouette_list.append(calculate_silhouette(point_list,cluster_list))

        bic_list.append(calculate_BIC(point_list, cluster_list))

    bic_average.append(copy.deepcopy(bic_list))

    try:

        bic_list = [entry for entry in bic_list if not math.isnan(entry)]
        bic_list_detrend = list(scipy.signal.detrend(bic_list))

    except:

        print bic_list

    best_k = bic_list_detrend.index(max(bic_list_detrend[1:])) + 1

    if count_correct:

        if best_k == 4:

            corr_clust_nbr += 1

    tot_clust_nbr += 1

    print "BEST K" , best_k

    found = False
    for cl in cluster_matrix:

        if len(cl) == best_k:

            found = True

            prev_cluster_list, corr_id, tot_id = plot_scatter(best_k, cl , colors,prev_cluster_list, corr_id, tot_id)

    assert found == True

    if plot_BIC:

        plt.figure(101)

        plt.clf()
        plt.title("BIC")
        plt.xlabel("Iteration")
        plt.xlabel("Value")
        plt.scatter(xrange(min_k, k+1),bic_list)

    if plot_silhouette:

        plt.figure(102)

        plt.title("Silhouette")
        plt.xlabel("Iteration")
        plt.xlabel("Value")
        plt.plot(xrange(min_k, k + 1), silhouette_list)

def mean(a):

    return sum(a) / len(a)

if plot_average_BIC:

    plt.figure(101)

    plt.clf()
    plt.title("BIC")
    plt.xlabel("Iteration")
    plt.xlabel("Value")

    plt.scatter(xrange(len(map(mean, zip(*bic_average)))), map(mean, zip(*bic_average)))

if plot_average_detrend_BIC:

    ffty=scipy.signal.detrend(list(map(mean, zip(*bic_average))))
    plt.scatter(xrange(len(map(mean, zip(*bic_average)))), ffty)
    plt.show()

if count_correct:

    print "CORRECT GUESS NUMBER CLUSTER RATIO: ", corr_clust_nbr / tot_clust_nbr
    print "CORRECT GUESS FINGERS RATIO: ", corr_id/tot_id

plt.show()
