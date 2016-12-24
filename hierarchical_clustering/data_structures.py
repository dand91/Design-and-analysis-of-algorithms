import matplotlib.cm as cm
import numpy as np
import operator
import math
import matplotlib.pyplot as plt
import sys

class cl(object):

    def __init__(self,id):

        self.color = None
        self.id = id
        self.points = []

    def get_trace_id(self):

        trace_id_list = []

        for point in self.points:

            trace_id_list.append(point.trace_id)

        return trace_id_list

    def get_p_l(self):

        return self.points

    def add(self,point):

        if isinstance(point, list):

            self.points.extend(point)

        else:

            self.points.append(point)

    def c_dist(self,other):

        ci_l = self.get_p_l()
        cj_l = other.get_p_l()

        d_l = []

        for i in xrange(0,len(ci_l)):
            for j in xrange(i,len(cj_l)):

                ci_p = ci_l[i]
                cj_p = cj_l[j]

                d_l.append(ci_p.p_dist(cj_p))

        #return max(d_l)
        return min(d_l)
        #return np.mean(d_l)

    def i_dist(self):

        c_l = self.get_p_l()

        d_l = []

        for i in xrange(0, len(c_l)):
            for j in xrange(i, len(c_l)):

                ci_p = c_l[i]
                cj_p = c_l[j]

                d_l.append(ci_p.p_dist(cj_p))

        # return max(d_l)
        # return min(d_l)
        # return np.mean(d_l)
        return np.var(d_l)

class c_p(object):

    def __init__(self,x,y):

        self.x = x
        self.y = y

    def p_dist(self,other):

        return math.sqrt( (self.x - other.x)**2 + (self.y - other.y)**2 )

class ff_p(object):

    def __init__(self,trace_id,x_mm,y_mm,vx_mm_s,vy_mm_s,pressure,A):

        self.trace_id = trace_id
        self.x_mm = x_mm
        self.y_mm = y_mm
        self.vx_mm_s = vx_mm_s
        self.vy_mm_s = vy_mm_s
        self.pressure = pressure
        self.A = A

    def p_dist(self,other):

        l = [1, 0.1, 0, 0.01]

        return math.sqrt(l[0] * (self.x_mm - other.x_mm) ** 2 + l[0] * (self.y_mm - other.y_mm) ** 2 + \
                         l[1] * (self.vx_mm_s - other.vx_mm_s) ** 2 + l[1] * (self.vy_mm_s - other.vy_mm_s) ** 2 + \
                         l[2] * math.fabs(self.pressure - other.pressure) + \
                         l[3] * math.fabs(self.A - other.A))

class c_e(object):

    def __init__(self,dist,ci,cj):

        self.dist = dist
        self.ci = ci
        self.cj = cj

    def __lt__(self, other):

        return other.dist > self.dist

def plot_and_count(k, cluster_list, colors, pcl, corr_id, tot_id, plot_best):

    color_cluster_history = []

    for p in pcl:

        d_l = []

        p_x = []
        p_y = []

        for point in p.get_p_l():

            p_x.append(point.x_mm)
            p_y.append(point.y_mm)

        for c in cluster_list:

            c_x = []
            c_y = []

            for point in c.get_p_l():

                c_x.append(point.x_mm)
                c_y.append(point.y_mm)

            dist = math.sqrt( (sum(c_x)/len(c_x) - sum(p_x)/len(p_x) ) **2 +  (sum(c_y)/len(c_y) - sum(p_y)/len(p_y) ) **2 )

            d_l.append((dist, c))

        if len(d_l) > 0:

            if min(d_l)[0] < 10:

                color_cluster_history.append( ( min(d_l),min(d_l)[1],p ))
                color_cluster = min(d_l)[1]

                if p.color is None:

                    if len(colors) == 0:

                        colors = ["red", "green", "blue", "yellow", "purple", "orange", "magenta", "darkkhaki",
                                  "orchid", "skyblue", "red"]

                    color_cluster.color = colors.pop()

                else:

                    color_cluster.color = p.color

    for d, c, p in color_cluster_history:

        c_id = c.get_trace_id()
        c_id.append(sys.maxint)
        p_id = p.get_trace_id()
        p_id.append(sys.maxint)

        corr_id += len(list(set(c_id).intersection(p_id))) - 1

    for c in cluster_list:

        tot_id += len(c.get_trace_id())

    for cluster in cluster_list:

        X = []
        Y = []

        for point in cluster.get_p_l():

            X.append(point.x_mm)
            Y.append(point.y_mm)

        if plot_best:

            plt.figure(100)

            plt.title("Scatter plot, k = " + str(k))
            plt.xlabel("X")
            plt.ylabel("Y")
            plt.scatter(X, Y, color=cluster.color, alpha=0.75)

            plt.pause(0.05)

    return cluster_list , colors, corr_id, tot_id
