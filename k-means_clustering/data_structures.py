import math

zero_guard = 0.000001

class cl_point(object):

    def __init__(self,x,y):

        self.x = x
        self.y = y

    def dist(self, other):

        return math.sqrt( math.pow( self.x - other.x , 2 ) + math.pow( self.y - other.y , 2 ) )

    def __cmp__(self, other):

        return self.x == other.x and self.y == other.y

class ff_point(object):

    def __init__(self, trace_id, x_mm, y_mm, vx_mm_s, vy_mm_s, pressure, a):

        self.trace_id = trace_id
        self.x_mm = x_mm
        self.y_mm = y_mm
        self.vx_mm_s = vx_mm_s
        self.vy_mm_s = vy_mm_s
        self.pressure = pressure
        self.a = a

    def dist(self, other):

        l = [1, 0.1, 0, 0.01]

        return math.sqrt( l[0]*(self.x_mm - other.x_mm) ** 2 + l[0]*(self.y_mm - other.y_mm) ** 2 +\
                l[1]*(self.vx_mm_s - other.vx_mm_s) ** 2 + l[1]*(self.vy_mm_s - other.vy_mm_s) ** 2 +\
                l[2]*math.fabs(self.pressure - other.pressure) + \
                l[3]*math.fabs(self.a - other.a))

class cl_ff(object):

    def __init__(self):

        self.center = ff_point(0, 0, 0, 0, 0, 0, 0)
        self.ff_point_list = []
        self.color = None

    def add_point(self, ff_point):

        self.ff_point_list.append(ff_point)

    def get_point(self, index):

        return self.ff_point_list[index]

    def get_point_list(self):

        return self.ff_point_list

    def get_trace_id(self):

        trace_list = []

        for point in self.ff_point_list:
            trace_list.append(point.trace_id)

        return trace_list

    def get_x_mm(self):

        x_list = []

        for point in self.ff_point_list:
            x_list.append(point.x_mm)

        return x_list

    def get_y_mm(self):

        y_list = []

        for point in self.ff_point_list:
            y_list.append(point.y_mm)

        return y_list

    def get_vx_mm_s(self):

        x_list = []

        for point in self.ff_point_list:
            x_list.append(point.vx_mm_s)

        return x_list

    def get_vy_mm_s(self):

        y_list = []

        for point in self.ff_point_list:
            y_list.append(point.vy_mm_s)

        return y_list

    def get_a(self):

        y_list = []

        for point in self.ff_point_list:
            y_list.append(point.a)

        return y_list

    def get_pressure(self):

        y_list = []

        for point in self.ff_point_list:
            y_list.append(point.pressure)

        return y_list

    def get_center(self):

        return self.center

    def set_center(self,center):

        self.center = center

    def len(self):

        return len(self.ff_point_list)

    def clear_list(self):

        self.ff_point_list = []

    def compute_center(self):

        X = self.get_x_mm()
        Y = self.get_y_mm()
        VX = self.get_vx_mm_s()
        VY = self.get_vy_mm_s()
        A = self.get_a()
        P = self.get_pressure()

        new_center = ff_point(0,sum(X)/ (len(X) + zero_guard)
                              ,sum(Y)/(len(Y) + zero_guard)
                              ,sum(VX)/(len(VX) + zero_guard)
                              ,sum(VY)/(len(VY) + zero_guard)
                              ,sum(A)/(len(A) + zero_guard)
                              ,sum(P)/(len(P) + zero_guard))

        self.center = new_center

        return new_center

    class cl(object):

        def __init__(self):

            self.center = cl_point(0, 0)
            self.c_point_list = []

        def add_point(self, c_point):

            self.c_point_list.append(c_point)

        def get_point(self, index):

            return self.c_point_list[index]

        def get_point_list(self):

            return self.c_point_list

        def get_x(self):

            x_list = []

            for point in self.c_point_list:
                x_list.append(point.x)

            return x_list

        def get_y(self):

            y_list = []

            for point in self.c_point_list:
                y_list.append(point.y)

            return y_list

        def get_center(self):

            return self.center

        def set_center(self, center):

            self.center = center

        def len(self):

            return len(self.c_point_list)

        def clear_list(self):

            self.c_point_list = []

        def compute_center(self):

            X = self.get_x()
            Y = self.get_y()

            new_center = cl_point(sum(X) / len(X), sum(Y) / len(Y))

            self.center = new_center

            return new_center