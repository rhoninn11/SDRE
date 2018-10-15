import matplotlib.pyplot as plt

class symulacja:
    M = 0
    C = 0
    G = 0
    M_obliczeniowa = 0
    C_obliczeniowa = 0
    G_obliczeniowa = 0
    q1 = 0
    q2 = 0
    q1_values = 0
    q2_values = 0
    q1_values_last = 0
    q2_values_last = 0
    Tn = 0.01
    plot_q11 = []
    plot_q12 = []
    plot_q13 = []
    plot_q21 = []
    plot_q22 = []
    plot_q23 = []
    plot_t = []

    def __init__(self, M, C, G, q1, q2):
        self.M = M
        self.C = C
        self.G = G
        self.q1 = q1
        self.q2 = q2

    def simulate(self, q1_initial, q2_initial):
        Tsim=0
        self.q1_values = self.q1_values_last = q1_initial
        self.q2_values = self.q2_values_last = q2_initial
        while Tsim < 4:
            self.insert_cords_values(self.q1_values, self.q2_values)
            self.q1_values = self.q2_values * self.Tn + self.q1_values_last
            self.q2_values = self.q2_values_last + self.Tn*(-self.M_obliczeniowa*self.C_obliczeniowa*self.q2_values - self.M_obliczeniowa*self.G_obliczeniowa)
            self.q1_values_last = self.q1_values
            self.q2_values_last = self.q2_values
            self.plot_q11.append(self.q1_values[0])
            self.plot_q12.append(self.q1_values[1])
            self.plot_q13.append(self.q1_values[2])
            self.plot_q21.append(self.q2_values[0])
            self.plot_q22.append(self.q2_values[1])
            self.plot_q23.append(self.q2_values[2])
            self.plot_t.append(Tsim)
            Tsim += self.Tn
            print(str(Tsim))

        plt.plot(self.plot_t, self.plot_q11, label='q11')
        plt.plot(self.plot_t, self.plot_q12, label='q12')
        plt.plot(self.plot_t, self.plot_q13, label='q13')
        plt.plot(self.plot_t, self.plot_q21, label='q21')
        plt.plot(self.plot_t, self.plot_q22, label='q22')
        plt.plot(self.plot_t, self.plot_q23, label='q23')
        plt.legend()
        plt.show()



    def insert_cords_values(self, x1, x2):
        for i in range(len(self.q1)):
            if i is 0:
                self.M_obliczeniowa = self.M.subs(self.q1[i], x1[i])
                self.C_obliczeniowa = self.C.subs(self.q2[i], x2[i])
                self.C_obliczeniowa = self.C_obliczeniowa.subs(self.q1[i], x1[i])
                self.G_obliczeniowa = self.G.subs(self.q1[i], x1[i])
            else:
                self.M_obliczeniowa = self.M_obliczeniowa.subs(self.q1[i], x1[i])
                self.C_obliczeniowa = self.C_obliczeniowa.subs(self.q2[i], x2[i])
                self.C_obliczeniowa = self.C_obliczeniowa.subs(self.q1[i], x1[i])
                self.G_obliczeniowa = self.G_obliczeniowa.subs(self.q1[i], x1[i])



