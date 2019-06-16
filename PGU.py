import math as m
import cmath as cm
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp



K = 6.
rMIN = 0.0
rMAX = 1.8
a = 0.6
b = 1.2
EXT = 1

L = 65
I = 120
d = 1
C = 0.8
c = 1.0
N = 2.5
h = (rMAX - rMIN) / I
t = C * h / c
nt = int(N / t)
r_ext = [rMIN + (i - 0.5)*h for i in range(I*2 + 2)]
r = [rMIN + (i - 0.5)*h for i in range(I + 2)]


#######Чтение коэффициантов из файла###########
f = open('/home/bobbi/PycharmProjects/Sofr/coef' + str(L) + '.dat', 'r') #Путь к файлам с коэффициентами
lines = [line.split() for line in f.readlines()]
alpha_r = []
alpha_im = []
beta_r = []
beta_im = []
for i in range(L):
    alpha_r.append(float(lines[i][0]))
    alpha_im.append(float(lines[i][1]))
    beta_r.append(float(lines[i][2]))
    beta_im.append(float(lines[i][3]))
alpha = [K**2 *c*(alpha_r[i] + 1j*alpha_im[i]) for i in range(len(alpha_r))]
beta = [K*c*(beta_r[i] + 1j*beta_im[i]) for i in range(len(beta_r))]


#####Свертка через полный интеграл#####
def integral(A, curr, next, nt):
    A = 0
    for n in range(nt+1):
        for l in range(L):
            A += alpha[l]*cm.exp(beta[l]*((nt - n)*t + t/2)) *0.5*(curr[n]+next[n])*t
    return A


####Нач функция####
def u0(I, Ext):
    h = (rMAX - rMIN) / I
    t = C * h / c
    nt = int(N / t)

    def v0(r):
        if (r < b) & (r > a):
            return m.exp((-4*(2*r - (a + b))**2) / ((b - a)**2 - (2*r - (a + b))**2))
        else:
            return 0

    r = [rMIN + (i - 0.5)*h for i in range(I*((Ext == 1) + 1) + 2)]
    u0 = [v0(r[i]) for i in range(I + 2)]
    return u0

####Решение#####
def u(I, K, Ext):
    h = (rMAX - rMIN) / I
    I = I * ((Ext == 1) + 1)
    t = C * h / c
    nt = int(N / t)

    def v0(r):
        if (r < b) & (r > a):
            return m.exp((-4*(2*r - (a + b))**2) / ((b - a)**2 - (2*r - (a + b))**2))
        else:
            return 0

    r = [rMIN + (i - 0.5)*h for i in range(I + 2)]

    a1 = [t**2 * c**2 * r[i]**(1-d) / h for i in range(I + 2)]
    a2 = [(r[i] + h/2)**(d-1) / h for i in range(I + 2)]
    a3 = [(r[i] - h/2)**(d-1) / h for i in range(I + 2)]

    u0 = [v0(r[i]) for i in range(I + 2)]
    u1 = list(u0)

    u_next = [0 for i in range(I + 2)]

    u11 = [0 for i in range(I + 2)]
    u12 = [0 for i in range(I + 2)]
    u13 = [0 for i in range(I + 2)]
    for i in range(1, I + 1):
        u12[i] = a1[i] * (a2[i] * (u0[i + 1] - u0[i]) - a3[i] * (u0[i] - u0[i - 1]))

    # u12[0] = u12[1]
    # u12[I + 1] = u12[I]

    for i in range(1, I + 2):
        u1[i] = u0[i] + 0.5 * u12[i] - (t*K)**2 / 2 * u0[i]

    u_prev = list(u0)
    u_curr = list(u1)

    sol = [[0 for i in range(I)] for n in range(nt)]
    sol[0] = u0[1:I+1]
    Conv = [0 for l in range(L)]
    # A = 0
    for n in range(1, nt):
        sol[n] = u_curr[1:I + 1]
        for i in range(1, I + 1):
            # u_next[i] = 2*u_curr[i] - u_prev[i] + a1[i] * (a2[i] * (u_curr[i + 1] - u_curr[i]) - a3[i] * (u_curr[i] - u_curr[i - 1])) - t**2*c**2*K**2*u_curr[i]
            u_next[i] = 2*u_curr[i] - u_prev[i] + (t*c/h)**2 * (u_curr[i + 1] - 2*u_curr[i] + u_curr[i - 1]) - (t*c*K)**2 * u_curr[i]

        #####Свертка итерациями####
        for l in range(L):
            Conv[l] = cm.exp(beta[l]*t)*Conv[l] + alpha[l]*t/12.*(cm.exp(beta[l]*t)*(u_curr[I] + u_prev[I+1]) + \
                                                                   4.*cm.exp(beta[l]*t/2)*(u_curr[I] + u_curr[I+1]) + \
                                                                   + (u_next[I] + u_curr[I+1]))

        #A = integral(A, curr, next, n)

        u_prev = list(u_curr)
        u_curr = list(u_next)

        u_curr[0] = -u_curr[1]

        u_curr[I+1] = h / (h + t*c) * (2*t*c*(np.real(sum(Conv))) + t*c/h*(u_curr[I] + u_prev[I] - u_prev[I+1]) + u_prev[I+1] + u_prev[I] - u_curr[I])

    return sol

def NT(I1, Ext):
    h1 = (rMAX - rMIN) / I1
    t1 = C * h1 / c
    nt1 = int(N / t1)
    return nt1

nt1 = NT(I, EXT)
nt2 = NT(I*3, EXT)
nt3 = NT(I*9, EXT)
nt = nt1 - 2

u1_ext = u(I, K, EXT)
# u2_ext = u(I*3, K, EXT)
# u3_ext = u(I*9, K, EXT)
u1 = u(I, K, EXT-1)
u2 = u(I*3, K, EXT-1)
u3 = u(I*9, K, EXT-1)

# for n in range(nt1):
#     plt.plot(r, u1[n])
#     # plt.plot(r[1:51], u2[50*n : 50*(n+1)])
#     # plt.plot(r[1:51], u3[50*n : 50*(n+1)])
#     plt.savefig('/home/bobbi/sols/' + str(n) + '.png')
#     plt.close()

d1 = [[0 for i in range(I)] for n in range(nt)]
d2 = [[0 for i in range(I)] for n in range(nt)]
for n in range(nt):
    for i in range(I):
        j = 3*i + 1
        k = 9*i + 4
        d1[n][i] = abs(u2[3*n][j] - u1[n][i])
        d2[n][i] = abs(u3[9*n][k] - u2[3*n][j])

u1m = [0 for i in range(nt)]
u2m = [0 for i in range(nt)]
for n in range(nt):
    u1m[n] = max(d1[n])
    u2m[n] = max(d2[n])

alpha1 = []
for n in range(1, nt):
    alpha1.append(u1m[n] / u2m[n])
    # alpha2.append(u2m[n] / u3m[n])

d3 = [[0 for i in range(I)] for n in range(nt)]
for n in range(nt):
    for i in range(I):
        d3[n][i] = abs(u1[n][i] - u1_ext[n][i])

time = [n*t for n in range(nt)]
for n in range(nt):
    plt.ylim(-3, 3)
    plt.plot(r_ext[1:I*2+1], u1_ext[n])
    plt.plot(r[1:I+1], u1[n])
    plt.savefig('/home/bobbi/PycharmProjects/Sofr/diff/{}.png'.format(n))
    plt.close()

plt.xlabel("$t, s.$")
plt.ylabel("$1 / alpha^2$")
plt.grid()
plt.ylim(0, 15)
plt.plot(time[1:], alpha1)
# plt.show()
plt.savefig('/home/bobbi/PycharmProjects/Sofr/PGU_solution' + str(int(K)) + '/differ.png')
for n in range(nt):
    # plt.xlabel("$r$")
    # plt.ylabel("$u_{ext}(r)-u_{ord}(r)$")
    # plt.subplot(211)
    plt.subplot(212)
    plt.ylabel("$||u_{ext}-u_{ord}||$")
    plt.grid()
    plt.ylim(0, 0.1)
    plt.plot(r[1:I+1], d3[n], label = "difference")
    plt.legend()
    plt.subplot(211)
    plt.title("$\kappa = ${0} $I = ${1} $L = ${2} $T = ${3:.2f} $s.$".format(K, I, L, n*t))
    plt.ylabel("$u(r)$")
    plt.grid()
    plt.ylim(-1.5, 1.5)
    plt.plot(r[1:I+1], u1[n], label = "ordinary solution")
    plt.plot(r_ext[1:I*2+1], u1_ext[n], label = "external solution")
    plt.legend()
    plt.savefig('/home/bobbi/PycharmProjects/Sofr/PGU_solution' + str(int(K)) + '/' + str(nt - n) + '.png')
    plt.close()
#
