import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def plot_availability_curve(mttf, delta_t, mu_c, mu_p, c, t_final):
    lambd = 1 / mttf
    n_final = int(t_final / delta_t)
    t = np.arange(0.0, t_final, delta_t)
    a_matrix = np.array([[1-3*lambd*delta_t, mu_c*delta_t, mu_p*delta_t, mu_c*delta_t, mu_c*delta_t],
                    [2*lambd*delta_t, 1-lambd*delta_t-mu_c*delta_t, 0, 0, 0],
                    [(1-c)*lambd*delta_t, 0, 1-2*lambd*delta_t - mu_p*delta_t, 0, 0],
                    [lambd*c*delta_t, 0, 0, 1-2*lambd*delta_t-mu_c*delta_t, 0],
                    [0, lambd*delta_t, 2*lambd*delta_t, 2*lambd*delta_t, 1-mu_c*delta_t]])

    p_zero = np.array([1, 0, 0, 0, 0])
    r = np.zeros(n_final)
    p = np.zeros(5)
    a_n = np.zeros((5, 5))

    for n in range(0, n_final):
        a_n = np.linalg.matrix_power(a_matrix, n)
        p = np.matmul(a_n, p_zero)
        r[n] = 1 - p[4]

    fig, ax = plt.subplots()
    ax.plot(t, r, label='Availability')

    xlim = ax.get_xlim()
    label = 'Asymptote = ' + str(round(r[n_final-1] - 0.00005, 4))
    ax.hlines(round(r[n_final-1] - 0.00005, 4), xlim[0], xlim[1], linestyles='dashed', color='black', label=label)

    ax.set(xlabel='time (h)', ylabel='Availability', title='Availability over time')

    ax.grid()

    fig.savefig('plots/availabity.png')
    plt.legend()
    plt.show()

def plot_reliability_curve(mttf, delta_t, mu_c, mu_p, c, t_final):
    lambd = 1 / mttf
    n_final = int(t_final / delta_t)
    t = np.arange(0.0, t_final, delta_t)
    a_matrix = np.array([[1-3*lambd*delta_t, mu_c*delta_t, mu_p*delta_t, mu_c*delta_t, 0],
                    [2*lambd*delta_t, 1-lambd*delta_t-mu_c*delta_t, 0, 0, 0],
                    [(1-c)*lambd*delta_t, 0, 1-2*lambd*delta_t - mu_p*delta_t, 0, 0],
                    [lambd*c*delta_t, 0, 0, 1-2*lambd*delta_t-mu_c*delta_t, 0],
                    [0, lambd*delta_t, 2*lambd*delta_t, 2*lambd*delta_t, 1]])

    p_zero = np.array([1, 0, 0, 0, 0])
    r = np.zeros(n_final)
    p = np.zeros(5)
    a_n = np.zeros((5, 5))

    for n in range(0, n_final):
        a_n = np.linalg.matrix_power(a_matrix, n)
        p = np.matmul(a_n, p_zero)
        r[n] = 1 - p[4]

    mttf_sys = np.trapz(r, t)
    print("MTTF system = " + str(mttf_sys))

    fig, ax = plt.subplots()
    ax.plot(t, r)

    ax.set(xlabel='time (h)', ylabel='Reliability', title='Reliability over time')

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.text(xlim[0] + 0.1, ylim[0] + 0.1, 'MTTFsystem = ' + str(round(mttf_sys, 2)))
    ax.grid()

    fig.savefig('plots/reliability.png')
    plt.show()
    
def main():
    mttf = 50
    delta_t = 0.01
    mu_c = 0.4
    mu_p = 0.3
    c = 0.6
    t_final = 50

    plot_availability_curve(mttf, delta_t, mu_c, mu_p, c, t_final)

    t_final = 1700

    plot_reliability_curve(mttf, delta_t, mu_c, mu_p, c, t_final)

if __name__ == "__main__":
    main()

