import numpy as np
import matplotlib.pyplot as plt

def SEIR_Low_High(initial_conditions, R0_low2low, R0_high2low, R0_high2high, R0_low2high, tau_i, tau_e, num_steps):
    S_low, S_high, E_low, E_high, I_low, I_high, R_low, R_high = np.zeros(num_steps), np.zeros(num_steps), np.zeros(num_steps), np.zeros(num_steps), np.zeros(num_steps), np.zeros(num_steps), np.zeros(num_steps), np.zeros(num_steps)
    S_low[0], S_high[0], E_low[0], E_high[0], I_low[0], I_high[0], R_low[0], R_high[0] = initial_conditions
    N_low = S_low[0]+E_low[0]+I_low[0]+R_low[0]
    N_high = S_high[0]+E_high[0]+I_high[0]+R_high[0]
    for t in range(num_steps-1):
        S_low[t+1] = S_low[t] - (R0_low2low/tau_i) * (S_low[t]/N_low) * I_low[t] - (R0_high2low/tau_i) * (S_low[t]/N_low) * I_high[t]
        S_high[t+1] = S_high[t] - (R0_high2high/tau_i) * (S_high[t]/N_high) * I_high[t] - (R0_low2high/tau_i) * (S_high[t]/N_high) * I_low[t]
        E_low[t+1] = E_low[t] + S_low[t] - S_low[t+1] - (1/tau_e)*E_low[t]
        E_high[t+1] = E_high[t] + S_high[t] - S_high[t+1] - (1/tau_e)*E_high[t]
        I_low[t+1] = I_low[t] - (1/tau_i)*I_low[t] + (1/tau_e)*E_low[t]
        I_high[t+1] = I_high[t] - (1/tau_i)*I_high[t] + (1/tau_e)*E_high[t]
        R_low[t+1] = R_low[t] + (1/tau_i)*I_low[t]
        R_high[t+1] = R_high[t] + (1/tau_i)*I_high[t]
    return np.vstack([S_low, S_high, E_low, E_high, I_low, I_high, R_low, R_high]).T


def show_SEIR_Low_high(initial_conditions, R0_low2low, R0_high2low, R0_high2high, R0_low2high, tau_i, tau_e, num_steps, sick_days,
        probability_severe_given_low, probability_severe_given_high):

    data = SEIR_Low_High(initial_conditions, R0_low2low, R0_high2low, R0_high2high, R0_low2high, tau_i, tau_e, num_steps)
    m_low = np.sum(data[0][::2])
    m_high = np.sum(data[0][1::2])

    subplot_base = 111
    ax1=plt.subplot(subplot_base)
    ax1.plot(np.arange(num_steps), 100*data[:,4]/m_low, 'b-', np.arange(num_steps), 100*data[:,5]/m_high, 'r-')
    ax1.legend(['% sick_low', '% sick_high'])
    ax1.grid(True)
    ax1.set_xlabel('days')
    plt.savefig('mixedSEIR_percent_sick.png')
    plt.clf()

    ax=plt.subplot(subplot_base)
    ax.plot(np.arange(num_steps), data[:,4], 'b-', np.arange(num_steps), data[:,5], 'r-')
    ax.legend(['#sick_low', '#sick_high'])
    ax.grid(True)
    ax.set_xlabel('days')
    plt.savefig('mixedSEIR_num_sick.png')
    plt.clf()

    ax=plt.subplot(subplot_base)
    ax.plot(np.arange(num_steps), probability_severe_given_low*data[:,4]*sick_days/tau_i, 'b-', np.arange(num_steps), probability_severe_given_high*data[:,5]*sick_days/tau_i, 'r-')
    ax.legend(['#severe_low', '#severe_high'])
    ax.grid(True)
    ax.set_xlabel('days')
    plt.savefig('mixedSEIR_num_high.png')
    plt.clf()


if __name__ == '__main__':
    initial_conditions = [6996050.0, 1997150.0, 2500.0, 1500.0, 1250.0, 750.0, 200.0, 600.0]
    show_SEIR_Low_high(initial_conditions, 1.4, 0.02, 0.7, 0.02, 2.9, 5.0, 250, 11, 0.0025, 0.04)
