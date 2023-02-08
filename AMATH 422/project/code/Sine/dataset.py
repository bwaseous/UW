"""generating input and target"""

import numpy as np
import torch.utils.data as data


class SineWave(data.Dataset):
    def __init__(self, time_length=50, freq_range=10):
        self.time_length = time_length
        self.freq_range = freq_range

    def __len__(self):
        return 200

    def __getitem__(self, item):
        freq = np.random.randint(4, self.freq_range + 4)
        const_signal = np.repeat(freq / self.freq_range + 0.25, self.time_length)
        const_signal = np.expand_dims(const_signal, axis=1)
        t = np.arange(0, self.time_length*0.005, 0.005)
        target = np.sin(freq * t)
        target = np.expand_dims(target, axis=1)
        return const_signal, target


'''
a = 0.7
b = 0.8
r = 0.08

v_init = 0
w_init = 0


dt = 0.3
nsteps = 500

for i in range(50):
    v = np.zeros([nsteps])
    w = np.zeros([nsteps])
    t = np.zeros([nsteps])

    v[0] = v_init
    w[0] = w_init
    t[0] = 0.0

    I_range = 10
    I = np.random.uniform(0.4, 0.4 + I_range/10)

    inputs = np.zeros([nsteps])
    outputs = np.zeros([nsteps])
    outputs[0]= v[0]
    inputs[0] = I

    for step in np.arange(nsteps-1):
        v[step+1] = v[step] + (v[step]- ((v[step])**3) /3 -w[step] + I)*dt
        w[step+1] = w[step] + ((v[step] +a - b*w[step]) * r)*dt
        t[step+1] = t[step] + dt
        outputs[step+1] = v[step+1]
        inputs[step+1] = I
    outputs[-1] = outputs[-2]
    inputs[-1] = inputs[-2]

'''