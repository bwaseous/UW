"""training models"""

import argparse
import os
import sys

import numpy as np
import torch
import torch.optim as optim

sys.path.append('../')

from torch.autograd import Variable

from dataset import SineWave
from model import RecurrentNeuralNetwork


def main(activation):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(device)
    print(torch.cuda.device_count())
    print(torch.cuda.get_device_name(torch.cuda.current_device()))

    os.makedirs('trained_model', exist_ok=True)
    save_path = f'trained_model/{activation}'
    os.makedirs(save_path, exist_ok=True)

    hidden_size = 256 #512#####################
    model = RecurrentNeuralNetwork(n_in=1, n_out=1, n_hid=hidden_size, device=device,
                                   activation=activation, sigma=0, use_bias=True).to(device)
                                    #n_hid=200
    train_dataset = SineWave(freq_range=20, time_length=100) #freq_range=20, time_length=200
    seed = np.random.seed()
    batch_size = 200 #200
    train_dataloader = torch.utils.data.DataLoader(train_dataset, batch_size=batch_size,
                                                   num_workers=2, shuffle=True,
                                                   worker_init_fn=lambda x: seed)

    print(model)

    optimizer = optim.Adam(filter(lambda p: p.requires_grad, model.parameters()),
                           lr=0.0001, weight_decay=0.0001) #lr = 0.0001##############
                #optim.LBFGS(filter(lambda p: p.requires_grad, model.parameters()),
                #            lr=0.001)

    for epoch in range(2001): #2001
        model.train()
        for i, data in enumerate(train_dataloader):
            inputs, target, = data
            #inputs, target = torch.nn.functional.normalize(inputs), torch.nn.functional.normalize(target)
            inputs, target, = inputs.float(), target.float()
            inputs, target = Variable(inputs).to(device), Variable(target).to(device)

            hidden = torch.zeros(batch_size, hidden_size)
            hidden = hidden.to(device)

            optimizer.zero_grad()
            hidden = hidden.detach()
            hidden_list, output, hidden = model(inputs, hidden)

            loss = torch.nn.MSELoss()(output, target)
            loss.backward()
            #torch.nn.utils.clip_grad_norm_(model.parameters(), 5)
            optimizer.step()

        if epoch > 0 and epoch % 1 == 0:
            print(f'Train Epoch: {epoch}, Loss: {loss.item():.6f}')
            #print('output', output[0, :, 0].cpu().detach().numpy())
            #print('target', target[0, :, 0].cpu().detach().numpy())
            torch.save(model.state_dict(), os.path.join(save_path, f'epoch_{epoch}.pth'))
            if loss.item() < 0.001:
                exit()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='PyTorch RNN training')
    parser.add_argument('--activation', type=str, default='relu')
    args = parser.parse_args()
    # print(args)
    main(args.activation)
