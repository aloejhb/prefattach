import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def threshold_network(r0,W,thresh,tstep):
    if W.shape[0] != W.shape[1] or W.shape[0] != r0.shape[0]:
        raise Exception('W should be a square matrix with row number same as length of r0!')
    r = np.empty((tstep, r0.shape[0]))
    r[0] = r0
    for t in range(tstep-1):
        r[t+1] = np.matmul(W, r[t]) > thresh
    return r


def get_weight_matrix(N, weight_type, param=None):
    if weight_type == 'random':
        W = [random.gauss(0, 1) for i in range(N**2)]
        W = np.reshape(W, (N, N))
    elif weight_type == 'random_bin':
        if param:
            excit_frac = param['excit_frac']
        else:
            excit_frac = 0.5
        weights = [1-excit_frac, excit_frac]
        W = random.choices([-1, 1], weights=weights, k=N**2)
        W = np.reshape(W, (N, N))
    return W


def show_movie(movie):
    fig, ax = plt.subplots()
    ims = []
    for i in range(movie.shape[2]):
        im = plt.imshow(movie[:, :, i], animated=True)
        ims.append([im])
    ani = animation.ArtistAnimation(fig, ims, interval=50, repeat_delay=1000)
    return ani


if __name__ == '__main__':
    N = 25**2
    tstep = 1000
    thresh = 1
    fplus = 0.5
    r0 = np.array([random.randrange(2) for i in range(N)])
    W = get_weight_matrix(N, 'random_bin', dict(excit_frac=fplus))
    r = threshold_network(r0, W, thresh, tstep)

    xx = int(np.sqrt(N))
    yy = int(N / xx)
    movie = np.reshape(r, (xx, yy, tstep))
    print('making movie ...')
    ani = show_movie(movie[:, :, :200])
    ani.save('dynamic_images.mp4')
    # plt.plot(r)
    # plt.show()
    # plt.plot(r[:,0:2])
    plt.show()
