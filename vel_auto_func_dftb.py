#!/usr/bin/env python
# coding=utf-8

# ====== Imports ======
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import spectrogram
#import timing


# ====== Init params (CHECK THIS BEFORE RUNNING!) ======
initialIter = 0  # Set this to the expected value for first iteration number
inicial_dt = 40e-9
deltaT = 0.5*(10^(-15))
samplingRate = 1 / deltaT

# ====== Processing methods ======
def parse_xyz(filename):

    with open(filename) as f:

        i = 0
        nsteps = 1
        coords = []
        chgs = []
        vels = []
        line = next(f)
        while True:

            try:

                if i == 0:
                    nat = int(line.strip().split()[0])
                    line = next(f)
                    line = next(f)
                    i += 1

                j = 0
                while j < nat:
                    data = line.strip().split()
                    data[1:] = list(map(float, data[1:]))
                    coors = data[1:4]
                    chg = data[4]
                    vel = data[5:]
                    coords.append(coors)
                    chgs.append(chg)
                    vels.append(vel)
                    line = next(f)
                    j += 1

                nsteps += 1
                i = 0

            except StopIteration:
                break

    coords = np.array(coords).reshape(nsteps, nat, -1)
    chgs = np.array(chgs).reshape(nsteps, nat, -1)
    vels = np.array(vels).reshape(nsteps, nat, -1)

    return coords, chgs, vels


def calculateFFT(dT, signalInfo): # TODO: Check what we return here
    """
    Calculates the Fast Fourier Transform for a given sample
    """
    t = np.arange(0, 1000e-6, dT)
    fscale = t / max(t)  # normalize t [0..1]
    y = (np.exp(-2 * np.pi * 2e6 * t) * np.cos( 2 * np.pi * 2e6 * t * fscale ))  # Fancy calculation
    
    #TODO: Do we use autocorrelation as signalData (before kaiser)?
    #signalData = signalInfo
    
    #y *= np.hanning(len(y))
    y *= np.kaiser( len(y), 6 )  # modifies FancyCalc array using this Bessel function
    signalData = np.concatenate( (y, ([0] * 10 * len(y) )) )  # kaiser works till signalLength-1. signalLength term is added
    
    # FFT calculation
    Fs = 1 / dT  # sampling rate, Fs = 500MHz = 1/2ns
    signalLength = len(signalData)  # length of the signal
    k = np.arange(signalLength)  # aritmethic progression. Step 1. Zero to length of signal data.
    T = signalLength / Fs  # Period of signal
    frequencyRange = k / T  # two sides frequency range. Nyquist!
    frequencyRange = frequencyRange[range(signalLength / 2)]  # one side frequency range
    Y = np.fft.fft(signalData) / signalLength  # FFT computing
    Y = Y[range(signalLength / 2)] / max(Y[range(signalLength / 2)])  # FFT normalization
    return y, frequencyRange, Y


def plotData(y, t, dT, freqRange, fftResult):  # TODO: Check validity of what we return here
    # plotting the data
    xlabelText = "Time (micro seconds)"
    ylabelText = "Amplitude"
    plt.subplot(3, 1, 1)
    plt.plot(t * 1e3, y, 'r')
    plt.xlabel(xlabelText)
    plt.ylabel(ylabelText)
    plt.grid()
    
    # plotting the spectrum
    plt.subplot(3, 1, 2)
    plt.plot(freqRange[0:600], abs(fftResult[0:600]), 'k')
    plt.xlabel('Freq (Hz)')
    plt.ylabel('|Y(freq)|')
    plt.grid()
    
    # plotting the specgram
    plt.subplot(3, 1, 3)
    Pxx, freqs, bins, im = plt.specgram(y, NFFT=512, Fs=1/dT, noverlap=10)
    plt.show()
    return Pxx, freqs, bins, im


def acf(series):
    '''Returns the autocorrelation function of a time series.'''
    N = len(series)
    avg = np.mean(series)
    c0 = np.sum((series - avg)**2) / N

    def r(j):
        return np.sum((series[:N - j] - avg) * (series[j:] - avg)) / (N - j)

    t = np.arange(N)
    acf_t = list(map(r, t))

    return acf_t #/ c0


if __name__ == '__main__':  # A main method as entrypoint for execution 
    """
    Main method.
    Get the data source file name from the call to the program. If it is not 
    available, a message with the correct use of the script is shown and the 
    program exits.
    """
    try: 
        dataSourceFileName = sys.argv[1];
    except:
        print ("========================================")
        print ("Please check your syntax. \n Usage: $ python", sys.argv[0], "resultsdftb.xyz")
        print (" resultsdftb.xyz is your data file with the following structure:")
        print ("+ ---- + --- +")
        print ("|        N   |")
        print ("+ ---- + --- +")
        print ("|MD iter: M  |")
        print ("+ ---- + --- + --- + --- + ------ + --- + --- + ----+")
        print ("| ATOM | X   |  Y  |  Z  | CHARGE | VX  |  VY |  VZ |")  
        print ("| ---- | --- | --- | --- | ------ | --- | --- | ----|") 
        print ("| N    | 0   |  0  |  0  |   5.0  | 0.0 | 0.1 |  0.1|")
        print ("+ ---- + --- + --- + --- + ------ + --- + --- + ----+")
        print ("========================================") 
        sys.exit(1)
        
    else:
        coors, chgs, vels = parse_xyz(dataSourceFileName)
       
        nstp, nat, ncoor = vels.shape
        vels = vels.reshape(nstp, nat*ncoor)

        acf = np.apply_along_axis(acf, 0, vels)
        acf = acf.reshape(nstp, nat, ncoor)
        acf = acf.sum(axis=2)
        acf = acf.sum(axis=1)
        x = np.arange(nstp)

#        Plotting Spectrogram!
        ts = 1e-10
        f, t, Sxx = spectrogram(acf, ts)
        plt.pcolormesh(t, f, Sxx)

#       Plotting fast fourier transform 
#        ft = np.fft.fft(acf)
#        freq = np.fft.fftfreq(nstp)
#        plt.plot(freq[:int(nstp/2.0)], np.abs(ft[:int(nstp/2.0)]**2))
        plt.show()

        # ++++ TESTING PURPOSES: See current dict content, this is how we access it ++++
        #print2DimDict("VelocitiesDict", velocitiesDict)
        #print2DimDict("Autocorrelation dict", autocorrelationDict)
        
    
