import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
from numpy.fft import fft, fftfreq

# this cell takes text files containing the julian date, $\omega_0$, time interval and complex data and creates a
# dictionary of these files.

file_extention = '.CPX'
# the number of lines that are not complex data
number_of_header_lines = 3

path_to_directory = os.getcwd() + '/data/data/'
# make the directory readable to OS
directory = os.fsencode(path_to_directory)
# create a list of the files in the directory
directory_list = os.listdir(directory)
# initiate a list in which to put the data from each file

data_files = []
# iterate over all files in the directory
for i in range(len(directory_list)):
    # make the filename human readable
    filename = os.fsdecode(directory_list[i])
    # choose only the text files with data in it
    if filename.endswith(file_extention):
        # open and read the file
        file_path = path_to_directory + filename
        f = open(file_path, 'r')
        # create an array of the complex data in the form [|real|, |imaginary|, |real|...]
        real_data = np.genfromtxt(file_path, skip_header=number_of_header_lines, usecols=0)
        im_data = np.genfromtxt(file_path, skip_header=number_of_header_lines, usecols=1)
        # slice the data into two arrays, one of the real part and one of the magnitude of the imaginary part
        # combine the real and imaginary parts into one complex array with one index per complex number
        complex_array = real_data + 1j * im_data
        pads = np.zeros(100000)
        complex_array = np.hstack((complex_array, pads))

        # create dictionary keys for each file
        # the name of the key is the same as the filename without the extention
        jd = Time(float(f.readline()), format='jd')
        Gregorian_date = jd.isot
        w0 = float(f.readline())
        time_interval = float(f.readline())
        end_time = len(complex_array) * time_interval
        file_dict = {
            "jd": jd,
            "w0": w0,
            "UTC_date": Gregorian_date,
            "time_interval": time_interval,
            "end_time": end_time,
            "t": np.arange(0, end_time, time_interval),
            "complex_data": complex_array
        }
        f.close()
        data_files.append(file_dict)

## plot all the timeseries in a subplot
nrows = 3
ncols = 3
# initiate a subplot
fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10, 10))
# loop over the subplot positions and label the x and y axes
for k, ax in enumerate(fig.axes):
        ax.set_xlabel('time (s)')
        ax.set_ylabel('amplitude')
        ax.ticklabel_format(style='sci', scilimits=(0, 0))
        # k is the plot number of a flattened subplot
        ax.plot(data_files[k]["t"], data_files[k]["complex_data"])

for m in range(len(data_files)):
    f = data_files[m]
    data_files[m]['G'] = abs(fft(f['complex_data']))**2
    data_files[m]['n_samples'] = len(f['complex_data'])
    data_files[m]['freq'] = fftfreq(data_files[m]['n_samples'], d=f['time_interval'])

fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10, 10))
# loop over the subplot positions and label the x and y axes
for c in range(nrows):
    for d in range(ncols):
        ax[c][d].set_xlabel('frequency (Hz)')
        ax[c][d].set_ylabel('power')
        ax[c][d].ticklabel_format(style='sci', scilimits=(0, 0))
        # h is the plot number of a flattened subplot
        h = c * ncols + d
        ax[c][d].set_xlim(left=-0.005, right=0.005)
        ax[c][d].plot(data_files[h]["freq"], data_files[h]["G"], '.')


for r in range(len(data_files)):
    zipped = zip(data_files[r]['G'],data_files[r]['freq'])
    zipped = sorted(zipped, key=lambda x: x[0])
    data_files[r]['peak_freq'] = zipped[-1][1]
    power = [power[0] for power in zipped]
    frequency = [frequency[1] for frequency in zipped]
    sfreq = sorted(frequency)
    for i in range(len(sfreq)-1):
        if sfreq[i] - sfreq [i-1] == sfreq[i+1] - sfreq[i]:
            same_interval = True
        else:
            same_interval = False        
    if same_interval == True:
        data_files[r]['peak_err'] = sfreq[1]-sfreq[0]
    else:
        print('the frequency interval is not constant')

peaks = []
err_peaks = []
timesjd = []
for p in range(len(data_files)):
    peaks.append(data_files[p]['peak_freq'])
    err_peaks.append(data_files[p]['peak_err'])
    timesjd.append(data_files[p]['jd'].jd)

fi, a = plt.subplots()
a.errorbar(timesjd, peaks, yerr = err_peaks, )

plt.show()
