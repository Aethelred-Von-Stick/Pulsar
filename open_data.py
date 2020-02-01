import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
import seaborn
from numpy.fft import fft, fftfreq

#this cell takes text files containing the julian date, $\omega_0$, time interval and complex data and creates a
#dictionary of these files.

file_extention = '.CPX'
#the number of lines that are not complex data
number_of_header_lines = 3


path_to_directory = '/media/sam/Windows/Users/samwa/Anaconda3/envs/irregularsatellites/myscripts/Pulsars/data/data/'
#make the directory readable to OS
directory = os.fsencode(path_to_directory)
#create a list of the files in the directory
directory_list = os.listdir(directory)
#initiate a list in which to put the data from each file

data_files = []
#iterate over all files in the directory
for i in range(len(directory_list)):
    #make the filename human readable
    filename = os.fsdecode(directory_list[i])
    #choose only the text files with data in it
    if filename.endswith(file_extention):
        #open and read the file
        file_path = path_to_directory + filename
        f = open(file_path, 'r')
        #create an array of the complex data in the form [|real|, |imaginary|, |real|...]
        real_data = np.genfromtxt(file_path, skip_header = number_of_header_lines, usecols = 0)
        im_data = np.genfromtxt(file_path, skip_header = number_of_header_lines, usecols = 1)
        #slice the data into two arrays, one of the real part and one of the magnitude of the imaginary part
        #combine the real and imaginary parts into one complex array with one index per complex number
        complex_array = real_data + 1j*im_data
        
        #create dictionary keys for each file
        #the name of the key is the same as the filename without the extention
        jd = Time(float(f.readline()), format = 'jd')
        Gregorian_date = jd.isot
        w0 = float(f.readline())
        time_interval = float(f.readline())
        end_time = len(complex_array)*time_interval
        file_dict = {
            "jd" : jd,
            "w0" : w0,
            "UTC_date" : Gregorian_date,
            "time_interval" : time_interval,
            "end_time" : end_time,
            "t" : np.arange(0,end_time,time_interval),
            "complex_data" : complex_array
        }
        f.close()
        data_files.append(file_dict)

#plot all the timeseries in a subplot 
nrows = 4
ncols = 2
#initiate a subplot
fig, ax = plt.subplots(nrows = nrows, ncols = ncols, figsize = (10,10))
#loop over the subplot positions and label the x and y axes
for a in range(nrows):
    for b in range(ncols):
        ax[a][b].set_xlabel('time (s)')
        ax[a][b].set_ylabel('amplitude')
        ax[a][b].ticklabel_format(style = 'sci', scilimits = (0,0))
        #k is the plot number of a flattened subplot
        k = a*ncols + b
        ax[a][b].plot(data_files[k]["t"], data_files[k]["complex_data"])

plt.show()

for m in range(len(data_files)):
    f = data_files[m]
    data_files[m]['G'] = fft(abs(f['complex_data'])**2)
    data_files[m]['n_samples'] = len(f['complex_data'])
    data_files[m]['freq'] = fftfreq(data_files[m]['n_samples'], d = f['time_interval'])

fig, ax = plt.subplots(nrows = nrows, ncols = ncols, figsize = (10,10))
#loop over the subplot positions and label the x and y axes
for c in range(nrows):
    for d in range(ncols):
        ax[c][d].set_xlabel('frequency (Hz)')
        ax[c][d].set_ylabel('power')
        ax[c][d].ticklabel_format(style = 'sci', scilimits = (0,0))
        #h is the plot number of a flattened subplot
        h = c*ncols + j
        ax[c][d].plot(data_files[h]["freq"], data_files[h]["G"])

plt.show()