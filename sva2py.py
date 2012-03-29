# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Author:
#     Karol Augustin <karol@augustin.pl>
# git repository: http://git.nimitz.pl

import numpy
import pylab
from lxml import etree
from scipy.signal import filtfilt, butter, freqz, lfilter, cheb1ord, cheby2
import scipy




class sva2py:
#Mozna podac jeden parametr jako istotny czlon nazwy plikow, rozszezenia zostana dodane automatycznie, lub dwie nazwy z rozszezeniami oddzielnie dla pliku raw i xml.
    def __init__(self, file_name, xml_file_name = False):
        if xml_file_name:
            self.s = numpy.fromfile(str(file_name),'double')
            tree = etree.parse(str(xml_file_name))
        else:
            self.s = numpy.fromfile(str(file_name)+'.raw','double')
            tree = etree.parse(str(file_name)+'.xml')
        ns = {'rs': 'http://signalml.org/rawsignal'}
        self.fs = tree.xpath('rs:samplingFrequency', namespaces = ns)[0].text
        self.cc = tree.xpath('rs:channelCount', namespaces = ns)[0].text
        self.sc = tree.xpath('rs:sampleCount', namespaces = ns)[0].text
        self.st = tree.xpath('rs:sampleType', namespaces = ns)[0].text
        self.bo = tree.xpath('rs:byteOrder', namespaces = ns)[0].text
        self.fsts = tree.xpath('rs:firstSampleTimestamp', namespaces = ns)[0].text
        self.cl = tree.xpath('rs:channelLabels/rs:label', namespaces = ns)
        self.file_name = file_name
    def __str__(self):
        intro =  'Svarog signal import: '+ str(self.file_name)
        cc = 'Number of channels: ' + str(self.cc)
        sf = 'Sampling frequency: ' + str(self.fs)
        sc = 'Sample count: ' + str(self.sc)
        st = 'Sample type: ' + str(self.st)
        bo = 'Byte order: ' + str(self.bo)
        fsts = 'First sample timestamp ' + str(self.fsts)
        cl = 'Channel labels: \n'
        v = 0
        for i in range(len(self.cl)):
            cl = cl + '   Channel '+ str(v)+ ': ' + self.cl[i].text + '\n'
            v = v + 1
        return intro + '\n' + cc + '\n' + sf + '\n' + sc + '\n' + st + '\n' + bo + '\n' + fsts + '\n' + cl + '\n'
    def channel(self, number, type='name'):
#Return signal from specified channel
        if type == 'int' or type == None:
            number = int(number)    
            if number >= int(self.cc):
                raise ValueError     #'Error: Channels available: 0 to ' + str(int(self.cc) - 1)            
            return self.s[int(number)::int(self.cc)]
        elif type == 'str' or type == 'name':
            i = 0
            while self.cl[i].text != number:
                i += 1
            return self.s[i::int(self.cc)] 
        else:
            return None # raise something
    
    def montage(self, name, type='linkedears', filtr='high', start=None, stop=None, mixed=False):
        if filtr:
            [b,a] = butter(3,1.0/(self.samplingFrequency()/2.0), btype='high')
        if filtr == 'alpha':
            Wn = [8.5/(128/2.0),14.5/(128/2.0)]
            [g,h] = cheby2(4, 20, Wn, btype='bandpass', analog=0, output='ba')
        if type == 'linkedears':
            s = self.channel(name, 'name')
            #s = filtfilt(b,a,s - (self.channel('A1', 'name') + self.channel('A2', 'name'))/2)
            if start and stop: 
                t = filtfilt(b,a,s[start:stop] - (self.channel('A1', 'name')[start:stop] + self.channel('A2', 'name')[start:stop])/2)
                if filtr == 'alpha': t = filtfilt(g,h,t)
            elif stop:
                t = filtfilt(b,a,s[:stop] - (self.channel('A1', 'name')[:stop] + self.channel('A2', 'name')[:stop])/2)
                if filtr == 'alpha': t = filtfilt(g,h,t)
            elif start:
                t = filtfilt(b,a,s[start:] - (self.channel('A1', 'name')[start:] + self.channel('A2', 'name')[start:])/2)
            else:
                t = filtfilt(b,a,s - (self.channel('A1', 'name') + self.channel('A2', 'name'))/2)
                if filtr == 'alpha': t = filtfilt(g,h,t)
            if mixed: t = self.mixer(t)
        return t

    def autocorrelate(self, c, montage='linkedears', filtr='high', start=None, stop=None, norm=True, mixed=False):
        if montage == 'direct':
            signal = c
        else:
            signal = self.montage(c, type=montage, filtr=filtr, start = start, stop = stop, mixed=mixed)
        correlation = numpy.correlate(signal, signal, 'full')
        if norm:
            return correlation/max(correlation)
        else:
            return correlation

    def correlate(self, c1, c2, montage='linkedears', filtr='high', start=None, stop=None, norm=True, mixed=False):
        if montage == 'direct':
            signal1 = c1
            signal2 = c2
        else:
            signal1 = self.montage(c1, type=montage, filtr=filtr, start=start, stop=stop, mixed=mixed)
            signal2 = self.montage(c2, type=montage, filtr=filtr, start=start, stop=stop, mixed=mixed)
        correlation = numpy.correlate(signal1, signal2, 'full')
        if norm: correlation = correlation/(numpy.std(signal1)*numpy.std(signal2))/len(signal1)
        return correlation
    
    def periodogram(self, c, filtr='high', start=None, stop=None, montage='linkedears', window='blackman', use=True):
        if use:
            s = self.montage(c, filtr=filtr, start=start, stop=stop, type=montage)
            okno = self.window(window, len(s))
        else:
            okno = window
            s = c
        F_samp = self.samplingFrequency()
        s = s*okno
        N_fft = len(s)
        S = numpy.fft.rfft(s,N_fft)#/numpy.sqrt(N_fft)
        P = S*S.conj()/numpy.sum(okno**2)
    
        P = P.real 
        F = numpy.linspace(0, F_samp/2., len(S))
        return (numpy.fft.fftshift(P),numpy.fft.fftshift(F))


    def fft(self, signal):
        s = numpy.fft.rfft(signal)
        f = numpy.linspace(0, self.samplingFrequency()/2., len(s))
        return numpy.abs(s), f

    def window(self, name, N):
        if name == 'blackman':
            window = numpy.blackman(N)
        elif name == 'hamming':
            window = numpy.hamming(N)
        return window

    def pwelch(self, c, window='blackman', filtr='high', przesuniecie=1/10., dlugosc=1/8., start=None, stop=None, montage='linkedears'):
        Fs = self.samplingFrequency()
        s = self.montage(c, filtr=filtr, type=montage, start=start, stop=stop)
        N = len(s)
        przesuniecie = przesuniecie*dlugosc*N
        okienko = self.window(window, dlugosc*N)
        N_s = len(okienko)
        start_fragmentow = numpy.arange(0,N-N_s+1,przesuniecie)
        ile_fragmentow = len(start_fragmentow)
        print 'Fragmentow: ', ile_fragmentow
        ile_przekrycia = N_s*ile_fragmentow/float(N)
        P_sredni = numpy.zeros(N_s/2+1)
        for i in range(ile_fragmentow):
            s_fragment = s[start_fragmentow[i]:start_fragmentow[i]+N_s]
            (P, F) = self.periodogram(s_fragment, window=okienko, use=False)
            #print numpy.shape(P), numpy.shape(P_sredni)
            P_sredni += P
        return P_sredni/ile_przekrycia, F

    def mixer(self, x):
        y = scipy.signal.irfft(numpy.abs(numpy.fft.rfft(x))*scipy.signal.exp(1j*numpy.random.random(len(x)/2+1)*2*numpy.pi))
        #y = numpy.real(scipy.signal.ifft(numpy.abs(scipy.signal.fft(x))*scipy.signal.exp(1j*numpy.random.random(len(x))*2*numpy.pi)))
        return y

    def signal(self):
#Return whole signal
        return self.s
    def samplingFrequency(self):
#Return sampling frequency
        return float(self.fs)
    def channelCount(self):
#Return channel count
        return float(self.cc)
    def firstSampleTimestamp(self):
#Return first sample timestamp
        return float(self.fsts)
    def sampleCount(self):
#Return sample count
        return float(self.sc)
    def plot(self, channel, filename = False, dpi = 400,  xaxis='sample', yaxis='voltage'):
#Plot signal from specified channel (to file)
        x = self.channel(channel)
        pylab.plot(x)
        if filename:
            pylab.savefig(filename, dpi = dpi)
        else:
            pylab.show()
        return None
    def trigger(self):
#Return first samples for each trigger press
        for i in range(0, len(self.cl)):
            if self.cl[i].text == 'TRIGGER':
                z = i
                break
        trigg = self.channel(z, type='int')
        triggers = []
        c = []
        for i in range(0, len(trigg)):
            if trigg[i] == 1:
                c.append(i)
        for k in c:
            if trigg[k-1] != 1:
                triggers.append(k)
        return triggers

    def trigger2(self, c, type='name'):
        z =[]
        q = []
        x = self.channel(c, type=type)
        x = x - min(x)
        #a = (max(x) - min(x))/2.
        a = numpy.average(x)
        for i in range(1,len(x)):
            if x[i] > a and x[i-1]<a:
                z.append(i)
#        for k in z:
#            if x[k-1] < a:
#                q.append(k)
        return z #q
    def trigger3(self, c1, c2, type='name'):
        r1 = []
        r2 = []
        t1 = self.trigger2(c1, type=type)
        t2 = self.trigger2(c2, type=type)
        if len(t1)>len(t2):
            t1b = t1
            t2b = t2
        else:
            t2b = t1
            t1b = t2

        for i in t1b:
            if t2b.count(i) == 1:
                r2.append(i)
            else:
                r1.append(i)
        return r1, r2

    def get_p300_tags(self, idx=1, rest=False, samples = True, Fs = None):
        from xml.dom import minidom

        """Returns tags with words from different groups
        
        Parameters:
        -----------
        idx [= 1]: int
            defines which tags to return
        samples : bool
            if true, positions will be returned as samples not in seconds
        Fs : float or None
            the sampling frequency used to convert positions to samples
        
        Returns:
        --------
        exp_list : list
            a list of positions of target
        """

        ftag = minidom.parse(str(self.file_name)+'.tag')
        tag_list = [e for e in ftag.getElementsByTagName('tag') \
                    if e.attributes['name'].value == 'blink']
        exp_list = []

        fsp = self.samplingFrequency()
        if(samples):
            if Fs != None:
                fsp = Fs
        else: fsp = 1.0
        for e in tag_list:
            index = e.getElementsByTagName('index')[0].firstChild.data
            if not rest:
                if int(index) == idx:
                    exp_list.append(float(e.attributes['position'].value))
            else:
                if int(index) !=  idx:
                    exp_list.append(float(e.attributes['position'].value))
        
        return numpy.array(exp_list) * fsp 
