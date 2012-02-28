'''Modul obsluguje sygnal zapisany w formacie SignalML (http://signalml.org/) wraz z plikiem xml opisujacym jego zawartosc. Poszczegolne metody zwracaja dane pozyskane z pliku xmlowego w formacie umozliwiajacym ich wykozystanie podczas analizy sygnalu. Mozliwe jest tez pozyskanie za pomoca klasy poszczegolnych kanalow zapisanych w pliku raw, jak rowniez np. numerow probek rozpoczynajacych kazde wystapienie bodzca (trigger). 


'''
import numpy
import pylab
from lxml import etree

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
    def channel(self, number):
#Return signal from specified channel
        if int(number)>= int(self.cc):
            print 'Error: Channels available: 0 to ' + str(int(self.cc) - 1)
            return None
        return self.s[int(number)::int(self.cc)]  
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
        trigg = self.channel(z)
        triggers = []
        c = []
        for i in range(0, len(trigg)):
            if trigg[i] == 1:
                c.append(i)
        for k in c:
            if trigg[k-1] != 1:
                triggers.append(k)
        return triggers


projekt = SignalML('ania-michalska-trigger1.raw', 'ania-michalska-trigger1.xml')
print projekt.channel(2)
