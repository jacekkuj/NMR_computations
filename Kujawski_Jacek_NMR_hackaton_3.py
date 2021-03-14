#!/usr/bin/env python
"""
Program ma odczytywać sygnały z pliku wynikowego *.log (obliczenia z użyciem Gaussian G16A).

Odpalamy np. z konsoli np.:

Kujawski_Jacek_NMR_hackathon_3.py paracetamol_DMSO_NMR.log paracetamol_DMSO_NMR.txt

Jako dane eksperymentalne proponuję wpdać (np. dla widma 1H NMR paracetamolu):
7.350
7.350
6.688
6.688
9.66
1.987
1.987
1.987
9.14
"""

import os, sys
import numpy as np
import sklearn
import sklearn.metrics as metrics
from sklearn.metrics import mean_absolute_error as mae
import matplotlib.pyplot as plt

class NMR:

	def starting_program(self):
		if len(sys.argv) < 3:
			sys.exit('USAGE: %s [input_file] [output_file]\n' % sys.argv[0])
		self.finput = open(sys.argv[1], 'r', encoding='utf-8')
#		self.foutput = open(sys.argv[2], 'w+', encoding='utf-8')
		self.lines = self.finput.readlines()
		self.lines = [self.lines[i].strip() for i in range(len(self.lines))]
		self.atom_numbers = []
		self.extracted_peaks = []
		print("Extracting Isotropic values from a *.log file...")
		for i in range(len(self.lines)):
			if str.find(self.lines[i], 'H    Isotropic =') != -1:
#				self.foutput.write(str(self.lines[i].split()[0]) + '\t' + str(self.lines[i].split()[4]) + '\n')
				self.atom_numbers.append(int(self.lines[i].split()[0]))
				self.extracted_peaks.append(float(self.lines[i].split()[4]))
		print("Extraction of Isotropic values done and written!")


	def ref_TMS(self):
		self.href = float(input("Enter 'H_ref' value in ppm (i.e. 31.6674): "))

	def empirical_peaks(self):
		print(u"Entering experimental signals from 1H NMR spectrum.\n")
		self.empiricalPeaks = []
		peak = 0
		for peak in range(len(self.atom_numbers)):
			print(f"Atom number {self.atom_numbers[peak]}. \nEnter experimental shift [ppm]: ")
			self.empiricalPeak = input()
			try:
				self.empiricalPeaks.append(float(self.empiricalPeak))
				peak += 1
			except ValueError:
				print(u'\nValue error!\n')
				continue

	def computed_peaks(self):
		self.computedPeaks = []
		for peak in self.extracted_peaks:
			self.computedPeaks.append(self.href - peak)

	def writing_results(self):
		self.result_file_name = f"{sys.argv[2][0:20]}_results.txt"
		self.result_file = open(self.result_file_name, 'w+', encoding='utf-8')
		print("\nResults given in [ppm]:\n")
		header = "Hydrogen\t%s\t%s\t%10s\t%s\n" % (u'Theoretical', '\tExperimental', '\tError', '\tRelative error')
		print(header)
		for i in range(len(self.computedPeaks)):
			print(u"%dH\t%19.4f\t%23.4f\t%10.4f\t%13.4f" % (self.atom_numbers[i], self.computedPeaks[i], self.empiricalPeaks[i], self.empiricalPeaks[i] - self.computedPeaks[i],
			(self.empiricalPeaks[i] - self.computedPeaks[i]) / self.computedPeaks[i]))
			self.result_file.write(u"%dH\t%19.4f\t%23.4f\t%10.4f\t%13.4f\n" % (self.atom_numbers[i], self.computedPeaks[i], self.empiricalPeaks[i], self.empiricalPeaks[i] - self.computedPeaks[i],
			(self.empiricalPeaks[i] - self.computedPeaks[i]) / self.computedPeaks[i]))
		from sklearn.metrics import mean_absolute_error as mae
		MAE = mae(self.empiricalPeaks, self.computedPeaks)
		print(f"MAE: {MAE} ppm")
		self.result_file.close()

# Do poprawienia jest wygenerowanie wykresu analizy regresji,
# tj. analiza statystyczna dla dwóch list:
# self.computedPeaks oraz self.empiricalPeaks
	def charts(self):
		self.computedPeaks = np.array(self.computedPeaks)
		self.empiricalPeaks = np.array(self.empiricalPeaks)
		peaks = np.array([self.computedPeaks, self.empiricalPeaks]).flatten()
		peaks2 = np.array([self.empiricalPeaks, self.computedPeaks]).flatten()
		peaksCorrCoef = np.corrcoef(peaks, peaks2)[0, 1]
		chart = plt.figure(spectrum[:-4])
		empiricalScatter = plt.scatter(self.computedPeaks, self.empiricalPeaks, marker='x', color='b')
		computedScatter = plt.scatter(self.computedPeaks, self.computedPeaks, marker='o', color='g')
		plt.legend((empiricalScatter, computedScatter),
				   ('Experimental', 'Theoretical'),
				   scatterpoints=1,
				   loc='lower right',
				   ncol=3,
				   fontsize=12)
		plt.xlabel('Experimental [ppm]')
		plt.ylabel('Theoretical [ppm]')
		polyCoeff = np.polyfit(peaks, peaks2, 1)
		fittedPoly = (self.computedPeaks[:, 1] * polyCoeff[0]) + polyCoeff[1]
		plt.plot(self.computedPeaks[:, 1], fittedPoly)
		plt.figtext(0.25, 0.75, 'r^2= %.5f' % (peaksCorrCoef ** 2), fontsize=12)
		plt.figtext(0.25, 0.8, 'y = %.3fx + %.3f' % (polyCoeff[0], polyCoeff[1]))
		plt.title("NMR correlation")
		plt.show()
		chart.savefig(spectrum[:-3] + 'png')

if __name__ == '__main__':
	nmr = NMR()
	nmr.starting_program()
	nmr.ref_TMS()
	nmr.empirical_peaks()
	nmr.computed_peaks()
	nmr.writing_results()
	nmr.charts()


