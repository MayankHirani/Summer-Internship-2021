
# Libraries
import functions
from importlib import reload
reload(functions)
from functions import *
import os
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression


d = {'Test':['A', 'B'], 'Another':['C', 'D'], 'Third':['E', 'F']}
df = pd.DataFrame(data=d)
# df.to_csv('test.csv', sep='\t')
print(df.columns[1])
print(list(df.loc[:, 'Test']))
print(df)


# Main class to read in all information
class Main():

	def __init__(self):
		# Dictionaries store key: chromosome and section, value: the DNA sequence at that part
		self.test_sequences = {}; self.train_sequences = {}
		# All the PWMs stored in the file, key: filename, value: the PWM
		self.pwm_data = {}
		# Data for each TF. Key: Name of PWM for binding, Value: A dictionary of each sequence and its high score
		self.score_data = {}; self.score_data_test = {}

	# Retrieve the sequences and format them as one-hot matrices
	def retrieve(self, tf_name):

		# Open the two files for test data and train data
		f = open('data/topic1-data/' + tf_name + '/' + tf_name + '-test-sequence.fa', 'r').readlines()
		g = open('data/topic1-data/' + tf_name + '/' + tf_name + '-train-sequence.fa', 'r').readlines()

		# print(len(f)); print(len(g))

		# Make a dictionary with keys being the sequences name and values the sequence
		for x in range(0, len(f), 2):
			self.test_sequences[f[x].strip('\n').strip('>')] = functions.onehot(f[x+1].strip('\n').upper())

		for x in range(0, len(g), 2):
			self.train_sequences[g[x].strip('\n').strip('>')] = functions.onehot(g[x+1].strip('\n').upper())

		# print(self.test_sequences)
		# print(self.train_sequences)

	# Retrieve the PWMs of the TF
	def getPWM(self, tf_name):

		# Go through the TF's folder and read through the WTMX files
		directory = 'data/topic1-data/' + tf_name
		for file in os.listdir(directory):
			if file.endswith('.wtmx'):
				
				# Read the file
				f = open(os.path.join(directory, file), 'r').readlines()

				# Add the PWM to the dictionary with the key as the file title
				pwm = functions.pwm_from_wtmx(f)
				pwm = functions.changeValues(pwm)
				self.pwm_data[file.split('.')[0]] = pwm

		# printing
		# for key in self.pwm_data.keys():
			# print(key)
			# print(self.pwm_data[key])
			# print()

	# Get the highest score of each sequence with each PWM
	def getSequenceScores(self):

		# Go through each PWM
		for key in self.pwm_data.keys():

			# -- TRAIN --
			pwm = self.pwm_data[key]

			# A dictionary of each sequence and its high score
			pwmScores = {}

			# Go through each sequence as one hot matrices
			for sequenceKey in self.train_sequences.keys():
				
				# Add the score and sequence name to the dictionary
				sequence = self.train_sequences[sequenceKey]
				pwmScores[sequenceKey] = functions.calculateBestScore(sequence, pwm)

			# Add it to the full dictionary with all the PWMs
			self.score_data[key] = pwmScores


			# -- TEST --
			pwm = self.pwm_data[key]

			# A dictionary of each sequence and its high score
			pwmScores = {}

			# Go through each sequence as one hot matrices
			for sequenceKey in self.test_sequences.keys():
				
				# Add the score and sequence name to the dictionary
				sequence = self.test_sequences[sequenceKey]
				pwmScores[sequenceKey] = functions.calculateBestScore(sequence, pwm)

			# Add it to the full dictionary with all the PWMs
			self.score_data_test[key] = pwmScores

		# print(self.score_data)
		# print(self.score_data_test)

	# Used to make a dataframe of the data (col1: Sequences, col2: ChIP values in .bed file, col3+: All the PWMs and their high scores one by one)
	def makeDataframe(self, tf_name):

		# At this point, score_data is complete, so we can pull our values directly from there.
		data = {'Sequences':[], 'ChIP Score':[]}; data_test = {'Sequences':[], 'ChIP Score':[]}

		# [TRAIN]
		# Add columns for the PWMs for the binding sites
		for pwm_name in self.pwm_data.keys():
			data[pwm_name] = []

		# Add the sequences and ChIP scores
		chip_scores = open('data/topic1-data/' + tf_name + '/' + tf_name + '-train.bed', 'r')
		for line in chip_scores.readlines():
			# Add the sequence names
			data['Sequences'].append(line.split('\t')[0])
			# Add the ChIP score
			data['ChIP Score'].append(line.split('\t')[1].strip('\n'));

		# Go through each sequence and add the appropriate high score for each PWM
		for sequence in data['Sequences']:
			# Iterate through each PWM
			for pwm_name in self.pwm_data.keys():
				# Add the high score of that PWM with that sequence
				data[pwm_name].append(self.score_data[pwm_name][sequence])

		# [TEST]
		# Add columns for the PWMs for the binding sites
		for pwm_name in self.pwm_data.keys():
			data_test[pwm_name] = []

		# Add the sequences and ChIP scores
		chip_scores = open('data/topic1-data/' + tf_name + '/' + tf_name + '-test.bed', 'r')
		for line in chip_scores.readlines():
			# Add the sequence names
			data_test['Sequences'].append(line.split('\t')[0])
			# Add the ChIP score
			data_test['ChIP Score'].append(line.split('\t')[1].strip('\n'));

		# Go through each sequence and add the appropriate high score for each PWM
		for sequence in data_test['Sequences']:
			# Iterate through each PWM
			for pwm_name in self.pwm_data.keys():
				# Add the high score of that PWM with that sequence
				data_test[pwm_name].append(self.score_data_test[pwm_name][sequence])


		# Create the dataframes and print them
		self.df_train = pd.DataFrame(data=data)
		self.df_test = pd.DataFrame(data=data_test)
		self.df_train.to_csv('train_data.csv', sep='\t')
		self.df_test.to_csv('test_data.csv', sep='\t')
		print(self.df_train)
		print('\n\n\n\n')
		print(self.df_test)

	# Make the linear regression model and analyze the scores
	def linearRegression(self):

		# Data for the dataframe. First column is the PWM name, second is the train score, third is the test score
		data = {'PWM':[], 'Regression Score (Train)':[], 'Regression Score (Test)':[]}

		# List of ChIP scores from the train sequences
		chip_scores_train = self.df_train.loc[:, 'ChIP Score']
		chip_scores_test = self.df_test.loc[:, 'ChIP Score']

		# For each PWM, get a linear regression model
		for PWM in self.df_train.columns[2:]:
			# The PWM scores
			pwm_scores_train = np.array(self.df_train.loc[:, PWM]).reshape(-1, 1)
			pwm_scores_test = np.array(self.df_test.loc[:, PWM]).reshape(-1, 1)
			# Create the Linear Regression model for the train scores
			model = LinearRegression().fit(pwm_scores_train, chip_scores_train)
			training_score = model.score(pwm_scores_train, chip_scores_train)
			testing_score = model.score(pwm_scores_test, chip_scores_test)
			# Add our values to the dataframe data
			data['PWM'].append(PWM)
			data['Regression Score (Train)'].append(training_score)
			data['Regression Score (Test)'].append(testing_score)

			print(PWM, training_score, testing_score)

		# Make the dataframe with the scores
		self.df_scores = pd.DataFrame(data=data)
		self.df_scores.to_csv('final_scores.csv', sep='\t')
		print(self.df_scores)


# TF name
tf_name = input('TF Name: ').upper()

# Run
# Functions.retrieve(tf_name)
reader = Main()
reader.retrieve(tf_name)
reader.getPWM(tf_name)
reader.getSequenceScores()
reader.makeDataframe(tf_name)
reader.linearRegression()
